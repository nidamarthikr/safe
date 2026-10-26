#!/usr/bin/env python3
"""
SAFE - Statistical Analysis of Structure Factors and Normalised E Values

run_instat.py

Prepare and execute INSTAT calculations from native STIN files or
supported crystallographic structure files.

    Author
    ------
    Kutumbarao Nidamarthi
    Beamline Scientist, MX Beamline
    Diamond Light Source, UK
    Email: kutumbarao.nidamarthi@diamond.ac.uk

Compile and run INSTAT in either of two modes.

Mode 1, existing INSTAT control file plus CIF-formatted reflection file:
    safe-run-instat --stin data/STIN --hkl data/hb5427-sup2.hkl

Mode 2, generate STIN from a structure CIF, SHELX INS, or SHELX RES file,
plus its reflection file:
    safe-run-instat --structure data/model.res --hkl data/model.hkl
    safe-run-instat --structure data/model.cif --hkl data/model.hkl

The bundled INSTAT Fortran source is located automatically
within the SAFE installation.
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from collections import OrderedDict
from pathlib import Path
from datetime import datetime

SCRIPT_DIR = Path(__file__).resolve().parent
SOURCE_FILE = SCRIPT_DIR / "src" / "instat.f"


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--stin",type=Path,
        help="Existing INSTAT STIN file. Composition and controls are read from it.",)
    mode.add_argument("--structure",type=Path,
        help="Structure CIF, SHELX INS, or SHELX RES file used to generate STIN.",)
    parser.add_argument("--hkl",type=Path,
        required=True,help="CIF-formatted reflection HKL or raw SHELX HKL file.",)
    parser.add_argument("--workdir",type=Path,
        default=None,help="Optional output directory. Default: results/<input>_<timestamp>",)
    parser.add_argument("--compiler",
            default="/usr/bin/gfortran",help="Fortran compiler. Default: /usr/bin/gfortran",)
    parser.add_argument("--composition",
        help="Override ASU composition, for example C:24,H:45,N:1,O:1.",)
    parser.add_argument("--nunit",type=int,
        help="Override number of asymmetric units in the unit cell.",)
    parser.add_argument("--title",
        help="Override the generated STIN title.",)
    parser.add_argument("--ispg", type=int, default=1)
    parser.add_argument("--iser1", type=int, default=1)
    parser.add_argument("--iser2", type=int, default=2)
    parser.add_argument("--jsym", type=int, default=9)
    parser.add_argument("--ncoef", type=int, default=40)
    parser.add_argument("--nbesk", type=int, default=20)
    parser.add_argument("--ilatt", type=int, default=1)
    return parser.parse_args()


def require_file(path: Path, label: str) -> Path:
    path = path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(f"{label} not found: {path}")
    return path


def strip_cif_value(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in "'\"":
        value = value[1:-1]
    return re.sub(r"\([^)]*\)$", "", value).strip()


def cif_scalar(text: str, tags: list[str]) -> str | None:
    for tag in tags:
        match = re.search(
            rf"(?im)^\s*{re.escape(tag)}\s+(.+?)\s*$",
            text,)
        if match:
            return strip_cif_value(match.group(1))
    return None


def parse_formula(formula: str) -> OrderedDict[str, int]:
    """Parse common crystallographic formula strings such as C24 H45 Au N O P S.
    Decimal, disordered, parenthesized, or polymeric formulae should be supplied
    explicitly with --composition.
    """
    formula = strip_cif_value(formula).replace("·", " ").replace(".", " ")
    tokens = re.findall(r"([A-Z][a-z]?)(?:\s*([0-9]+(?:\.[0-9]+)?))?", formula)
    if not tokens:
        raise ValueError(f"Could not parse chemical formula: {formula!r}")
    composition: OrderedDict[str, int] = OrderedDict()
    for element, count_text in tokens:
        count_float = float(count_text) if count_text else 1.0
        if not count_float.is_integer():
            raise ValueError(
                "Non-integer formula count found. Use --composition with integer ASU counts.")
        composition[element.upper()] = composition.get(element.upper(), 0) + int(count_float)
    return composition


def parse_composition_override(value: str) -> OrderedDict[str, int]:
    composition: OrderedDict[str, int] = OrderedDict()
    for item in value.split(","):
        try:
            element, count = item.split(":", 1)
            element = element.strip().upper()
            number = int(count.strip())
        except ValueError as exc:
            raise ValueError(
                "Composition must look like C:24,H:45,N:1,O:1"
            ) from exc
        if not re.fullmatch(r"[A-Z][A-Z]?", element) or number < 0:
            raise ValueError(f"Invalid composition item: {item!r}")
        if number:
            composition[element] = number
    if not composition:
        raise ValueError("Composition override is empty.")
    return composition


def parse_structure_cif(path: Path) -> tuple[str, list[float], OrderedDict[str, int], int]:
    text = path.read_text(encoding="utf-8", errors="replace")
    title = cif_scalar(text, ["_publ_section_title", "_chemical_name_common", "_shelx_title"])
    if not title:
        data_match = re.search(r"(?im)^\s*data_(\S+)", text)
        title = data_match.group(1) if data_match else path.stem
    cell_tags = ["_cell_length_a", "_cell_length_b", "_cell_length_c",
        "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma",]
    values = [cif_scalar(text, [tag]) for tag in cell_tags]
    if any(value is None for value in values):
        missing = [tag for tag, value in zip(cell_tags, values) if value is None]
        raise ValueError("Missing CIF cell tag(s): " + ", ".join(missing))
    cell = [float(value) for value in values if value is not None]
    formula = cif_scalar(text, ["_chemical_formula_sum", "_chemical_formula_moiety"])
    if formula is None:
        raise ValueError("Structure CIF has no chemical formula. Supply --composition C:...,H:...")
    composition = parse_formula(formula)
    z_text = cif_scalar(text, ["_cell_formula_units_Z"])
    nunit = int(float(z_text)) if z_text is not None else 1
    return title, cell, composition, nunit


def logical_shelx_lines(path: Path) -> list[str]:
    physical = path.read_text(encoding="latin-1", errors="replace").splitlines()
    logical: list[str] = []
    pending = ""
    for raw in physical:
        line = raw.strip()
        if not line:
            continue
        pending = f"{pending} {line}".strip()
        if pending.endswith("="):
            pending = pending[:-1].rstrip()
            continue
        logical.append(pending)
        pending = ""
    if pending:
        logical.append(pending)
    return logical


def parse_shelx(path: Path) -> tuple[str, list[float], OrderedDict[str, int], int]:
    lines = logical_shelx_lines(path)
    title = path.stem
    cell: list[float] | None = None
    elements: list[str] = []
    unit_counts: list[int] | None = None
    z_value: int | None = None
    for line in lines:
        fields = line.split()
        keyword = fields[0].upper()
        if keyword == "TITL":
            title = line[4:].strip() or path.stem
        elif keyword == "CELL":
            if len(fields) < 8:
                raise ValueError("SHELX CELL requires wavelength and six cell parameters.")
            cell = [float(value) for value in fields[2:8]]
        elif keyword == "ZERR" and len(fields) >= 2:
            z_value = int(round(float(fields[1])))
        elif keyword == "SFAC":
            candidates = fields[1:]
            # Element-list form: SFAC C H N O. Numerical SFAC tables are not inferred.
            if candidates and all(re.fullmatch(r"[A-Za-z]{1,2}", x) for x in candidates):
                elements = [x.upper() for x in candidates]
        elif keyword == "UNIT":
            unit_counts = [int(round(float(value))) for value in fields[1:]]
    if cell is None:
        raise ValueError("No CELL instruction found in SHELX INS/RES file.")
    if not elements:
        raise ValueError("No element-list SFAC instruction found; use --composition.")
    if unit_counts is None or len(unit_counts) != len(elements):
        raise ValueError("UNIT counts do not match the SFAC element list; use --composition.")
    # SHELX UNIT normally describes the unit-cell contents. Convert to ASU
    # counts when Z is present. Otherwise retain counts with nunit=1.
    nunit = z_value or 1
    composition: OrderedDict[str, int] = OrderedDict()
    for element, total in zip(elements, unit_counts):
        if total % nunit != 0:
            raise ValueError(f"UNIT count {total} for {element} is not divisible by Z={nunit}. "
                "Use --composition and --nunit explicitly.")
        composition[element] = total // nunit
    return title, cell, composition, nunit


INSTAT_FAMILIES = {
    # (point group, lattice centring):
    ("1", "P"):    (1, 1, 2),
    ("-1", "P"):   (1, 1, 2),

    ("2", "P"):    (3, 3, 5),
    ("2", "C"):    (3, 4, 6),
    ("m", "P"):    (4, 3, 5),
    ("m", "C"):    (4, 4, 6),
    ("2/m", "P"):  (5, 3, 5),
    ("2/m", "C"):  (5, 4, 6),

    ("222", "P"):  (6, 7, 14),
    ("222", "C"):  (6, 8, 15),
    ("222", "F"):  (6, 9, 16),

    ("mm2", "P"):  (7, 10, 14),
    ("mm2", "C"):  (7, 11, 15),
    ("mm2", "F"):  (7, 12, 16),

    ("mmm", "P"):  (8, 7, 14),
    ("mmm", "C"):  (8, 8, 15),
    ("mmm", "F"):  (8, 9, 16),}


def extract_space_group_from_title(title: str) -> str:
    """Extract a space-group symbol from a SHELX TITL description."""
    match = re.search(r"\bin\s+([PABCIFR][A-Za-z0-9_/\-\s]*)$",
        title.strip(), re.IGNORECASE,)
    if not match:
        raise ValueError(f"Could not extract space group from TITL: {title!r}")
    return match.group(1).strip()

def infer_instat_controls(space_group_symbol: str) -> list[int]:
    """allographic space group to an INSTAT comparison family."""
    try:
        import gemmi
    except ImportError as exc:
        raise RuntimeError(
            "Gemmi is required for automatic space-group mapping. "
            "Install it with: pip install gemmi"
        ) from exc
    try:
        sg = gemmi.SpaceGroup(space_group_symbol)
    except RuntimeError as exc:
        raise ValueError(
            f"Unrecognised space-group symbol: {space_group_symbol}"
        ) from exc
    crystal_system = sg.crystal_system_str()
    point_group = sg.point_group_hm()
    centring = sg.centring_type()
    if crystal_system not in {
        "triclinic",
        "monoclinic",
        "orthorhombic",}:
        raise ValueError(f"INSTAT does not support the {crystal_system} crystal system "
            f"for space group {sg.xhm()}.")
    key = (point_group, centring)
    if key not in INSTAT_FAMILIES:
        raise ValueError("No INSTAT family mapping is available for "
            f"point group {point_group} with {centring} centring.")
    ispg, iser_acentric, iser_centric = INSTAT_FAMILIES[key]
    return [ispg,iser_acentric,iser_centric,
        9,       # compare both candidates
        40,      # Fourier coefficients
        20,      # Bessel terms
        {"P": 1, "A": 2, "B": 2, "C": 2, "I": 2, "F": 4}.get(centring, 1),]

def generate_stin(structure: Path,reflection_name: str,
    args: argparse.Namespace,) -> str:
    """Generate INSTAT STIN text from structure CIF, SHELX INS, or SHELX RES."""
    suffix = structure.suffix.lower()
    if suffix == ".cif":
        title, cell, composition, detected_nunit = parse_structure_cif(structure)
    elif suffix in {".ins", ".res"}:
        title, cell, composition, detected_nunit = parse_shelx(structure)
    else:
        raise ValueError("--structure must be a .cif, .ins, or .res file.")
    if args.title:
        title = args.title
    if args.composition:
        composition = parse_composition_override(args.composition)
    nunit = args.nunit if args.nunit is not None else detected_nunit
    if not 1 <= len(composition) <= 10:
        raise ValueError("INSTAT supports between 1 and 10 atom types.")
    if nunit < 1:
        raise ValueError("nunit must be at least 1.")
    lines = [title[:80],
        " ".join(f"{value:.8g}" for value in cell),
        str(len(composition)),]
    lines.extend(f" {element:<2} {count:6d}" for element, count in composition.items())
    space_group_symbol = extract_space_group_from_title(title)
    controls = infer_instat_controls(space_group_symbol)
    lines.extend([str(nunit),reflection_name,
    " ".join(str(value) for value in controls),])
    return "\n".join(lines) + "\n"


def replace_stin_reflection_name(stin_text: str, reflection_name: str) -> str:
    """Replace the HKL filename in an existing STIN with the supplied basename."""
    lines = stin_text.splitlines()
    if len(lines) < 5:
        raise ValueError("STIN is too short.")
    try:
        nat = int(lines[2].split()[0])
        hkl_index = 4 + nat
    except (ValueError, IndexError) as exc:
        raise ValueError("Could not locate atom count in STIN.") from exc
    if hkl_index >= len(lines):
        raise ValueError("Could not locate reflection filename in STIN.")
    lines[hkl_index] = reflection_name
    return "\n".join(lines) + "\n"

def is_cif_reflection_file(path: Path) -> bool:
    """Return True if the reflection file has the tag required by INSTAT."""
    text = path.read_text(
        encoding="latin-1",
        errors="replace",)
    return "_refln_observed_status" in text.lower()

def convert_raw_shelx_hkl(raw_hkl: Path,output_hkl: Path,) -> None:
    """
    Convert raw SHELX HKLF 4 reflections to the CIF-HKL layout
    expected by INSTAT.
    Raw input:
        h  k  l  Fo_squared  sigma_Fo_squared
    INSTAT reads:
        h  k  l  Fc_squared  Fo_squared
    Fc_squared is unavailable in a raw HKLF 4 file, but INSTAT
    reads and then ignores it. Therefore, zero is written for Fc².
    """
    reflection_count = 0
    with raw_hkl.open( "r",encoding="latin-1",
        errors="replace",) as input_handle, output_hkl.open(
        "w",encoding="ascii",) as output_handle:
        output_handle.write("data_reflections\n")
        output_handle.write("#\n")
        output_handle.write("loop_\n")
        output_handle.write("_refln_index_h\n")
        output_handle.write("_refln_index_k\n")
        output_handle.write("_refln_index_l\n")
        output_handle.write("_refln_F_squared_calc\n")
        output_handle.write("_refln_F_squared_meas\n")
        output_handle.write("_refln_F_squared_sigma\n")
        output_handle.write("_refln_observed_status\n")
        for line_number, raw_line in enumerate(
            input_handle,start=1,):
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split()
            if len(fields) < 5:
                continue
            try:
                h = int(fields[0])
                k = int(fields[1])
                l = int(fields[2])
                fo_squared = float(fields[3])
                sigma = float(fields[4])
            except ValueError as exc:
                raise ValueError(
                    f"Invalid raw HKL record at line "
                    f"{line_number}: {raw_line.rstrip()}"
                ) from exc
            # Standard SHELX end-of-reflections record.
            if h == 0 and k == 0 and l == 0:
                break
            fc_squared = 0.0
            status = "o"
            # The first five fields exactly match INSTAT's:
            # FORMAT(3I4,2F12.2)
            output_handle.write(
                f"{h:4d}"
                f"{k:4d}"
                f"{l:4d}"
                f"{fc_squared:12.2f}"
                f"{fo_squared:12.2f}"
                f"{sigma:12.2f} "
                f"{status}\n")
            reflection_count += 1
    if reflection_count == 0:
        raise ValueError(f"No valid reflections were found in {raw_hkl}")
    print(f"Converted {reflection_count} raw SHELX reflections "
        f"to {output_hkl}")

def prepare_workdir(args: argparse.Namespace,) -> tuple[Path, Path]:
    if args.workdir is not None:
        workdir = args.workdir.expanduser().resolve()
    else:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        if args.stin is not None:
            input_name = args.stin.stem
        else:
            input_name = args.structure.stem

        workdir = (SCRIPT_DIR / "results" / f"{input_name}_{timestamp}")
    workdir.mkdir(parents=True,exist_ok=False,)
    supplied_hkl = require_file(args.hkl,"Reflection file",)
    # Existing CIF-HKL can be copied directly.
    if is_cif_reflection_file(supplied_hkl):
        staged_hkl = workdir / supplied_hkl.name
        shutil.copy2(supplied_hkl,staged_hkl,)
        print("Reflection format: CIF-HKL")
    # Raw SHELX HKL must be converted before INSTAT reads it.
    else:
        staged_hkl = (workdir / "MXDATA.hkl")
        convert_raw_shelx_hkl(supplied_hkl, staged_hkl,)
        print("Reflection format: raw SHELX HKL")
    if args.stin is not None:
        stin_path = require_file(args.stin,"STIN file",)
        stin_text = stin_path.read_text(
            encoding="latin-1",
            errors="replace",)
        stin_text = replace_stin_reflection_name(
            stin_text,staged_hkl.name,)
    else:
        structure = require_file(
            args.structure,"Structure file",)
        stin_text = generate_stin(
            structure,staged_hkl.name,args,)
    generated_stin = workdir / "STIN"
    generated_stin.write_text(
        stin_text,encoding="ascii",errors="replace",)
    return workdir, staged_hkl

def compile_instat(executable: Path, compiler: str) -> None:
    source = require_file(SOURCE_FILE, "INSTAT Fortran source")
    compiler_path = shutil.which(compiler) if not Path(compiler).is_file() else compiler
    if not compiler_path:
        raise RuntimeError(f"Fortran compiler not found: {compiler}")
    command = [str(compiler_path),"-O2","-std=legacy","-ffixed-line-length-72",
        str(source),"-o",str(executable),]
    print("COMMAND")
    print(" ".join(command))
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    (executable.parent / "compile.log").write_text(
        "COMMAND\n" + " ".join(command) + "\n\nSTDOUT\n" + result.stdout
        + "\nSTDERR\n" + result.stderr,encoding="utf-8",)
    if result.stdout:
        print("\nSTDOUT\n" + result.stdout)
    if result.stderr:
        print("\nSTDERR\n" + result.stderr)
    if result.returncode != 0:
        raise RuntimeError("Fortran compilation failed. See compile.log.")


def run_instat(workdir: Path, executable: Path) -> None:
    result = subprocess.run([str(executable)],cwd=workdir,
        text=True, capture_output=True, check=False,)
    (workdir / "python_stdout.txt").write_text(
        result.stdout, encoding="utf-8",)
    (workdir / "python_stderr.txt").write_text(
        result.stderr, encoding="utf-8",)
    # Record whether the INSTAT run succeeded or failed
    status_file = workdir / "RUN_STATUS.txt"
    if result.returncode == 0:
        status_file.write_text(
            "SUCCESS\n", encoding="utf-8",)
    else:
        status_file.write_text(
            f"FAILED: INSTAT return code {result.returncode}\n",
            encoding="utf-8",)
    if result.stdout:
        print("\nINSTAT STDOUT\n" + result.stdout)
    if result.stderr:
        print("\nINSTAT STDERR\n" + result.stderr)
    if result.returncode != 0:
        raise RuntimeError(
            f"INSTAT failed with return code {result.returncode}. "
            "See python_stderr.txt and SOUT in the work directory.")

def main() -> int:
    args = parse_arguments()
    try:
        workdir, staged_hkl = prepare_workdir(args)
        executable = workdir / ("instat.exe" if sys.platform.startswith("win") else "instat")
        compile_instat(executable, args.compiler)
        run_instat(workdir, executable)
        print(f"\nINSTAT completed successfully.")
        print(f"STIN:   {workdir / 'STIN'}")
        print(f"HKL:    {staged_hkl}")
        print(f"Output: {workdir}")
        return 0
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
