#!/usr/bin/env python3
"""
SAFE - Statistical Analysis of Structure Factors and Normalised E Values

prepare_sm.py

Prepare small-molecule crystallographic data for INSTAT analysis.

    Author
    ------
    Kutumbarao Nidamarthi
    Beamline Scientist, MX Beamline
    Diamond Light Source, UK
    Email: kutumbarao.nidamarthi@diamond.ac.uk

Prepare Small molecule files without modifying run_instat.py or the original inputs.

Creates:
  <prefix>_fixed.stin       native INSTAT control file pointing to fixed CIF-HKL
  <prefix>_fixed.hkl        CIF-HKL with fixed-width numeric reflection rows
  <prefix>_shelx.ins        separate SHELX instruction file
  <prefix>_shelx.hkl        raw SHELX HKLF-4 reflection file

Example:
    safe.prepare-sm \
      --stin data/sm.ins \
      --cif-hkl data/sm.hkl \
      --wavelength 0.71073 \
      --prefix data/sm

The wavelength must match the experiment. It is not inferred or guessed.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stin", required=True, type=Path)
    parser.add_argument("--cif-hkl", required=True, type=Path)
    parser.add_argument("--wavelength", required=True, type=float)
    parser.add_argument("--prefix", required=True, type=Path)
    parser.add_argument("--space-group",default="P 1",
        help="SHELX space group; default: P 1",)
    parser.add_argument("--unit-multiplier",type=int,default=None,
        help=("Multiplier for STIN element counts when writing SHELX UNIT. "
            "Default: use the NUNIT value from STIN."),)
    return parser.parse_args()

def read_stin(path: Path):
    lines = path.read_text(encoding="latin-1", errors="replace").splitlines()
    if len(lines) < 8:
        raise ValueError("STIN file is shorter than expected")
    title = lines[0].strip()
    cell_fields = lines[1].split()
    if len(cell_fields) < 6:
        raise ValueError("Could not read six cell parameters from STIN line 2")
    cell = [float(x) for x in cell_fields[:6]]
    n_elements = int(lines[2].split()[0])
    elements = []
    counts = []
    start = 3
    for line in lines[start:start + n_elements]:
        fields = line.split()
        if len(fields) < 2:
            raise ValueError(f"Invalid composition line: {line}")
        elements.append(fields[0].title())
        counts.append(int(fields[1]))
    nunit_index = start + n_elements
    nunit = int(lines[nunit_index].split()[0])
    hkl_index = nunit_index + 1
    if hkl_index >= len(lines):
        raise ValueError("Could not locate reflection filename in STIN")
    return {
        "lines": lines,
        "title": title,
        "cell": cell,
        "elements": elements,
        "counts": counts,
        "nunit": nunit,
        "hkl_index": hkl_index,}

def parse_cif_hkl(path: Path):
    lines = path.read_text(encoding="latin-1", errors="replace").splitlines()
    tags = []
    tag_start = None
    data_start = None
    for i, raw in enumerate(lines):
        text = raw.strip()
        if text.lower() == "loop_":
            candidate = []
            j = i + 1
            while j < len(lines) and lines[j].strip().startswith("_"):
                candidate.append(lines[j].strip())
                j += 1
            lower = [tag.lower() for tag in candidate]
            if "_refln_index_h" in lower and "_refln_index_k" in lower and "_refln_index_l" in lower:
                tags = candidate
                tag_start = i
                data_start = j
                break
    if not tags or data_start is None:
        raise ValueError("Could not find the _refln reflection loop")
    lookup = {tag.lower(): i for i, tag in enumerate(tags)}
    def find_tag(candidates):
        for tag in candidates:
            if tag.lower() in lookup:
                return lookup[tag.lower()]
        return None
    h_col = lookup["_refln_index_h"]
    k_col = lookup["_refln_index_k"]
    l_col = lookup["_refln_index_l"]
    fc_col = find_tag(["_refln_f_calc", "_refln_f_squared_calc"])
    fo_col = find_tag(["_refln_f_meas", "_refln_f_squared_meas"])
    sig_col = find_tag(["_refln_f_sigma", "_refln_f_squared_sigma"])
    status_col = find_tag(["_refln_observed_status", "_refln_status"])
    if fo_col is None:
        raise ValueError("Could not find measured F-squared column")
    reflections = []
    for line_number, raw in enumerate(lines[data_start:], start=data_start + 1):
        text = raw.strip()
        if not text or text.startswith("#"):
            continue
        if text.startswith("_") or text.lower() in {"loop_", "stop_"} or text.lower().startswith("data_"):
            break
        fields = text.split()
        if len(fields) < len(tags):
            continue
        try:
            h = int(fields[h_col])
            k = int(fields[k_col])
            l = int(fields[l_col])
            fc = float(fields[fc_col]) if fc_col is not None else 0.0
            fo = float(fields[fo_col])
            sig = float(fields[sig_col]) if sig_col is not None else 0.0
            status = fields[status_col] if status_col is not None else "o"
        except ValueError as exc:
            raise ValueError(f"Invalid reflection at line {line_number}: {raw}") from exc
        reflections.append((h, k, l, fc, fo, sig, status))
    if not reflections:
        raise ValueError("No reflection rows found")
    return lines, data_start, reflections


def write_fixed_cif(path: Path, lines, data_start, reflections):
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for line in lines[:data_start]:
            handle.write(line.rstrip("\r\n") + "\n")
        for h, k, l, fc, fo, sig, status in reflections:
            # The first 36 characters exactly match FORMAT(3I4,2F12.2).
            handle.write(f"{h:4d}{k:4d}{l:4d}{fc:12.2f}{fo:12.2f}"
                f"{sig:12.2f} {status}\n")

def write_fixed_stin(path: Path, stin, fixed_hkl_name: str):
    lines = list(stin["lines"])
    lines[stin["hkl_index"]] = fixed_hkl_name
    path.write_text("\n".join(lines) + "\n", encoding="ascii")

def write_shelx_hkl(path: Path, reflections):
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for h, k, l, _fc, fo, sig, _status in reflections:
            handle.write(f"{h:4d}{k:4d}{l:4d}{fo:12.2f}{sig:12.2f}\n")
        handle.write(f"{0:4d}{0:4d}{0:4d}{0.0:12.2f}{0.0:12.2f}\n")

def write_shelx_ins(path: Path, stin, wavelength: float, space_group: str, multiplier: int):
    a, b, c, alpha, beta, gamma = stin["cell"]
    sg = re.sub(r"\s+", " ", space_group.strip()).upper()
    if sg != "P 1":
        raise ValueError("This small helper currently writes symmetry instructions only for P 1. "
            "Use --space-group 'P 1' for AP1 or edit the generated INS carefully.")
    unit_counts = [count * multiplier for count in stin["counts"]]
    with path.open("w", encoding="ascii", newline="\n") as handle:
        handle.write(f"TITL {stin['title']} - SHELX copy\n")
        handle.write(f"CELL {wavelength:.6f} {a:.5f} {b:.5f} {c:.5f} "
            f"{alpha:.5f} {beta:.5f} {gamma:.5f}\n")
        handle.write("ZERR 1 0 0 0 0 0 0\n")
        handle.write("LATT -1\n")
        handle.write("SFAC " + " ".join(stin["elements"]) + "\n")
        handle.write("UNIT " + " ".join(str(x) for x in unit_counts) + "\n")
        handle.write("TREF\n")
        handle.write("HKLF 4\n")
        handle.write("END\n")

def main():
    args = parse_args()
    stin = read_stin(args.stin)
    lines, data_start, reflections = parse_cif_hkl(args.cif_hkl)
    prefix = args.prefix.expanduser().resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    fixed_hkl = Path(str(prefix) + "_fixed.hkl")
    fixed_stin = Path(str(prefix) + "_fixed.stin")
    shelx_hkl = Path(str(prefix) + "_shelx.hkl")
    shelx_ins = Path(str(prefix) + "_shelx.ins")
    write_fixed_cif(fixed_hkl, lines, data_start, reflections)
    write_fixed_stin(fixed_stin, stin, fixed_hkl.name)
    write_shelx_hkl(shelx_hkl, reflections)
    multiplier = args.unit_multiplier if args.unit_multiplier is not None else stin["nunit"]
    write_shelx_ins(shelx_ins, stin, args.wavelength, args.space_group, multiplier)
    print(f"Reflections: {len(reflections)}")
    print(f"INSTAT STIN copy: {fixed_stin}")
    print(f"INSTAT fixed CIF-HKL: {fixed_hkl}")
    print(f"SHELX INS copy: {shelx_ins}")
    print(f"SHELX HKL copy: {shelx_hkl}")
    print(f"SHELX UNIT multiplier: {multiplier}")
    print("Original files were not modified.")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}")
        raise SystemExit(1)
