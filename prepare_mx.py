#!/usr/bin/env python3
"""
SAFE - Statistical Analysis of Structure Factors and Normalised E Values

prepare_mx.py

Prepare protein, DNA and RNA crystallographic datasets for
INSTAT analysis from MTZ or structure-factor CIF files.

    Author
    ------
    Kutumbarao Nidamarthi
    Beamline Scientist, MX Beamline
    Diamond Light Source, UK
    Email: kutumbarao.nidamarthi@diamond.ac.uk

Prepare simple macromolecular inputs for run_instat.py.

Usage:
    python prepare_mx_instat.py mx_input.json

The JSON file supplies the MTZ, FASTA, number of ASU copies, and optional
extra elemental composition. This script writes:

    <output_dir>/<name>.ins
    <output_dir>/<name>.hkl
    <output_dir>/preparation_summary.json

The generated INS and raw HKL can be passed to the existing wrapper.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
import traceback
from collections import Counter
from pathlib import Path

import gemmi
import numpy as np


# Elemental composition of amino-acid residues inside a peptide chain.
# One H2O is added separately for the two chain termini.
RESIDUE_COMPOSITION = {
    "A": {"C": 3, "H": 5, "N": 1, "O": 1},
    "R": {"C": 6, "H": 12, "N": 4, "O": 1},
    "N": {"C": 4, "H": 6, "N": 2, "O": 2},
    "D": {"C": 4, "H": 5, "N": 1, "O": 3},
    "C": {"C": 3, "H": 5, "N": 1, "O": 1, "S": 1},
    "E": {"C": 5, "H": 7, "N": 1, "O": 3},
    "Q": {"C": 5, "H": 8, "N": 2, "O": 2},
    "G": {"C": 2, "H": 3, "N": 1, "O": 1},
    "H": {"C": 6, "H": 7, "N": 3, "O": 1},
    "I": {"C": 6, "H": 11, "N": 1, "O": 1},
    "L": {"C": 6, "H": 11, "N": 1, "O": 1},
    "K": {"C": 6, "H": 12, "N": 2, "O": 1},
    "M": {"C": 5, "H": 9, "N": 1, "O": 1, "S": 1},
    "F": {"C": 9, "H": 9, "N": 1, "O": 1},
    "P": {"C": 5, "H": 7, "N": 1, "O": 1},
    "S": {"C": 3, "H": 5, "N": 1, "O": 2},
    "T": {"C": 4, "H": 7, "N": 1, "O": 2},
    "W": {"C": 11, "H": 10, "N": 2, "O": 1},
    "Y": {"C": 9, "H": 9, "N": 1, "O": 2},
    "V": {"C": 5, "H": 9, "N": 1, "O": 1},
}


# Nucleotide residues within nucleic-acid chains.
# RNA uses A, C, G, U.
RNA_RESIDUE_COMPOSITION = {
    "A": {"C": 10, "H": 12, "N": 5, "O": 6, "P": 1},
    "C": {"C": 9,  "H": 12, "N": 3, "O": 7, "P": 1},
    "G": {"C": 10, "H": 12, "N": 5, "O": 7, "P": 1},
    "U": {"C": 9,  "H": 11, "N": 2, "O": 8, "P": 1},
}

# DNA uses DA, DC, DG, DT internally.
DNA_RESIDUE_COMPOSITION = {
    "A": {"C": 10, "H": 12, "N": 5, "O": 5, "P": 1},
    "C": {"C": 9,  "H": 12, "N": 3, "O": 6, "P": 1},
    "G": {"C": 10, "H": 12, "N": 5, "O": 6, "P": 1},
    "T": {"C": 10, "H": 13, "N": 2, "O": 7, "P": 1},
}


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_fasta_general(path: Path) -> list:
    """Read FASTA without assuming protein residue codes."""
    sequences = []
    current = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if current:
                sequences.append(
                    "".join(current).upper())
                current = []
        else:
            current.append(
                re.sub(r"\s+", "", line))
    if current:
        sequences.append("".join(current).upper())
    if not sequences:
        raise ValueError(
            f"No sequence found in {path}")
    return sequences

def sequence_composition(sequences: list[str], copies: int) -> Counter:
    """Calculate elemental composition for all protein chains in the ASU."""
    total = Counter()
    for sequence in sequences:
        for residue in sequence:
            total.update(RESIDUE_COMPOSITION[residue])
        # Add H2O for the N- and C-terminal groups of each chain.
        total.update({"H": 2, "O": 1})
    for element in list(total):
        total[element] *= copies
    return total

def nucleotide_composition(sequences: list[str],copies: int,
    polymer_type: str,) -> Counter:
    """Calculate RNA or DNA elemental composition."""
    polymer_type = polymer_type.lower()
    if polymer_type == "rna":
        residue_table = RNA_RESIDUE_COMPOSITION
    elif polymer_type == "dna":
        residue_table = DNA_RESIDUE_COMPOSITION
    else:
        raise ValueError("Nucleic-acid type must be 'rna' or 'dna'.")
    total = Counter()
    for sequence in sequences:
        invalid = sorted(set(sequence) - set(residue_table))
        if invalid:
            raise ValueError(f"Unsupported {polymer_type.upper()} residue code(s): "
                + ", ".join(invalid))
        for residue in sequence:
            total.update(residue_table[residue])
        # Terminal groups for one nucleic-acid chain.
        total.update({"H": 2, "O": 1})
    for element in list(total):
        total[element] *= copies
    return total


def add_extra_composition(total: Counter, extra: dict) -> None:
    """Add optional elemental counts supplied directly in the JSON file."""
    for element, count in extra.items():
        count = int(count)
        if count < 0:
            raise ValueError(f"Negative elemental count for {element}")
        total[element.upper()] += count


def choose_mtz_columns(mtz: gemmi.Mtz, config: dict):
    """Choose I/SIGI or F/SIGF columns, with optional JSON overrides."""
    labels = mtz.column_labels()
    intensity_pairs = [
        ("IMEAN", "SIGIMEAN"),
        ("I", "SIGI"),
        ("IOBS", "SIGIOBS"),]
    amplitude_pairs = [
        ("F", "SIGF"),
        ("FP", "SIGFP"),
        ("FOBS", "SIGFOBS"),]
    if config.get("value_column") and config.get("sigma_column"):
        value_label = config["value_column"]
        sigma_label = config["sigma_column"]
        data_kind = config.get("data_kind", "intensity").lower()
    else:
        value_label = sigma_label = data_kind = None
        for value, sigma in intensity_pairs:
            if value in labels and sigma in labels:
                value_label, sigma_label, data_kind = value, sigma, "intensity"
                break
        if value_label is None:
            for value, sigma in amplitude_pairs:
                if value in labels and sigma in labels:
                    value_label, sigma_label, data_kind = value, sigma, "amplitude"
                    break
    if value_label not in labels or sigma_label not in labels:
        raise ValueError("Could not select reflection columns. Available columns: "
            + ", ".join(labels)
            + ". Set value_column, sigma_column and data_kind in the JSON file.")
    return value_label, sigma_label, data_kind


def extract_reflections(mtz: gemmi.Mtz, value_label: str,
                        sigma_label: str, data_kind: str):
    hkl = np.asarray(mtz.make_miller_array(), dtype=int)
    values = np.asarray(mtz.column_with_label(value_label).array, dtype=float)
    sigmas = np.asarray(mtz.column_with_label(sigma_label).array, dtype=float)
    reflections = []
    for indices, value, sigma in zip(hkl, values, sigmas):
        if not np.isfinite(value) or not np.isfinite(sigma):
            continue
        if data_kind == "amplitude":
            # I = F^2 and sigma(I) is propagated as 2 F sigma(F).
            intensity = value * value
            sigma_i = abs(2.0 * value * sigma)
        elif data_kind == "intensity":
            intensity = value
            sigma_i = sigma
        else:
            raise ValueError("data_kind must be 'intensity' or 'amplitude'")
        reflections.append(
            (int(indices[0]), int(indices[1]), int(indices[2]),
             float(intensity), float(sigma_i)))
    if not reflections:
        raise ValueError("No usable reflections were found in the MTZ file")
    return reflections


def cif_value(block, tags: list[str]) -> str | None:
    """Return the first defined scalar CIF value from a list of tags."""
    for tag in tags:
        value = block.find_value(tag)
        if value is None:
            continue
        value = str(value).strip().strip("'\"")
        if value not in {"", ".", "?"}:
            return value
    return None


def cif_float(block, tags: list[str]) -> float:
    """Read a required floating-point scalar from an SF-CIF block."""
    value = cif_value(block, tags)
    if value is None:
        raise ValueError("Required SF-CIF item was not found: "
            + ", ".join(tags))
    # Remove an optional CIF standard uncertainty:
    # 93.154(2) becomes 93.154.
    value = re.sub(r"\([0-9]+\)$","",value,)
    return float(value)

def find_cif_loop(block, required_tag: str):
    """Return the CIF loop containing a required reflection tag."""
    column = block.find_loop(required_tag)
    if not column:
        raise ValueError(f"Required SF-CIF reflection tag was not found: "
            f"{required_tag}")
    loop = column.get_loop()
    if loop is None:
        raise ValueError(f"SF-CIF tag is not contained in a loop: "
            f"{required_tag}")
    return loop


def choose_sf_cif_columns(loop,config: dict,):
    """
    Choose amplitude or intensity columns from an SF-CIF reflection loop.
    Amplitudes are preferred unless value_column, sigma_column and
    data_kind are supplied in the JSON configuration.
    """
    tag_lookup = {str(tag).lower(): index
        for index, tag in enumerate(loop.tags)}
    if (config.get("value_column")
        and config.get("sigma_column")):
        value_tag = str(
            config["value_column"]).lower()
        sigma_tag = str(
            config["sigma_column"]).lower()
        data_kind = config.get("data_kind",
            "amplitude",).lower()
        # Permit short tags in the JSON file.
        if not value_tag.startswith("_refln."):
            value_tag = "_refln." + value_tag
        if not sigma_tag.startswith("_refln."):
            sigma_tag = "_refln." + sigma_tag
    else:
        candidates = [("_refln.f_meas_au",
                "_refln.f_meas_sigma_au",
                "amplitude",),
            ("_refln.f_meas",
                "_refln.f_meas_sigma",
                "amplitude",),
            ("_refln.intensity_meas",
                "_refln.intensity_sigma",
                "intensity",),
            ("_refln.intensity_meas_au",
                "_refln.intensity_sigma_au",
                "intensity",),]
        value_tag = None
        sigma_tag = None
        data_kind = None
        for candidate_value, candidate_sigma, candidate_kind in candidates:
            if (candidate_value in tag_lookup
                and candidate_sigma in tag_lookup):
                value_tag = candidate_value
                sigma_tag = candidate_sigma
                data_kind = candidate_kind
                break

    if (value_tag not in tag_lookup
        or sigma_tag not in tag_lookup):
        raise ValueError(
            "Could not select SF-CIF reflection columns. "
            "Available reflection tags: "
            + ", ".join(loop.tags)
            + ". Set value_column, sigma_column and data_kind "
            "in the JSON file if necessary.")
    if data_kind not in {
        "amplitude",
        "intensity",}:
        raise ValueError(
            "data_kind must be 'amplitude' or 'intensity'")
    return value_tag, sigma_tag, data_kind


def extract_sf_cif_reflections(path: Path,config: dict,):
    """
    Read reflections and crystallographic metadata from an SF-CIF file.
    Returns:
        reflections
        metadata
        value_tag
        sigma_tag
        data_kind
        available_tags
    """
    document = gemmi.cif.read_file(str(path))
    if len(document) == 0:
        raise ValueError(f"The SF-CIF file contains no data blocks: {path}")
    block = document.sole_block()
    loop = find_cif_loop(block,"_refln.index_h",)
    tag_lookup = {
        str(tag).lower(): index
        for index, tag in enumerate(loop.tags)}
    required_index_tags = [
        "_refln.index_h",
        "_refln.index_k",
        "_refln.index_l",]

    missing = [tag
        for tag in required_index_tags
        if tag not in tag_lookup]
    if missing:
        raise ValueError("Missing SF-CIF Miller-index tag(s): "
            + ", ".join(missing))
    (value_tag, sigma_tag, data_kind,) = choose_sf_cif_columns(
        loop,config,)
    width = loop.width()
    values = list(loop.values)
    reflections = []
    def parse_number(text):
        text = str(text).strip()
        if text in {
            "",
            ".",
            "?",}:
            return None
        text = re.sub(r"\([0-9]+\)$","",text,)
        return float(text)
    h_column = tag_lookup["_refln.index_h"]
    k_column = tag_lookup["_refln.index_k"]
    l_column = tag_lookup["_refln.index_l"]
    value_column = tag_lookup[value_tag]
    sigma_column = tag_lookup[sigma_tag]
    for start in range(0,len(values),width,):
        row = values[start:start + width]
        if len(row) != width:
            continue
        try:
            h = int(float(row[h_column]))
            k = int(float(row[k_column]))
            l = int(float(row[l_column]))
            value = parse_number(row[value_column])
            sigma = parse_number(row[sigma_column])
        except (TypeError, ValueError):
            continue
        if value is None or sigma is None:
            continue
        if not (np.isfinite(value)
            and np.isfinite(sigma)):
            continue
        if data_kind == "amplitude":
            intensity = value * value
            sigma_i = abs(2.0 * value * sigma)
        else:
            intensity = value
            sigma_i = sigma
        reflections.append((h,k,l,float(intensity),float(sigma_i),))
    if not reflections:
        raise ValueError("No usable reflections were found "
            "in the SF-CIF file.")
    cell = {
        "a": cif_float(
            block,
            ["_cell.length_a"],
        ),
        "b": cif_float(
            block,
            ["_cell.length_b"],
        ),
        "c": cif_float(
            block,
            ["_cell.length_c"],
        ),
        "alpha": cif_float(
            block,
            ["_cell.angle_alpha"],
        ),
        "beta": cif_float(
            block,
            ["_cell.angle_beta"],
        ),
        "gamma": cif_float(
            block,
            ["_cell.angle_gamma"],
        ),}
    space_group = cif_value(
        block,
        [
            "_symmetry.space_group_name_H-M",
            "_space_group.name_H-M_alt",
            "_space_group.name_H-M_full",
        ],)
    space_group_number = cif_value(
        block,
        [
            "_symmetry.Int_Tables_number",
            "_space_group.IT_number",
        ],)
    if space_group is None:
        raise ValueError(
            "The SF-CIF file does not contain a "
            "recognised space-group name.")
    metadata = {
        "cell": cell,
        "space_group": space_group,
        "space_group_number": (
            int(space_group_number)
            if space_group_number is not None
            else None),
        "data_block": block.name,}
    return (reflections,metadata,value_tag,
        sigma_tag,data_kind,list(loop.tags),)

def write_metadata_ins_sf_cif(path: Path,name: str,
    metadata: dict,composition: Counter,) -> None:
    """Write INS metadata using cell and symmetry read from SF-CIF."""
    cell = metadata["cell"]
    space_group = metadata["space_group"]
    # Convert the hexagonal-axis H notation used by some SF-CIF
    # files to the R notation expected by run_instat.py.
    if space_group.upper().startswith("H "):
        space_group = "R " + space_group[2:]
    elements = sorted(composition)
    counts = [composition[element]
        for element in elements]
    with path.open("w",encoding="ascii",) as handle:
        handle.write(f"TITL {name} in {space_group}\n")
        handle.write(
            "CELL 1.00000 "
            f"{cell['a']:.5f} "
            f"{cell['b']:.5f} "
            f"{cell['c']:.5f} "
            f"{cell['alpha']:.5f} "
            f"{cell['beta']:.5f} "
            f"{cell['gamma']:.5f}\n")
        handle.write("SFAC " + " ".join(elements) + "\n")
        handle.write("UNIT " + " ".join(str(value)
                for value in counts) + "\n")
        handle.write("HKLF 4\n")
        handle.write("END\n")


def write_raw_hkl(path: Path, reflections) -> None:
    """Write raw HKLF-4 style h k l I sigma(I) records."""
    with path.open("w", encoding="ascii") as handle:
        for h, k, l, intensity, sigma in reflections:
            handle.write(f"{h:4d}{k:4d}{l:4d}{intensity:12.2f}{sigma:12.2f}\n")
        handle.write(f"{0:4d}{0:4d}{0:4d}{0.0:12.2f}{0.0:12.2f}\n")


def write_metadata_ins(path: Path, name: str, mtz: gemmi.Mtz,
                       composition: Counter) -> None:
    """Write the simple INS metadata format already accepted by run_instat.py."""
    cell = mtz.cell
    space_group = mtz.spacegroup.short_name()
    elements = sorted(composition)
    counts = [composition[element] for element in elements]
    with path.open("w", encoding="ascii") as handle:
        handle.write(f"TITL {name} in {space_group}\n")
        handle.write(
            "CELL 1.00000 "
            f"{cell.a:.5f} {cell.b:.5f} {cell.c:.5f} "
            f"{cell.alpha:.5f} {cell.beta:.5f} {cell.gamma:.5f}\n")
        handle.write("SFAC " + " ".join(elements) + "\n")
        handle.write("UNIT " + " ".join(str(x) for x in counts) + "\n")
        handle.write("HKLF 4\n")
        handle.write("END\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    #group = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument("config",type=Path,help="Preparation JSON file")
    parser.add_argument("--dataset",help="Dataset name, e.g. protein_dna_example")
    args = parser.parse_args()
    config_path = args.config.expanduser().resolve()
    config = read_json(config_path)
    base = config_path.parent
    if args.dataset:
        datasets = config.get("datasets", [])
        matches = [
            item
            for item in datasets
            if item.get("name") == args.dataset]
        if not matches:
            available = [
            item.get("name", "<unnamed>")
            for item in datasets]
            raise ValueError(
                f"Dataset {args.dataset!r} not found in "
                f"{config_path}. "
                f"Available datasets: "
                f"{', '.join(available)}")
    config = matches[0]
    #print("\nSelected dataset:")
    #print(json.dumps(config, indent=2))
    #mtz_path = (base / config["mtz_file"]).resolve()
    if config.get("mtz_file"):
        reflection_format = "mtz"
        reflection_path = (base / config["mtz_file"]).resolve()
    elif config.get("sf_cif_file"):
        reflection_format = "sf_cif"
        reflection_path = (base / config["sf_cif_file"]).resolve()
    else:
        raise ValueError("The dataset must contain either "
            "'mtz_file' or 'sf_cif_file'.")
    protein_files = [(base / filename).resolve()
        for filename in config.get("protein_sequence", [])]
    dna_files = [(base / filename).resolve()
        for filename in config.get("dna_sequence", [])]
    rna_files = [(base / filename).resolve()
        for filename in config.get("rna_sequence", [])]
    copies = int(config.get("copies_per_asu", 1))
    # Use the dataset name when provided; otherwise use the MTZ filename.
    # Define paths first
    #mtz_path = (base / config["mtz_file"]).resolve()
    # Use dataset name when available, otherwise MTZ stem
    #name = config.get("name") or mtz_path.stem
    #output_directory_name = (config.get("output_dir")
    #    or f"{name}_prepared")
    #output_dir = (mtz_path.parent / output_directory_name).resolve()
    # Validate input
    #if not mtz_path.is_file():
    #    raise FileNotFoundError(mtz_path)
    # Create output directory
    # Use dataset name when available, otherwise reflection-file stem.
    name = (config.get("name") or reflection_path.stem)
    output_directory_name = (config.get("output_dir") or f"{name}_prepared")
    output_dir = (reflection_path.parent / output_directory_name).resolve()
    if not reflection_path.is_file():
        raise FileNotFoundError(reflection_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    # Load MTZ object
    #mtz = gemmi.read_mtz_file(str(mtz_path))
    #if mtz.spacegroup is None:
    #    raise ValueError(
    #    "The MTZ file has no space-group information")
    # Validate inputs
    # MTZ is loaded here. SF-CIF is read later after composition
    # generation, using the independent SF-CIF reader.
    mtz = None
    if reflection_format == "mtz":
        mtz = gemmi.read_mtz_file(str(reflection_path))
    #if mtz.spacegroup is None: raise ValueError(
    #        "The MTZ file has no space-group information")
    if reflection_format == "mtz":
        if mtz.spacegroup is None:
            raise ValueError("The MTZ file has no space-group information")
    if copies < 1:
        raise ValueError(
        "copies_per_asu must be at least 1")
    if not reflection_path.is_file():
        raise FileNotFoundError(reflection_path)
    #if not mtz_path.is_file():
    #    raise FileNotFoundError(mtz_path)
    all_sequence_files = protein_files + dna_files + rna_files
    if not all_sequence_files:
        raise ValueError(
        "No protein, DNA, or RNA sequence files were supplied.")
    for sequence_file in all_sequence_files:
        if not sequence_file.is_file():
            raise FileNotFoundError(sequence_file)
    # Create output directory
    #output_dir.mkdir(parents=True,exist_ok=True,)
    # Read MTZ
    #mtz = gemmi.read_mtz_file(str(mtz_path))
    #if mtz.spacegroup is None:
    #    raise ValueError(
    #        "The MTZ file has no space-group information")
    # Protein composition
    protein_composition = Counter()
    for fasta_path in protein_files:
        sequences = read_fasta_general(fasta_path)
        protein_composition.update(sequence_composition(sequences,copies,))
    # DNA composition
    dna_composition = Counter()
    for fasta_path in dna_files:
        sequences = read_fasta_general(fasta_path)
        dna_composition.update(nucleotide_composition(sequences,copies,"dna",))
    # RNA composition
    rna_composition = Counter()
    for fasta_path in rna_files:
        sequences = read_fasta_general(fasta_path)
        rna_composition.update(nucleotide_composition(sequences,copies,"rna",))
    # Combine all compositions
    final_composition = Counter()
    final_composition.update(protein_composition)
    final_composition.update(dna_composition)
    final_composition.update(rna_composition)
    # Add optional metals or other elemental counts
    add_extra_composition(final_composition,config.get("additional_components", {}),)
    # Select MTZ columns and extract reflections
    #value_label, sigma_label, data_kind = choose_mtz_columns(mtz,config,)
    #reflections = extract_reflections(mtz,value_label,sigma_label,data_kind,)
    #hkl_output = output_dir / f"{name}.hkl"
    #ins_output = output_dir / f"{name}.ins"
    #write_raw_hkl(hkl_output, reflections)
    #write_metadata_ins(ins_output, name, mtz, final_composition)
    hkl_output = (output_dir / f"{name}.hkl")
    ins_output = (output_dir / f"{name}.ins")
    if reflection_format == "mtz":
        value_label, sigma_label, data_kind = (choose_mtz_columns(mtz,config,))
        reflections = extract_reflections(mtz,value_label,sigma_label,data_kind,)
        write_raw_hkl(hkl_output,reflections,)
        write_metadata_ins(ins_output,name,mtz,final_composition,)
        space_group_name = (mtz.spacegroup.xhm())
        cell_values = [mtz.cell.a,mtz.cell.b,mtz.cell.c,mtz.cell.alpha,mtz.cell.beta,mtz.cell.gamma,]
        resolution_low = (mtz.resolution_low())
        resolution_high = (mtz.resolution_high())
        available_columns = (mtz.column_labels())
        reflection_metadata = {"input_mtz": str(reflection_path)}
    else:
        (reflections,sf_cif_metadata,value_label,sigma_label,
            data_kind,available_columns,) = extract_sf_cif_reflections(reflection_path,config,)
        write_raw_hkl(hkl_output,reflections,)
        write_metadata_ins_sf_cif(ins_output,name,sf_cif_metadata,final_composition,)
        space_group_name = (sf_cif_metadata["space_group"])
        cell = sf_cif_metadata["cell"]
        cell_values = [cell["a"],cell["b"],cell["c"],cell["alpha"],cell["beta"],cell["gamma"],]
        # Resolution is calculated from the cell and Miller indices.
        unit_cell = gemmi.UnitCell(cell["a"],cell["b"],cell["c"],cell["alpha"],cell["beta"],cell["gamma"],)
        resolution_values = []
        for h, k, l, _, _ in reflections:
            if h == 0 and k == 0 and l == 0:
                continue
            d_value = unit_cell.calculate_d([h, k, l])
        if (np.isfinite(d_value) and d_value > 0):
            resolution_values.append(float(d_value))
        resolution_low = (max(resolution_values) if resolution_values else None)
        resolution_high = ( min(resolution_values) if resolution_values else None)
        reflection_metadata = {"input_sf_cif": str(reflection_path),
                                "sf_cif_data_block": (sf_cif_metadata["data_block"]),
                                "space_group_number": (sf_cif_metadata[ "space_group_number"]),}
    summary = {"name": name,
                "system_type": config.get("system_type", "unspecified"),
                #"input_mtz": str(mtz_path),
                "reflection_format": reflection_format,"input_reflections": str(reflection_path),
                        **reflection_metadata,
                "protein_sequence_files": [str(path) for path in protein_files],
                "dna_sequence_files": [str(path) for path in dna_files],
                "rna_sequence_files": [str(path) for path in rna_files],
                "copies_per_asu": copies,
                "protein_sequence_records": len(protein_files),
                "dna_sequence_records": len(dna_files),
                "rna_sequence_records": len(rna_files),
                #"space_group": mtz.spacegroup.xhm(),
                #"cell": [mtz.cell.a,mtz.cell.b,mtz.cell.c,
                #            mtz.cell.alpha,mtz.cell.beta,mtz.cell.gamma,],
                #"resolution_low": mtz.resolution_low(),
                #"resolution_high": mtz.resolution_high(),
                #"available_mtz_columns": mtz.column_labels(),
                "space_group": space_group_name,
                "cell": cell_values,
                "resolution_low": resolution_low,
                "resolution_high": resolution_high,
                "available_reflection_columns": available_columns,
                "selected_value_column": value_label,
                "selected_sigma_column": sigma_label,
                "data_kind": data_kind,
                "reflections_written": len(reflections),
                "protein_composition": dict(protein_composition),
                "dna_composition": dict(dna_composition),
                "rna_composition": dict(rna_composition),
                "additional_components": config.get(
                "additional_components",
                {},),
                "final_composition": dict(final_composition),
                "generated_ins": str(ins_output),
                "generated_hkl": str(hkl_output),
                "notes": [("Protein, DNA and RNA compositions were calculated "
                            "from the supplied FASTA files."),
                            ("Sequence compositions were multiplied by "
                                "copies_per_asu."),
                            ("additional_components currently represents elemental "
                                "counts for the complete asymmetric unit."),
                            ("Solvent, glycans, modified residues and unlisted "
                                "ligands are not inferred automatically."),],}
   
    (output_dir / "preparation_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8")
    print(f"Prepared INS: {ins_output}")
    print(f"Prepared HKL: {hkl_output}")
    print(f"Reflections: {len(reflections)}")
    #print(f"Space group: {mtz.spacegroup.xhm()}")
    #print(f"MTZ columns: {value_label}, {sigma_label} ({data_kind})")
    print(f"Reflection format: "f"{reflection_format}")
    print(f"Space group: "f"{space_group_name}")
    print(f"Reflection columns: "f"{value_label}, "f"{sigma_label} "f"({data_kind})")
    print(f"Summary: {output_dir / 'preparation_summary.json'}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception:
    #import traceback
        traceback.print_exc()
    #except Exception as exc:
    #    print(f"ERROR: {exc}", file=sys.stderr)
    #    raise SystemExit(1)
