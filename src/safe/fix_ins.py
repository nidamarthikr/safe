#!/usr/bin/env python3

"""
SAFE - Statistical Analysis of Structure Factors and Normalised E Values

fix_ins.py

safe-fix-ins input.ins output.ins

Inspect and repair supported INS/STIN formatting issues prior
to statistical analysis.

    Author
    ------
    Kutumbarao Nidamarthi
    Beamline Scientist, MX Beamline
    Diamond Light Source, UK
    Email: kutumbarao.nidamarthi@diamond.ac.uk
"""


from pathlib import Path

KEYWORDS = {
    "TITL", "CELL", "ZERR", "LATT", "SYMM",
    "SFAC", "UNIT", "FVAR", "HKLF", "END"
}

def clean_ins(infile, outfile):
    lines_out = []
    with open(infile, "r", encoding="latin-1") as f:
        for line in f:
            line = line.rstrip()
            if not line.strip():
                continue
            parts = line.split()
            key = parts[0].upper()
            if key in KEYWORDS:
                line = key + " " + " ".join(parts[1:])
            lines_out.append(line)
    if not any(x.startswith("END") for x in lines_out):
        lines_out.append("END")
    with open(outfile, "w") as f:
        for line in lines_out:
            f.write(line + "\n")


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(
        description="Clean and repair SHELX INS formatting for SAFE.")
    parser.add_argument("input_file",type=Path,
        help="Input SHELX INS file.",)
    parser.add_argument("output_file", type=Path,
        help="Output cleaned INS file.",)
    args = parser.parse_args()
    if not args.input_file.is_file():
        parser.error(f"Input file not found: {args.input_file}")
    clean_ins(args.input_file, args.output_file)
    print(f"Input:  {args.input_file.resolve()}")
    print(f"Output: {args.output_file.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


#if __name__ == "__main__":
#    import sys
#
#    if len(sys.argv) != 3:
#        print("Usage:")
#        raise SystemExit(1)
#
#    clean_ins(sys.argv[1], sys.argv[2])
