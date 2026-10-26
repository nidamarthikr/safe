#!/usr/bin/env python3

#!/usr/bin/env python3
"""
SAFE - Statistical Analysis of Structure Factors and Normalised E Values

fix_ins.py

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


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 3:
        print("Usage:")
        print("  python fix_ins.py input.ins output.ins")
        raise SystemExit(1)

    clean_ins(sys.argv[1], sys.argv[2])
