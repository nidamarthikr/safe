#!/usr/bin/env python3
"""
SAFE - Statistical Analysis of Structure Factors and Normalised E Values

    analyse_instat.py

Analyse INSTAT output files and summarise statistical indicators,
space-group preferences and distribution metrics.

    Author
    ------
    Kutumbarao Nidamarthi
    Beamline Scientist, MX Beamline
    Diamond Light Source, UK
    Email: kutumbarao.nidamarthi@diamond.ac.uk

Analyse an INSTAT result directory and create CentroMK-style outputs.

Usage:
    safe-analyse-instat results/lcys_20260907_124138

Expected files:
    INSTPG  numeric plot table from INSTAT
    SOUT    INSTAT text report
    STIN    optional, used for metadata

Default INSTPG column interpretation:
    0: |E| bin midpoint
    1: observed histogram
    2: exact acentric curve
    3: exact centric curve

If your INSTPG order differs, change it explicitly, for example:
    safe-analyse-instat RESULTS --columns 0,1,3,2
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt



INSTAT_SPACE_GROUPS = {
    1: ("P1", "acentric"),
    2: ("P-1", "centric"),
    3: ("P2", "acentric"),
    4: ("C2", "acentric"),
    5: ("P2/m", "centric"),
    6: ("C2/m", "centric"),
    7: ("P222", "acentric"),
    8: ("C222", "acentric"),
    9: ("F222", "acentric"),
    10: ("Pmm2", "acentric"),
    11: ("Cmm2", "acentric"),
    12: ("Fmm2", "acentric"),
    13: ("Fdd2", "acentric"),
    14: ("Pmmm", "centric"),
    15: ("Cmmm", "centric"),
    16: ("Fmmm", "centric"),
    17: ("Fddd", "centric"),
}


IDEAL_ACENTRIC_VALUES = {
    "mean_abs_E": 0.88622693,
    "mean_abs_E2_minus_1": 0.73575888,
    "percent_E_below_0_25": 6.05869,
    "percent_E_above_2": 1.83156,}

IDEAL_CENTRIC_VALUES = {
    "mean_abs_E": 0.79788456,
    "mean_abs_E2_minus_1": 0.96788290,
    "percent_E_below_0_25": 19.74127,
    "percent_E_above_2": 4.55003,}



def normalize_to_density(values, e):
    """Normalize histogram or curve so its numerical integral is one."""
    values = np.clip(np.asarray(values, dtype=float),0.0,None,)
    widths = bin_widths(e)
    integral = np.sum(values * widths)
    if integral <= 0:
        raise ValueError(
            "Cannot normalize a distribution with zero integral.")
    return values / integral


def parse_stin_candidates(path: Path):
    """Read ISER1 and ISER2 from the final control line of STIN."""

    if not path.is_file():
        return []

    lines = [
        line.strip()
        for line in path.read_text(
            encoding="latin-1", errors="replace"
        ).splitlines()
        if line.strip()]
    if not lines:
        return []
    controls = lines[-1].split()
    if len(controls) < 3:
        return []
    try:
        iser1 = int(controls[1])
        iser2 = int(controls[2])
    except ValueError:
        return []
    candidates = []
    for iser in (iser1, iser2):
        if iser == 0:
            continue
        name, character = INSTAT_SPACE_GROUPS.get(
            iser,
            (f"ISER {iser}", "unknown"),)
        candidates.append({
            "iser": iser,
            "space_group": name,
            "character": character,})
    return candidates

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result_dir", type=Path, help="Successful INSTAT result directory")
    parser.add_argument("--columns", default="0,1,2,3",
        help="INSTPG columns: E,histogram,exact_acentric,exact_centric; default 0,1,2,3",)
    parser.add_argument("--output-dir", type=Path, default=None,
        help="Output directory; default: RESULT_DIR/analysis",)
    parser.add_argument("--title", default=None, help="Optional plot title")
    return parser.parse_args()


def numeric_rows(path: Path) -> np.ndarray:
    rows = []
    for raw in path.read_text(encoding="latin-1", errors="replace").splitlines():
        text = raw.replace("D", "E").replace("d", "e")
        values = re.findall(r"(?<![A-Za-z])[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?", text)
        if len(values) >= 4:
            try:
                rows.append([float(v) for v in values])
            except ValueError:
                continue
    if not rows:
        raise ValueError(f"No numeric table was found in {path}")
    width = max(len(row) for row in rows)
    rows = [row for row in rows if len(row) == width]
    return np.asarray(rows, dtype=float)


def parse_instpg(path: Path, columns: str):
    data = numeric_rows(path)
    indices = [int(x.strip()) for x in columns.split(",")]
    if len(indices) != 4:
        raise ValueError("--columns must contain four indices")
    if max(indices) >= data.shape[1]:
        raise ValueError(f"INSTPG has {data.shape[1]} numeric columns, but --columns requests {max(indices)}")
    e, observed, exact_a, exact_c = (data[:, i] for i in indices)
    order = np.argsort(e)
    return e[order], observed[order], exact_a[order], exact_c[order], data.shape[1]


def parse_sout_r_factors(path: Path):
    text = path.read_text(encoding="latin-1", errors="replace",)
    blocks = re.split(r"(?=HISTOGRAM\s+VS)",text,flags=re.IGNORECASE,)
    found = []
    for block in blocks:
        if not re.match(
            r"HISTOGRAM\s+VS",
            block,flags=re.IGNORECASE,):
            continue
        heading_match = re.search(
            r"HISTOGRAM\s+VS\s+([^\r\n]+)",
            block,re.IGNORECASE,)
        r_match = re.search(
            r"R[- ]?FACTOR\s*\(.*?\)\s*"
            r"([-+]?\d*\.?\d+)",
            block,re.IGNORECASE,)
        if heading_match and r_match:
            heading = heading_match.group(1).strip(" -")
            found.append((heading, float(r_match.group(1))))
    result = {"raw": found,"centric": None,"acentric": None,
            "acentric_ideal": None,}
    for heading, value in found:
        heading_upper = heading.upper()
        if "IDEAL" in heading_upper:
            if ("NON-CENTRO" in heading_upper
                or "ACENT" in heading_upper):
                result["acentric_ideal"] = value
        elif ("NON-CENTRO" in heading_upper
            or "ACENT" in heading_upper):
            result["acentric"] = value
        elif "CENTRO" in heading_upper:
            result["centric"] = value
    return result

def bin_widths(x):
    if len(x) == 1:
        return np.ones(1)
    edges = np.empty(len(x) + 1)
    edges[1:-1] = 0.5 * (x[:-1] + x[1:])
    edges[0] = max(0.0, x[0] - 0.5 * (x[1] - x[0]))
    edges[-1] = x[-1] + 0.5 * (x[-1] - x[-2])
    return np.diff(edges)


def normalize_mass(values, widths=None):
    y = np.clip(np.asarray(values, dtype=float), 0.0, None)
    mass = y if widths is None else y * widths
    total = mass.sum()
    if total <= 0:
        raise ValueError("Distribution has no positive mass")
    return mass / total


def wilson_pdfs(e):
    e = np.asarray(e, dtype=float)
    centric = math.sqrt(2.0 / math.pi) * np.exp(-0.5 * e**2)
    acentric = 2.0 * e * np.exp(-e**2)
    return acentric, centric


def empirical_metrics(e, observed):
    p = normalize_mass(observed)
    n = float(np.sum(np.clip(observed, 0, None)))
    return {
        "reflection_count_from_histogram": int(round(n)),
        "mean_abs_E": float(np.sum(p * e)),
        "mean_abs_E2_minus_1": float(np.sum(p * np.abs(e**2 - 1.0))),
        "percent_E_below_0_25": float(100.0 * p[e < 0.25].sum()),
        "percent_E_above_2": float(100.0 * p[e > 2.0].sum()),
        "E_min": float(e.min()),
        "E_max": float(e.max()),
        "nonzero_bins": int(np.count_nonzero(observed > 0)),
        "number_of_bins": int(len(e)),
        "metrics_source": "histogram-bin approximation",}


def distances(e, observed, model_curve):
    widths = bin_widths(e)
    p = normalize_mass(observed)
    q = normalize_mass(model_curve, widths)
    cp = np.cumsum(p)
    cq = np.cumsum(q)
    ks = float(np.max(np.abs(cp - cq)))
    cvm = float(np.sum(q * (cp - cq)**2))
    l1 = float(np.sum(np.abs(p - q)))
    return {"KS_histogram": ks, "CvM_histogram": cvm, "L1_probability": l1}, cp, cq


def density_scaled_to_histogram(curve, e, observed):
    widths = bin_widths(e)
    target = np.sum(np.clip(observed, 0, None))
    integral = np.sum(np.clip(curve, 0, None) * widths)
    return curve * target / integral if integral > 0 else curve


def write_csv(path, e, observed, exact_a, exact_c, ideal_a, ideal_c,
              obs_cdf, exact_a_cdf, exact_c_cdf, ideal_a_cdf, ideal_c_cdf):
    headers = [
        "abs_E", "observed_histogram", "exact_acentric", "exact_centric",
        "ideal_acentric", "ideal_centric", "observed_cdf",
        "exact_acentric_cdf", "exact_centric_cdf", "ideal_acentric_cdf",
        "ideal_centric_cdf",]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(headers)
        for row in zip(e, observed, exact_a, exact_c, ideal_a, ideal_c,
                       obs_cdf, exact_a_cdf, exact_c_cdf, ideal_a_cdf, ideal_c_cdf):
            writer.writerow([f"{v:.10g}" for v in row])


def save_figure(fig, output_dir, stem):
    fig.savefig(output_dir / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(output_dir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def pdf_to_expected_counts(pdf_values, e, observed):
    """Convert probability-density values to expected bin counts."""

    widths = bin_widths(e)
    total_count = np.sum(np.clip(observed, 0.0, None))
    pdf_values = np.clip(np.asarray(pdf_values, dtype=float),
                    0.0,None,)
    probability = pdf_values * widths
    normalization = probability.sum()
    if normalization <= 0:
        raise ValueError("PDF normalization is zero.")
    probability /= normalization
    return total_count * probability


def main():
    args = parse_args()
    result_dir = args.result_dir.expanduser().resolve()
    instpg = result_dir / "INSTPG"
    sout = result_dir / "SOUT"
    stin = result_dir / "STIN"
    candidates = parse_stin_candidates(stin)
    if not instpg.is_file():
        raise FileNotFoundError(f"Missing INSTPG: {instpg}")
    if not sout.is_file():
        raise FileNotFoundError(f"Missing SOUT: {sout}")
    output_dir = (args.output_dir.expanduser().resolve() if args.output_dir
                  else result_dir / "analysis")
    output_dir.mkdir(parents=True, exist_ok=True)
    e, observed, exact_a, exact_c, source_columns = parse_instpg(instpg, args.columns)
    ideal_a_pdf, ideal_c_pdf = wilson_pdfs(e)
    # Convert all distributions to probability densities.
    observed_density = normalize_to_density(observed, e)
    exact_a_density = normalize_to_density(exact_a, e)
    exact_c_density = normalize_to_density(exact_c, e)
    # Analytical Wilson curves are already probability densities.
    ideal_a_density = ideal_a_pdf
    ideal_c_density = ideal_c_pdf
    metrics = empirical_metrics(e, observed)
    rdata = parse_sout_r_factors(sout)
    d_exact_a, obs_cdf, exact_a_cdf = distances(e, observed, exact_a)
    d_exact_c, _, exact_c_cdf = distances(e, observed, exact_c)
    d_ideal_a, _, ideal_a_cdf = distances(e, observed, ideal_a_pdf)
    d_ideal_c, _, ideal_c_cdf = distances(e, observed, ideal_c_pdf)
    comparison = {
    "exact_acentric": {
        **d_exact_a,
        "INSTAT_R": rdata["acentric"],},
    "exact_centric": {
        **d_exact_c,
        "INSTAT_R": rdata["centric"],},
    "ideal_acentric": {
        **d_ideal_a,},
    "ideal_centric": {
        **d_ideal_c,},}
    ra, rc = rdata["acentric"], rdata["centric"]
    if ra is not None and rc is not None:
        preferred = "acentric" if ra < rc else "centric" if rc < ra else "tie"
        delta_r = abs(ra - rc)
        normalized_preference = delta_r / (ra + rc) if (ra + rc) > 0 else None
    else:
        preferred = "unavailable from SOUT"
        delta_r = None
        normalized_preference = None
    # Add the new candidate mapping here
    candidate_space_groups = parse_stin_candidates(stin) 
    possible_candidate = None
    if len(candidate_space_groups) == 2:
        if preferred == "acentric":
            possible_candidate = candidate_space_groups[0]
        elif preferred == "centric":
            possible_candidate = candidate_space_groups[1]
    elif len(candidate_space_groups) == 1:
         possible_candidate = candidate_space_groups[0]
    assessment = {
        "candidate_space_groups_tested": candidate_space_groups,
        "possible_space_group_character": preferred,
        "possible_INSTAT_space_group_family": possible_candidate,
        "basis": "lowest INSTAT R-factor",
        "R_difference_absolute": delta_r,
        "R_preference_normalized": normalized_preference,
            "warning": ("Notes\n"
                    "-----\n"
                    "The reported INSTAT symbol is the symmetry family tested "
                    "by INSTAT, not necessarily the full original space group. "
                    "Systematic absences, screw axes, glide planes, twinning, "
                    "pseudosymmetry and data quality should also be examined."
                ),}
    report = {
        "result_directory": str(result_dir),
        "INSTPG_numeric_columns_detected": source_columns,
        "INSTPG_columns_used": args.columns,
        "centroMK_style_estimators": metrics,
        "INSTAT_R_factors_raw": rdata["raw"],
        "comparison": comparison,
        "assessment": assessment,
        "diagnostic_note": (
            "KS and Cramer-von Mises values are computed from binned histogram probabilities, "
            "not individual reflection-level |E| values."
        ),}
    (output_dir / "analysis_summary.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    with (output_dir / "comparison_table.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["model", "INSTAT_R", "KS_histogram", "CvM_histogram", "L1_probability"])
        for model, values in comparison.items():
            writer.writerow([model, values.get("INSTAT_R"), values["KS_histogram"],
                             values["CvM_histogram"], values["L1_probability"]])

    write_csv(output_dir / "plot_data.csv", e, observed, exact_a_density, exact_c_density,
              ideal_a_density, ideal_c_density, obs_cdf, exact_a_cdf, exact_c_cdf,
              ideal_a_cdf, ideal_c_cdf)
    widths = bin_widths(e)
    title = args.title or result_dir.name
    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    widths = bin_widths(e)
    ax.bar(e,observed_density,width=0.9 * widths,color="0.82",edgecolor="0.25",
            linewidth=0.6,label="Observed |E| histogram",)
    ax.plot(e,exact_c_density,color="#1f77b4",linewidth=2.0,label="Exact centric",)
    ax.plot(e,exact_a_density,color="#d62728",linewidth=2.0,label="Exact acentric",)
    ax.plot(e,ideal_c_density,color="#1f77b4",linestyle="--",linewidth=1.5,label="Ideal centric",)
    ax.plot(e,ideal_a_density,color="#d62728",linestyle="--",linewidth=1.5,label="Ideal acentric",)
    ax.set_xlabel("|E|")
    ax.set_ylabel("P(|E|)")
    ax.set_title(title)
    ax.set_xlim(left=0)
    ax.legend(frameon=False,loc="upper center",bbox_to_anchor=(0.48, 1.0),ncol=2,)     
    #ax.bar(e, observed, width=0.9 * widths, color="0.82", edgecolor="0.25",
    #       linewidth=0.6, label="Observed |E| histogram")
    #ax.plot(e, exact_c_scaled, color="#1f77b4", linewidth=2.0, label="Exact centric")
    #ax.plot(e, exact_a_scaled, color="#d62728", linewidth=2.0, label="Exact acentric")
    #ax.plot(e, ideal_c_scaled, color="#1f77b4", linestyle="--", linewidth=1.5,
    #        label="Ideal centric")
    #ax.plot(e, ideal_a_scaled, color="#d62728", linestyle="--", linewidth=1.5,
    #        label="Ideal acentric")
    if ra is not None and rc is not None:
        annotation = (f"$R_{{acentric}}$ = {ra:.3f}\n"
            f"$R_{{centric}}$ = {rc:.3f}\n"
            f"$\\Delta R$ = {delta_r:.3f}\n"
            f"Possible character: {preferred}")
        if possible_candidate is not None:
            annotation += ("\nPossible INSTAT family: "
                f"{possible_candidate['space_group']}")
    ax.text(0.98,0.97,annotation,transform=ax.transAxes,
            ha="right",va="top",fontsize=9,
            bbox={"boxstyle": "round,pad=0.4",
                "facecolor": "white",
                "edgecolor": "0.5",
                "alpha": 0.9,},)
    #ax.set_xlabel("|E|")
    #ax.set_ylabel("Scaled frequency")
    #ax.set_title(title)
    ax.legend(frameon=False,loc="upper center",bbox_to_anchor=(0.5, 1.0),ncol=2,)
    #ax.set_xlim(left=0)
    save_figure(fig, output_dir, "E_histogram_and_pdfs")

    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    ax.step(e, obs_cdf, where="mid", color="black", linewidth=2.0, label="Observed")
    ax.plot(e, exact_c_cdf, color="#1f77b4", linewidth=2.0, label="Exact centric")
    ax.plot(e, exact_a_cdf, color="#d62728", linewidth=2.0, label="Exact acentric")
    ax.plot(e, ideal_c_cdf, color="#1f77b4", linestyle="--", linewidth=1.5,
            label="Ideal centric")
    ax.plot(e, ideal_a_cdf, color="#d62728", linestyle="--", linewidth=1.5,
            label="Ideal acentric")
    ax.set_xlabel("|E|")
    ax.set_ylabel("Cumulative probability")
    ax.set_title(title + " cumulative distributions")
    ax.set_xlim(left=0)
    ax.set_ylim(0, 1.02)
    ax.legend(frameon=False)
    save_figure(fig, output_dir, "E_cumulative_distributions")  
    if preferred == "acentric":
        theory = IDEAL_ACENTRIC_VALUES
        theory_name = "ideal acentric"
    else:
        theory = IDEAL_CENTRIC_VALUES
        theory_name = "ideal centric"
    distance_from_theory = {
        key: abs(metrics[key] - theory[key])
        for key in theory}
    r_ratio = None
    if ra is not None and rc is not None and ra > 0:
        r_ratio = rc / ra
    with (output_dir / "assessment.txt").open("w",encoding="utf-8",) as handle:
        handle.write(f"Dataset: {title}\n")
        handle.write(f"Reflections analysed: "f"{metrics['reflection_count_from_histogram']}\n")
        handle.write(f"<|E|>: {metrics['mean_abs_E']:.5f}\n")
        handle.write(f"<|E^2 - 1|>: "
            f"{metrics['mean_abs_E2_minus_1']:.5f}\n")
        handle.write(f"% |E| < 0.25: "
            f"{metrics['percent_E_below_0_25']:.3f}\n")
        handle.write(f"% |E| > 2: "
            f"{metrics['percent_E_above_2']:.3f}\n")
        handle.write(f"INSTAT R, centric: {rc}\n")
        handle.write(f"INSTAT R, acentric: {ra}\n\n")
        handle.write("Candidate space groups:\n")
        if candidates:
            for candidate in candidates:
                handle.write(
                    f"  {candidate['space_group']} "
                    f"(ISER {candidate['iser']}, "
                    f"{candidate['character']})\n")
        else:
            handle.write("  unavailable\n")
        handle.write("\n")
        handle.write(f"Possible space-group character: {preferred}\n")
        if possible_candidate:
            handle.write(
                "Suggested candidate space group: "
                f"{possible_candidate['space_group']}\n")
        else:
            handle.write("Suggested candidate space group: undetermined\n")
        handle.write("Assessment basis: lowest INSTAT R-factor\n")
        handle.write(f"Absolute R difference: {delta_r:.3f}\n")
        handle.write(f"Normalized R preference: "
                f"{normalized_preference:.3f}\n\n")
        if r_ratio is not None:
            handle.write(f"R-factor ratio (centric/acentric): "
                f"{r_ratio:.2f}\n")
        handle.write(f"\nDistance from {theory_name} theoretical distribution\n")
        handle.write("-----------------------------------\n")
        handle.write(f"<|E|> Delta = "
            f"{distance_from_theory['mean_abs_E']:.5f}\n")
        handle.write(f"<|E^2-1|> Delta = "
            f"{distance_from_theory['mean_abs_E2_minus_1']:.5f}\n")
        handle.write(f"%|E|<0.25 Delta = "
            f"{distance_from_theory['percent_E_below_0_25']:.3f}\n")
        handle.write(f"%|E|>2 Delta = "
            f"{distance_from_theory['percent_E_above_2']:.3f}\n")
        #handle.write(assessment["warning"] + "\n")
        #handle.write("\nModel comparison:\n")
        #handle.write(f"Exact acentric KS: "
        #    f"{comparison['exact_acentric']['KS_histogram']:.5f}\n")
        #handle.write(f"Exact centric KS: "
        #    f"{comparison['exact_centric']['KS_histogram']:.5f}\n")
        #handle.write(f"Ideal acentric KS: "
        #    f"{comparison['ideal_acentric']['KS_histogram']:.5f}\n")
        #handle.write(f"Ideal centric KS: "
        #    f"{comparison['ideal_centric']['KS_histogram']:.5f}\n")         
        handle.write("\nModel discrimination metrics\n")
        handle.write("----------------------------\n")
        handle.write(f"Delta R : {delta_r:.3f}\n")
        if r_ratio is not None:
            handle.write(f"Rratio  : {r_ratio:.2f}\n")
        if delta_r > 0.50:
            interpretation = (f"Strong statistical preference for the "
                f"{preferred} model.")
        elif delta_r > 0.20:
            interpretation = (f"Moderate statistical preference for the "
                f"{preferred} model.")
        else:
            interpretation = (f"Weak statistical preference for the "
                f"{preferred} model.")
        handle.write("\nInterpretation:\n")
        handle.write(f"{interpretation}\n")
        handle.write("\nKS\n")
        handle.write(f"acentric : "
            f"{comparison['exact_acentric']['KS_histogram']:.3f}\n")
        handle.write(f"centric  : "
            f"{comparison['exact_centric']['KS_histogram']:.3f}\n")
        handle.write("\nCvM\n")
        handle.write(f"acentric : "
            f"{comparison['exact_acentric']['CvM_histogram']:.5f}\n")
        handle.write(f"centric  : "
            f"{comparison['exact_centric']['CvM_histogram']:.5f}\n")
        handle.write("\nL1\n")
        handle.write(f"acentric : "
            f"{comparison['exact_acentric']['L1_probability']:.3f}\n")
        handle.write(f"centric  : "
            f"{comparison['exact_centric']['L1_probability']:.3f}\n")
        handle.write(assessment["warning"] + "\n")        
        print(f"Analysis written to: {output_dir}")
        print(f"Possible space-group character: {preferred}")
        if possible_candidate is not None:
            print("Suggested candidate space-group family: "
                f"{possible_candidate['space_group']}")
        else:
           print("Suggested candidate space-group family: undetermined")
        if delta_r is not None:
           print(f"Absolute R difference: {delta_r:.3f}")
        if normalized_preference is not None:
           print(f"Normalized R preference: "
                 f"{normalized_preference:.3f}")
        if r_ratio is not None:
            print(f"R-factor ratio (centric/acentric): "
                    f"{r_ratio:.2f}")
        print(f"KS (acentric/centric): "
            f"{comparison['exact_acentric']['KS_histogram']:.3f} / "
            f"{comparison['exact_centric']['KS_histogram']:.3f}")
        print(f"CvM (acentric/centric): "
            f"{comparison['exact_acentric']['CvM_histogram']:.5f} / "
            f"{comparison['exact_centric']['CvM_histogram']:.5f}")
        print(f"L1 (acentric/centric): "
            f"{comparison['exact_acentric']['L1_probability']:.3f} / "
            f"{comparison['exact_centric']['L1_probability']:.3f}")    
        return 0

if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}")
        raise SystemExit(1)
