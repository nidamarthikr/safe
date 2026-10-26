SAFE
SAFE (Statistical Analysis of Structure Factors and Normalised E Values) provides command-line workflows for preparing, running and analysing INSTAT calculations for small-molecule and macromolecular crystallography.
The initial release supports:
existing INSTAT control files (`STIN`) with CIF-formatted reflection files;
SHELX `INS`/`RES` structures with reflection files;
macromolecular reflection data in MTZ or structure-factor CIF format;
protein, DNA and RNA sequence-based composition preparation;
compilation and execution of the original INSTAT Fortran program;
analysis of the generated INSTAT output.
Installation
A Python environment with `gemmi` is required. A Fortran compiler, normally `gfortran`, is required to compile `src/instat.f`.
```bash
python -m pip install -e .
```
Check that the commands are available:
```bash
safe-run-instat --help
safe-analyse-instat --help
safe-prepare-sm --help
safe-prepare-mx --help
safe-fix-ins --help
```

Command-line tools
```text
safe-prepare-sm       Prepare small-molecule input files
safe-prepare-mx       Prepare protein, DNA or RNA input files
safe-run-instat       Compile and run INSTAT
safe-analyse-instat   Analyse an INSTAT result directory
safe-fix-ins          Check or repair supported INS/control-file formatting
```
1. Small-molecule workflow
Existing INSTAT control file and CIF-HKL file
Use `--stin` when the input file is already in the native INSTAT control-file format.
```bash
safe-run-instat \
    --stin data/STIN \
    --hkl data/sample.hkl
```
A native INSTAT control file contains the title, unit-cell parameters, elemental composition, `NUNIT`, reflection filename and INSTAT control values. It is not a SHELX INS file.
Example structure:
```text
sample data
  10.0000  11.0000  12.0000  90.000  90.000  90.000
  4
 C   20
 H   24
 N    2
 O    4
   1
 sample.hkl
 1  1  2  8   40   20   1
```
Keep the reflection filename short because the original Fortran code may impose a legacy filename-length limit.
SHELX INS or RES with reflection data
Use `--structure` when the structure file is a SHELX `INS` or `RES` file containing records such as `TITL`, `CELL`, `LATT`, `SFAC` and `UNIT`.
```bash
safe-run-instat \
    --structure data/sample.ins \
    --hkl data/sample.hkl
```
The `TITL` line should contain a recognizable space group, for example:
```text
TITL sample in P -1
```
The `SFAC` and `UNIT` records must describe the composition used for normalisation. Confirm whether `UNIT` represents the asymmetric-unit or unit-cell contents before interpreting the E statistics.
Preparing or checking a small-molecule dataset
```bash
safe-prepare-sm --help
safe-fix-ins --help
```
Use these commands to prepare supported legacy reflection files, create separated INSTAT/SHELX copies where required, or check the input format before running INSTAT. Preserve original archived files and write converted files under new names.
2. Macromolecular workflow
Macromolecular datasets are configured in a JSON file. The input may be an MTZ file or a structure-factor CIF file.
```bash
safe-prepare-mx mx_input_example.json --dataset protein1
```
The command creates a prepared directory containing:
```text
sample.ins
sample.hkl
preparation_summary.json
```
Run INSTAT using the prepared files:
```bash
safe-run-instat \
    --structure path/to/sample.ins \
    --hkl path/to/sample.hkl
```
Editing the JSON configuration
A typical dataset entry contains:
```json
{
  "name": "protein1",
  "system_type": "protein",
  "mtz_file": "data/protein1/protein1.mtz",
  "copies_per_asu": 1,
  "protein_sequence": [
    "data/protein1/protein1.fasta"
  ],
  "dna_sequence": [],
  "rna_sequence": [],
  "additional_components": {}
}
```
For structure-factor CIF input, replace `mtz_file` with the configured SF-CIF field used by `prepare_mx.py`, for example:
```json
"sf_cif_file": "data/rna1/rna1-sf.cif"
```
Set the molecular type with one of:
```text
protein
dna
rna
```
Use `copies_per_asu` carefully. The value directly affects the composition used for normalisation and can alter the centric/acentric interpretation.
Additional ordered components may be specified when their identity and number in the asymmetric unit are known:
```json
"additional_components": {
  "MG": 2,
  "ZN": 1
}
```
Do not add every crystallisation-buffer component automatically. Include only components known or reasonably established to be present in the crystal model.
Selecting MTZ columns
SAFE attempts to identify conventional intensity or amplitude columns automatically. For non-standard MTZ files, specify the columns explicitly in the dataset configuration:
```json
"value_column": "F_xdsproc",
"sigma_column": "SIGF_xdsproc",
"data_kind": "amplitude"
```
Supported observation types are experimental amplitudes with uncertainties or experimental intensities with uncertainties. Derived map coefficients and model structure factors should not be used as primary observations.
Sequence files
Protein sequences should contain standard one-letter amino-acid codes. DNA and RNA sequences should contain standard nucleotide codes.
```text
DNA: A C G T
RNA: A C G U
```
Unknown or modified residue codes such as `N` or `X` require explicit treatment. Do not silently remove modified residues when they contribute atoms to the crystal composition.
3. Running INSTAT
SAFE supports two principal modes.
Mode 1: existing native STIN
```bash
safe-run-instat --stin data/STIN --hkl data/sample.hkl
```
Mode 2: generate STIN from a structure file
```bash
safe-run-instat --structure data/sample.ins --hkl data/sample.hkl
```
By default, a timestamped directory is written under `results/`. A specific output directory may be selected with:
```bash
safe-run-instat \
    --structure data/sample.ins \
    --hkl data/sample.hkl \
    --workdir results/sample_run
```
Typical result files include:
```text
STIN
SOUT
FFIL
EFIL
FSTLSQ
RUN_STATUS.txt
python_stdout.txt
python_stderr.txt
```
The original INSTAT implementation supports only the crystal systems implemented in `src/instat.f`. Unsupported systems are rejected rather than silently converted to unrelated symmetry.
4. Analysing an INSTAT run
Pass the completed result directory to the analysis command:
```bash
safe-analyse-instat results/sample_YYYYMMDD_HHMMSS
```
The analysis reports quantities including:
```text
possible centric or acentric character
candidate space-group family
absolute and normalised R preference
centric/acentric R-factor ratio
KS distance
Cramer-von Mises distance
L1 distance
```
The output is written to an `analysis/` directory inside the selected result directory.
The statistical assignment should be treated as supporting evidence. Ambiguous cases should also be assessed using conventional crystallographic checks, structure solution, refinement and chemical reasonableness.
5. Input checks and common issues
Use `--stin` only for a native INSTAT control file.
Use `--structure` for a SHELX `INS`/`RES` file or supported structure CIF.
Confirm that CIF reflection columns represent the intended `F`, `F^2` or intensity quantities.
Confirm that the composition and `NUNIT`/copy count are consistent with the asymmetric unit and reported `F(000)` where available.
Keep original reflection files unchanged and save reformatted copies separately.
Use short staged reflection filenames where required by the legacy Fortran source.
Check that the selected MTZ columns are experimental observations, not map coefficients.
Treat modified nucleotides, ligands and ordered metals explicitly in the composition model.
Current scope
The initial SAFE release focuses on the validated INSTAT preparation and analysis workflow. Lévy analysis, alternative scaling procedures, expanded symmetry support and modified-residue handling are planned as later extensions.
Citation
If SAFE contributes to published work, please cite the software using the metadata in `CITATION.cff` and cite the original statistical and crystallographic methods used in the analysis.
License
SAFE is distributed under the MIT License. Redistribution of third-party source code, including the original INSTAT Fortran source, may be subject to separate permissions.
