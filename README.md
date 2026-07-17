# Fast Glycan Masking

A fast conformer-library approach for generating ensembles of **Man5-glycosylated antigen structures** for antibody docking.

This project replaces the computationally expensive **Rosetta GlycanTreeModeler** ensemble generation step with a reusable glycan conformer library. Instead of rebuilding and optimizing glycans for every antigen, the pipeline generates a library of physically reasonable Man5 conformations once and rapidly places them onto one or more glycosylation sites.

The resulting glycosylated antigen ensembles can be used for downstream antibody docking and glycan masking studies.

---

## Motivation

The original glycan masking workflow consists of:

```
Protein
    │
    ▼
Rosetta FastDesign
    │
    ▼
SimpleGlycosylateMover
    │
    ▼
GlycanTreeModeler
    │
    ▼
100 glycosylated antigen models
    │
    ▼
Antibody docking
```

While accurate, **GlycanTreeModeler** is computationally expensive because every glycosylation site requires repeated glycan conformational sampling and optimization.

The goal of this project is **not** to reproduce the Rosetta energy minimum.

Instead, the goal is to efficiently generate **physically reasonable glycan conformations** suitable for antibody docking.

The new workflow is

```
Reference Man5 structure
          │
          ▼
Generate reusable conformer library
          │
          ▼
FastDesign (optional)
          │
          ▼
Place glycan conformers
          │
          ▼
Protein/Glycan clash filtering
          │
          ▼
100 glycosylated antigen models
          │
          ▼
Antibody docking
```

The conformer library is generated once and reused for all proteins.

---

# Repository structure

```
fast-glycan-masking/
│
├── README.md
├── LICENSE
├── requirements.txt
├── .gitignore
│
├── configs/
│   └── man5_sampling.yaml
│
├── src/
│   └── fast_glycan_masking/
│       ├── __init__.py
│       ├── rotamer_library_generator.py
│       ├── glycan_rotamer_generator.py
│       ├── glycan_library_placer.py
│       └── fast_glycan_masking.py
│
├── scripts/
│   ├── generate_man5_library.sh
│   └── place_multiple_glycans.sh
│
├── examples/
│   ├── reference_glycoprotein.pdb
│   ├── designed_antigen.pdb
│   └── README.md
│
├── docs/
│   └── pipeline_overview.md
│
└── tests/
```

---

`rotamer_library_generator.py` contains the generic geometry and torsion-rotation utilities used by `glycan_rotamer_generator.py`. Both files must remain inside `src/fast_glycan_masking/`. The import in `glycan_rotamer_generator.py` is package-relative:

```python
from .rotamer_library_generator import ...
```

Run the modules from the repository root with `python -m`; do not execute files inside `src/` directly.

---

# Installation

Clone the repository

```bash
git clone https://github.com/sjiang29/fast-glycan-masking.git
cd fast-glycan-masking
```

Install dependencies

```bash
pip install -r requirements.txt
```

Add the source directory to your Python path

```bash
export PYTHONPATH="$PWD/src:${PYTHONPATH:-}"
```

---

# Workflow

The workflow consists of three preparation stages followed by antibody docking.

## Stage 1 — Generate a reusable Man5 conformer library

This stage only needs to be performed once.

```bash
python -m fast_glycan_masking.glycan_rotamer_generator \
    --pdb examples/Glyc_Des_head_6uig.pdb \
    --chain A \
    --resid 200 \
    --config configs/man5_sampling.yaml \
    --n-conformers 10000 \
    --max-attempts 200000 \
    --seed 2026 \
    --out-prefix libraries/Man5_library
```

Output

```
libraries/Man5_library.npz
libraries/Man5_library.pdb
```

The generated library can be reused for all future glycan placement jobs.

`--n-conformers` controls the number of glycan-only conformers stored in the reusable library. This is different from `--n-models` in the placement stage, which controls the number of complete glycoprotein structures written. Because self-clashing samples are rejected, `--max-attempts` should normally be larger than `--n-conformers`.

For a quick test run:

```bash
python -m fast_glycan_masking.glycan_rotamer_generator \
    --pdb examples/Glyc_Des_head_6uig.pdb \
    --chain A \
    --resid 200 \
    --config configs/man5_sampling.yaml \
    --n-conformers 100 \
    --max-attempts 5000 \
    --out-prefix test_outputs/Man5_library_test
```

To inspect the detected glycan topology and torsions without generating conformers:

```bash
python -m fast_glycan_masking.glycan_rotamer_generator \
    --pdb examples/Glyc_Des_head_6uig.pdb \
    --chain A \
    --resid 200 \
    --config configs/man5_sampling.yaml \
    --report-only
```

---

## Stage 2 — Prepare glycosylation sites

New glycosylation sites must satisfy the canonical N-linked sequon

```
N-X-S/T
```

where

* X can be any residue except Proline.

If the target site already satisfies this requirement, this step can be skipped.

Otherwise, use the existing Rosetta FastDesign workflow to introduce the N-linked sequon before glycan placement.

---

## Stage 3 — Place glycans

### Single-site placement

```bash
python -m fast_glycan_masking.glycan_library_placer \
    --library libraries/Man5_library.npz \
    --protein examples/designed_antigen.pdb \
    --site A:87 \
    --n-models 100 \
    --out-dir outputs/A87 \
    --prefix Glyc_Des_A87
```

### Multiple-site placement

Native and engineered glycosylation sites can be sampled simultaneously.

```bash
python -m fast_glycan_masking.glycan_library_placer \
    --library libraries/Man5_library.npz \
    --protein examples/designed_antigen.pdb \
    --site A:87 \
    --site A:200 \
    --n-models 100 \
    --out-dir outputs/A87_A200 \
    --prefix Glyc_Des_A87_A200
```

For each output model, the program

1. randomly samples one Man5 conformer for each glycosylation site,
2. aligns each glycan onto its target ASN,
3. checks protein–glycan clashes,
4. checks glycan–glycan clashes,
5. writes only clash-free glycoprotein models.

---

## Stage 4 — Antibody docking

The generated glycosylated antigen ensemble can be directly used in downstream antibody docking workflows.

---

# Library format

Two files are generated.

### Man5_library.pdb

A multi-model PDB containing all generated Man5 conformers.

Useful for

- PyMOL
- ChimeraX
- visualization
- debugging

### Man5_library.npz

A compressed NumPy archive containing

- Cartesian coordinates
- atom names
- residue names
- residue indices
- glycan topology
- attachment frame
- sampled torsion angles

This binary format is significantly faster than repeatedly parsing PDB files during glycan placement.

---

# Sampling strategy

Unlike uniform torsion sampling, glycosidic torsion angles are sampled according to linkage-specific conformational preferences reported in the carbohydrate literature.

Each glycosidic linkage type (β1→4, α1→3, α1→6, etc.) can use independent φ/ψ/ω distributions defined in

```
configs/man5_sampling.yaml
```

This allows the conformer library to approximate experimentally observed carbohydrate conformations while avoiding computationally expensive Rosetta optimization.

---

# Troubleshooting

## `ModuleNotFoundError: No module named 'fast_glycan_masking'`

Run commands from the repository root and set:

```bash
export PYTHONPATH="$PWD/src:${PYTHONPATH:-}"
```

## `ModuleNotFoundError: No module named 'rotamer_library_generator'`

Confirm that the file is located at:

```text
src/fast_glycan_masking/rotamer_library_generator.py
```

and that `glycan_rotamer_generator.py` uses:

```python
from .rotamer_library_generator import ...
```

Run the generator as a package module:

```bash
python -m fast_glycan_masking.glycan_rotamer_generator --help
```

Do not run:

```bash
python src/fast_glycan_masking/glycan_rotamer_generator.py
```

## `unrecognized arguments: --n-conformers`

Replace `glycan_rotamer_generator.py` with the updated version in this repository. The current interface uses `--n-conformers`; the former `--n-models` spelling is retained only as a backward-compatible hidden alias.

---

# Current limitations

- Currently supports **N-linked Man5** glycans.
- Designed for antibody docking rather than atomistic free-energy calculations.
- No Rosetta minimization is performed after glycan placement.
- Placement assumes an ASN attachment geometry.
- Sampling distributions should use linkage-specific torsion conventions consistent with the chosen carbohydrate definition.

---

# Future work

Planned extensions include

- Man9 support
- Hybrid glycans
- Complex glycans
- Flexible post-placement refinement
- GPU-accelerated clash detection
- Automatic torsion distributions derived from experimentally observed glycans

---

# Citation

If this software contributes to your work, please cite

1. RosettaCarbohydrates
2. GlycanTreeModeler
3. Relevant carbohydrate conformational sampling literature

---

# Acknowledgements

This project was developed in the Meiler Laboratory at Vanderbilt University.

The original generic rotamer generation framework for small molecules was developed by a laboratory colleague and adapted here for branched glycan conformational sampling and glycan masking applications.
