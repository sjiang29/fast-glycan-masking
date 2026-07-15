# Fast multi-site Man5 masking pipeline

This package keeps the **FastDesign sequon-creation step** from the original Rosetta workflow and replaces only the slow `GlycanTreeModeler` ensemble-generation step.

## Final workflow

```text
Input antigen
    ↓
Rosetta FastDesign for new masking sites that are not yet N-X-S/T
    ↓
Validate every requested native and engineered site
    ↓
Strip existing carbohydrate coordinates by default
    ↓
Place one sampled Man5 conformer at every requested site
    ↓
Check protein–glycan and glycan–glycan clashes jointly
    ↓
Write 100 complete multi-glycosylated antigen models
    ↓
Antibody docking
```

A native site such as `A:200` and a new masking site such as `A:87` can therefore appear in the **same output structure**, with both glycans sampled independently.

## Files

- `glycan_rotamer_generator.py` — creates a reusable Man5 conformer library from one Rosetta glycoprotein reference.
- `glycan_library_placer.py` — places the library at one or more already valid N-X-S/T sites.
- `fast_glycan_masking.py` — recommended wrapper: runs FastDesign when needed, validates sites, then calls the multi-site placer.
- `man5_sampling.yaml` — linkage-specific torsion sampling configuration.
- `rotamer_library_generator.py` — colleague-provided geometry helpers used by the generator.
- `requirements_fast_man5.txt` — Python dependencies.

## Installation

```bash
python -m pip install -r requirements_fast_man5.txt
```

---

## Step 1: generate the reusable Man5 library once

Use a Rosetta-generated structure that already contains a complete Man5 glycan at a known ASN, for example `A:200`:

```bash
python glycan_rotamer_generator.py \
  --pdb Glyc_Des_head_6uig_cut_ABC_A_0001_0001.pdb \
  --chain A \
  --resid 200 \
  --config man5_sampling.yaml \
  --n-models 10000 \
  --max-attempts 1000000 \
  --seed 2026 \
  --clash-cutoff 2.0 \
  --exclude-bonds 3 \
  --out-prefix Man5_library
```

Outputs:

```text
Man5_library.pdb   # visualization/debugging
Man5_library.npz   # main binary library used by the placer
```

The NPZ stores the glycan coordinates, atom/residue metadata, torsion values, sugar linkage topology, and the source ASN `CB/CG/ND2` attachment frame.

---

## Option A: place glycans when all sequons already exist

Use `glycan_library_placer.py` directly. Repeat `--site` for every glycan that should be present in each output model.

### One native site

```bash
time python glycan_library_placer.py \
  --library Man5_library.npz \
  --protein Des_head_6uig_cut_ABC_A_0001.pdb \
  --site A:200 \
  --n-models 100 \
  --clash-mode vdw \
  --vdw-scale 0.70 \
  --out-dir A200_fast_man5 \
  --prefix Glyc_Des_A200
```

### Native A200 plus new A87 simultaneously

This command assumes A87 already satisfies N-X-S/T:

```bash
time python glycan_library_placer.py \
  --library Man5_library.npz \
  --protein designed_antigen.pdb \
  --site A:200 \
  --site A:87 \
  --n-models 100 \
  --max-attempts 100000 \
  --clash-mode vdw \
  --vdw-scale 0.70 \
  --out-dir A200_A87_fast_man5 \
  --prefix Glyc_Des_A200_A87
```

For each final model, the program independently chooses one conformer for A200 and one for A87. It does **not** enumerate the full Cartesian product. It rejects a proposed combination if either glycan clashes with the protein or if the two glycans clash with each other.

---

## Option B: create a new sequon with FastDesign, then place all glycans

This is the recommended replacement for the old complete workflow.

```bash
time python fast_glycan_masking.py \
  --protein head_6uig_cut_ABC_A.pdb \
  --library Man5_library.npz \
  --new-site A:87 \
  --native-site A:200 \
  --n-models 100 \
  --fastdesign-template glycan_template_files/template_FastDesign.xml \
  --fastdesign-script glycan_template_files/fast_design.sh \
  --work-dir pos_87_fastdesign \
  --out-dir pos_87_fast_man5 \
  --prefix Glyc_Des_A87_A200
```

The wrapper does the following:

1. Checks the new site.
2. Creates a Rosetta resfile when mutations are required:
   - position `i` → ASN;
   - position `i+1` → ALA only when it is PRO;
   - position `i+2` → THR when it is not already SER or THR.
3. Updates the `<ReadResfile filename="...">` entry in the FastDesign XML.
4. Runs `fast_design.sh`.
5. Validates N-X-S/T after FastDesign.
6. Places and jointly filters Man5 at both A87 and A200.

Multiple new and native sites are supported:

```bash
--new-site A:87 \
--new-site A:145 \
--native-site A:200 \
--native-site B:200
```

## Sequon validation

Every requested position is checked as:

```text
position i     = N
position i+1   ≠ P
position i+2   = S or T
```

`glycan_library_placer.py` stops with an error if a site is invalid. Do not use `--skip-sequon-validation` for production structures.

## Existing/native glycans

By default the placer removes carbohydrate residues already present in the input and rebuilds all requested sites from the conformer library. This mirrors the old Rosetta XML behavior:

```xml
strip_existing="1"
```

This is appropriate when native A200 should vary across the same 100-model ensemble as the new masking glycan.

Use `--keep-existing-glycans` only when you intentionally want existing carbohydrate coordinates to remain fixed. Be careful not to request placement at a site that is already occupied, because that would duplicate the glycan.

## Clash filtering

Recommended initial setting:

```bash
--clash-mode vdw --vdw-scale 0.70
```

A clash is defined as:

```text
distance < vdw_scale × (radius atom 1 + radius atom 2)
```

Smaller values are more permissive. A practical calibration range is approximately `0.60–0.75`.

Fixed-distance mode is also available:

```bash
--clash-mode distance --clash-cutoff 2.0
```

The requested ASN residues are excluded from the protein clash test by default so that the inherited covalent attachment geometry is not rejected. All other protein atoms and all pairs of placed glycans are checked.

## Outputs

For 100 requested models:

```text
output_directory/
    Glyc_Des_..._0001.pdb
    Glyc_Des_..._0002.pdb
    ...
    Glyc_Des_..._0100.pdb
    Glyc_Des_..._placed.npz
    Glyc_Des_..._manifest.json
```

Each PDB contains the full antigen plus all requested Man5 trees. The manifest records the site list, selected library conformer for every site/model, clash settings, and acceptance statistics.

## Important assumptions and validation

- FastDesign remains responsible for creating a chemically and structurally reasonable ASN site.
- The fast method does not perform Rosetta minimization or nearby side-chain repacking after glycan placement.
- A clash-free conformation is not necessarily the global energy minimum.
- Before replacing production results, compare antibody docking from the fast and GlycanTreeModeler ensembles at several exposed and partially buried sites.
- The current implementation targets Rosetta-style N-linked Man5 trees represented by `LINK` records.
