# docs/DESIGN.md — Implementation Guide

> This document guides the step-by-step implementation of calculations described in `THEORY.md`.
> It is written for Claude Code (or any developer) to follow as a linear pipeline.

---

## Project Overview

**Goal:** Empirically parameterize a three-objective optimization model of antibody maturation using a 5M paired antibody sequence database (2.5M naïve + 2.5M memory), then generate the static (Pareto) and dynamic (Hamiltonian) analyses.

**Input data:**
- ~5M paired (VH + VL) antibody sequences
- Metadata per sequence: donor ID, isotype, cell compartment (naïve/memory)
- Clonal family assignments (pre-computed for memory compartment)
- Germline V(D)J gene assignments per sequence

**Core dependency:** `iggnition` (github.com/bnemoz/iggnition) for IMGT nucleotide-level alignment and numbering of all sequences. All positional analyses use iggnition-aligned coordinates.

**Language:** Python (Polars for data, NumPy/SciPy for numerics, matplotlib/plotly for visualization)

**Output:** Figures, tables, and parameter estimates that populate the Results section of the manuscript, plus serialized intermediate data for reproducibility.

---

## Repository Structure

```
antibody-maturation-optimization/
├── THEORY.md                    # Mathematical framework
├── docs/
│   └── DESIGN.md                # This file
├── data/
│   ├── raw/                     # Raw paired sequence files (parquet/tsv)
│   └── processed/               # iggnition-aligned, annotated sequences
├── src/
│   ├── __init__.py
│   ├── align.py                 # Step 0: iggnition alignment wrapper
│   ├── mutability.py            # Step 1: V_mut calculations (V1-V4)
│   ├── phi_structure.py         # Step 2: Φ_S calculations (S1-S3)
│   ├── phi_affinity.py          # Step 3: Φ_A calculations (A1-A3)
│   ├── phi_reactivity.py        # Step 4: Φ_R calculations (R1-R2)
│   ├── lagrange.py              # Step 5: λ estimation (L1-L3)
│   ├── pareto.py                # Step 6: Pareto front mapping (P1-P4)
│   ├── dynamics.py              # Step 7: Hamiltonian / lineage dynamics (D1-D4)
│   ├── biomarkers.py            # Step 8: α-β biomarkers and synthesis
│   └── utils.py                 # Shared helpers (IMGT region maps, amino acid properties, I/O)
├── notebooks/
│   ├── 01_alignment_qc.ipynb
│   ├── 02_vmut_analysis.ipynb
│   ├── 03_phi_calibration.ipynb
│   ├── 04_pareto_fronts.ipynb
│   ├── 05_dynamics.ipynb
│   └── 06_biomarkers.ipynb
├── results/
│   ├── figures/
│   └── tables/
├── tests/
│   └── test_*.py
└── pyproject.toml
```

---

## Step 0: Data Preparation and Alignment

**File:** `src/align.py`
**Purpose:** Ensure all sequences are IMGT-numbered at nucleotide resolution using `iggnition`, then build the master analysis table.

### 0.1 Run iggnition alignment

```
Input:  Raw paired sequences (VH_nt, VL_nt, VH_aa, VL_aa per row)
Output: IMGT-numbered position columns for every nucleotide and amino acid position
```

- Run `iggnition` on all VH and VL nucleotide sequences
- Output: per-sequence, per-position nucleotide identity at each IMGT-numbered position
- Store as a Polars DataFrame with columns: `seq_id`, `donor`, `isotype`, `compartment`, `clonal_family`, `v_gene_H`, `v_gene_L`, `j_gene_H`, `j_gene_L`, followed by positional columns `nt_H_1`, `nt_H_2`, ..., `nt_L_1`, `nt_L_2`, ..., and corresponding amino acid columns `aa_H_1`, ..., `aa_L_1`, ...

### 0.2 Build germline reference table

```
Input:  IMGT germline database + V/D/J gene assignments
Output: Germline nucleotide and amino acid at each IMGT position, per gene
```

- For each assigned V gene, retrieve the germline nucleotide sequence from IMGT
- Align germline to the same IMGT numbering scheme
- Store as a lookup table: `v_gene` → position → germline_nt / germline_aa

### 0.3 Compute mutation calls

```
Input:  Aligned sequences + germline reference
Output: Per-sequence, per-position binary mutation indicator + mutation type
```

For each sequence and each IMGT position:
- `is_mutated[i]` = 1 if `aa_i(x) != aa_i(g)`, else 0
- `mut_type[i]` = "replacement" if amino acid changed, "silent" if nucleotide changed but amino acid preserved, "none" otherwise
- `region[i]` = FWR1/CDR1/FWR2/CDR2/FWR3/CDR3/FWR4 (from IMGT scheme)

Compute summary columns:
- `n_mutations_H`, `n_mutations_L`, `n_mutations_total`
- `n_mut_CDR_H`, `n_mut_FWR_H`, `n_mut_CDR_L`, `n_mut_FWR_L`
- `n_replacement_CDR`, `n_silent_CDR`, `n_replacement_FWR`, `n_silent_FWR`
- `cdrh3_length` (amino acid count in CDRH3)

### 0.4 Quality control

- Remove sequences with ambiguous nucleotides (N) in V region
- Remove sequences with stop codons
- Remove sequences where V gene assignment confidence is below threshold
- Log: total sequences retained, breakdown by compartment/isotype/donor
- Save: `data/processed/aligned_master.parquet`

**Validation:** Compare mutation count distributions (naïve should peak at 0; memory should peak at 5-15 for VH).

---

## Step 1: Mutational Potential V_mut (Calculations V1–V4)

**File:** `src/mutability.py`
**Purpose:** Compute the intrinsic mutability landscape for each germline gene.

### 1.1 S5F Mutability Model (V1)

**What:** For each germline gene, compute the per-nucleotide-position mutability score from the 5-mer context.

**How:**

1. Load the S5F targeting/substitution model (from `shazam` R package data or re-derive from the naïve sequences)
   - The S5F model provides: for each 5-mer context centered on position i, the relative mutability and the substitution profile (probability of each nucleotide change)
   - If re-deriving: use ONLY synonymous mutations in memory sequences to avoid selection bias

2. For each germline gene g:
   - Extract the nucleotide sequence
   - At each position i, extract the 5-mer context c_i = (b_{i-2}, b_{i-1}, b_i, b_{i+1}, b_{i+2})
   - Look up S5F mutability: μ_i(g) = S5F_mutability(c_i)
   - Record which hotspot type applies: WRC, RGYW, WA, TW, SYC (coldspot), or neutral

3. Store: `results/tables/mutability_profiles.parquet` with columns: `v_gene`, `imgt_position`, `nucleotide`, `fivemer_context`, `mutability_score`, `hotspot_type`, `region`

**Validation:** Plot μ_i across positions for IGHV3-23 — should show peaks in CDR1/CDR2 at known overlapping AGCT hotspots.

### 1.2 AID Hotspot Enrichment (V2)

**What:** Per-germline ratio of AID hotspot density in CDRs vs FWRs.

**How:**

```python
for each germline g:
    wrc_cdr = count positions i where hotspot_type[i] in {WRC, RGYW} and region[i] in CDR
    wrc_fwr = count positions i where hotspot_type[i] in {WRC, RGYW} and region[i] in FWR
    len_cdr = count positions in CDR
    len_fwr = count positions in FWR
    R_AID(g) = (wrc_cdr / len_cdr) / (wrc_fwr / len_fwr)
```

**Output:** `results/tables/evolvability_index.csv` — one row per germline, columns: `v_gene`, `R_AID`, `n_wrc_cdr`, `n_wrc_fwr`, `total_mutability_cdr`, `total_mutability_fwr`

**Figure:** Bar plot of R_AID ranked across all IGHV genes. Highlight known "evolvable" genes (IGHV3-23, IGHV1-2, IGHV1-69).

### 1.3 CDR3 Length Effect (V3)

**What:** Does CDRH3 length correlate with total V-gene SHM load?

**How:**
- Bin memory sequences by CDRH3 length (in 2-aa bins: 8-9, 10-11, ..., 24-25, 26+)
- For each bin, compute mean and std of `n_mutations_H` (V-region only, excluding CDR3)
- Fit linear regression: `n_mutations_H_v_region ~ cdrh3_length`
- Stratify by germline to control for confounding

**Output:** `results/figures/cdrh3_length_vs_shm.png`

### 1.4 Mutation Accumulation by Isotype (V4)

**What:** Use isotype as a proxy for maturation depth.

**How:**
- Group memory sequences by isotype: IgM, IgG (aggregate subclasses initially), IgA, IgE
- Compute mean and distribution of `n_mutations_H` per isotype
- Fit: d(x) ~ Poisson(λ · t(ι)) where t(IgM) < t(IgG) < t(IgA) ≤ t(IgE)
- Estimate relative t(ι) per isotype by MLE of the Poisson model

**Output:** `results/tables/isotype_mutation_depth.csv`, `results/figures/mutation_by_isotype.png`

---

## Step 2: Structural Penalty Φ_S (Calculations S1–S3)

**File:** `src/phi_structure.py`
**Purpose:** Compute the structural constraint profile across the antibody V region.

### 2.1 Position-Specific Selection Pressure ω_i (S1)

**What:** dN/dS ratio at each IMGT position, using S5F as the neutral expectation.

**How:**

For each IMGT amino acid position i, across all memory sequences using a given germline g:

1. Count observed replacement mutations: `dN_obs(i)` = number of sequences with a replacement mutation at position i
2. Count observed silent mutations: `dS_obs(i)` = number of sequences with a silent mutation at position i  
3. Compute expected counts under neutrality:
   - `E[dN(i)]` = sum over all possible nucleotide changes at the codon containing position i, weighted by S5F mutability, that would cause an amino acid change
   - `E[dS(i)]` = same, for synonymous changes
4. Compute: `ω_i = (dN_obs / E[dN]) / (dS_obs / E[dS])`

Handle positions with zero silent mutations using a pseudocount or Bayesian shrinkage estimator.

**Output:** `results/tables/omega_per_position.parquet` — columns: `v_gene`, `imgt_position`, `region`, `omega`, `dN_obs`, `dS_obs`, `dN_exp`, `dS_exp`

**Figure:** ω profile across IMGT positions for top 10 germlines, with CDR/FWR regions shaded.

**Interpretation:**
- ω < 1 → purifying selection (structural constraint active) → high Φ_S
- ω ≈ 1 → neutral
- ω > 1 → positive selection (affinity-driven)

### 2.2 Forbidden Mutations (S2)

**What:** Positions where mutations are easy to generate but never observed (structural filter).

**How:**

```python
for each germline g, for each position i:
    expected_freq = mutability_score[i] / sum(mutability_score)  # from S5F
    observed_freq = n_mutated_at_i / n_total_sequences  # from memory data
    phi_S_filter[i] = log(expected_freq / (observed_freq + epsilon))
```

High `phi_S_filter` = position where SHM generates mutations that selection removes → strong structural constraint.

**Output:** `results/tables/forbidden_mutations.parquet`
**Figure:** Scatter plot of expected vs observed mutation frequency per position, colored by region (CDR/FWR). Points far below the diagonal are "forbidden."

### 2.3 VH/VL Interface Co-variation (S3)

**What:** Mutual information between VH and VL mutations at the interface.

**How:**

1. Define Vernier zone positions (VH: IMGT 2, 47, 48, 67, 69, 71, 78, 93, 94; VL: analogous positions — use Chothia/Al-Lazikani definitions)
2. For all paired memory sequences, extract amino acid identities at these positions
3. For each pair (i_H, j_L) of VH and VL interface positions:
   - Compute joint probability p(a, b) from observed frequencies
   - Compute marginal probabilities p(a), p(b)
   - MI(i_H, j_L) = Σ_{a,b} p(a,b) log[p(a,b) / (p(a)·p(b))]
4. Apply correction for finite sample size (Miller-Madow or Grassberger estimator)

**Output:** `results/tables/vh_vl_mutual_information.csv` — pairwise MI matrix
**Figure:** Heatmap of MI between VH and VL interface positions.

**Interpretation:** High MI = structurally coupled positions that must co-evolve. This is Φ_S^interface.

---

## Step 3: Affinity Deficit Φ_A (Calculations A1–A3)

**File:** `src/phi_affinity.py`

### 3.1 CDR Replacement Enrichment (A1)

**What:** R/S ratio in CDRs normalized by neutral expectation.

**How:**

For each memory sequence x:
```python
R_CDR = n_replacement_CDR(x)
S_CDR = n_silent_CDR(x)
RS_neutral = compute_neutral_RS(germline(x))  # from S5F model + codon table
if S_CDR > 0:
    delta_RS(x) = (R_CDR / S_CDR) / RS_neutral
else:
    delta_RS(x) = NaN  # or use Bayesian estimate
```

Compute Φ_A(x) = -log(delta_RS(x)) for each sequence (negative because high delta_RS = low affinity deficit).

**Output:** `results/tables/affinity_proxy.parquet` — per-sequence Φ_A values
**Figure:** Distribution of Φ_A across all memory sequences; overlay by isotype.

### 3.2 Convergent/Public Clonotype Analysis (A2)

**What:** Identify CDRH3s shared across ≥3 donors. Compare to private clonotypes.

**How:**

1. Extract CDRH3 amino acid sequences for all memory sequences
2. Cluster by CDRH3 identity (or 85% identity threshold for near-public)
3. Label clonotypes as "public" (found in ≥3 donors) or "private"
4. Compare:
   - Mean mutation count (public vs private)
   - Φ_S, Φ_A, Φ_R distributions
   - Germline usage enrichment

**Output:** `results/tables/public_clonotypes.csv`, `results/figures/public_vs_private.png`

### 3.3 Isotype-Stratified Affinity Proxy (A3)

**What:** Within clonal families, identify mutations enriched in IgG vs IgM.

**How:**

1. Filter clonal families that contain both IgM and IgG members
2. For each such family:
   - Identify mutations present in ≥50% of IgG members but ≤20% of IgM members
   - These are candidate "affinity-selected" mutations
3. Aggregate: what IMGT positions are enriched for affinity-selected mutations?
4. Compare to the Φ_S forbidden mutation map — affinity mutations should be at positions with low Φ_S (structurally permissive).

**Output:** `results/tables/affinity_selected_mutations.csv`
**Figure:** Positional enrichment of affinity-selected mutations across IMGT positions.

---

## Step 4: Reactivity Risk Φ_R (Calculations R1–R2)

**File:** `src/phi_reactivity.py`

### 4.1 Naïve-to-Memory Calibration (R1)

**What:** Calibrate Φ_R weights using the naïve→memory transition as a natural experiment.

**How:**

1. For all sequences (naïve + memory), compute four CDRH3 features:
   - `H_H3`: mean Kyte-Doolittle hydrophobicity of CDRH3 residues
   - `Q_H3`: net charge at pH 7.4 (K, R = +1; D, E = -1; H = +0.5)
   - `L_H3`: CDRH3 amino acid length
   - `Y_H3`: fraction of aromatic residues (F, W, Y) in CDRH3

2. Fit logistic regression:
   ```
   P(memory | x) = sigmoid(α_R · H_H3 + β_R · Q_H3 + γ_R · L_H3 + δ_R · Y_H3 + intercept)
   ```
   - The fitted coefficients (α_R, β_R, γ_R, δ_R) represent the log-odds effect of each feature on surviving GC selection.
   - **Negative coefficients** indicate features that reduce survival probability → features under negative selection → should have positive weight in Φ_R.

3. Define Φ_R(x) using the negative of the logistic regression linear predictor:
   ```
   Φ_R(x) = -(α_R · H_H3 + β_R · Q_H3 + γ_R · L_H3 + δ_R · Y_H3)
   ```
   so that higher Φ_R = higher risk of polyreactivity = more penalized.

**IMPORTANT:** Control for germline effects (some germlines like IGHV4-34 are inherently autoreactive). Include `v_gene` as a covariate, or fit separately per germline family.

**Output:** `results/tables/phi_R_coefficients.csv`, `results/figures/naive_vs_memory_phi_R.png`
**Validation:** The memory distribution of Φ_R should be shifted toward lower values compared to naïve.

### 4.2 Per-Germline Reactivity (R2)

**What:** Baseline Φ_R by germline.

**How:**
- Compute mean Φ_R per IGHV gene, separately for naïve and memory
- Compute the shift Δ⟨Φ_R⟩ = ⟨Φ_R⟩_naïve − ⟨Φ_R⟩_memory per germline
- Large Δ indicates strong anti-polyreactivity selection acting on that germline

**Output:** `results/tables/phi_R_by_germline.csv`
**Figure:** Ranked bar plot of Δ⟨Φ_R⟩ per germline. IGHV4-34 should rank high.

---

## Step 5: Lagrange Multiplier Estimation (Calculations L1–L3)

**File:** `src/lagrange.py`
**Purpose:** Estimate the constraint weights λ_S and λ_R from the observed mutation data.

### 5.1 Per-Mutation Gradient Decomposition (L1)

**What:** For each observed mutation (germline → mature), compute its effect on all three objectives. Then fit the stationarity condition to estimate λ.

**How:**

1. For each memory sequence x with germline g:
   - For each mutated position i ∈ Δ(x):
     - Compute ΔΦ_S(i): change in structural penalty when position i mutates from g_i to x_i
       - Use the ω_i value from Step 2.1: ΔΦ_S(i) ≈ -log(ω_i) (positions with ω << 1 have large ΔΦ_S)
       - Or use the forbidden mutation score from Step 2.2
     - Compute ΔΦ_A(i): change in affinity proxy
       - If position is in CDR and mutation is replacement: ΔΦ_A < 0 (affinity improving)
       - If position is in FWR and mutation is replacement: ΔΦ_A ≈ 0 (neutral for affinity)
       - Weight by position-specific affinity enrichment from Step 3.3
     - Compute ΔΦ_R(i): change in reactivity risk
       - If position is in CDRH3: compute change in H, Q+, Y contributions
       - If position is outside CDRH3: ΔΦ_R ≈ 0

2. Fit linear regression on the population of all accepted mutations:
   ```
   ΔΦ_A(i) = -λ_S · ΔΦ_S(i) - λ_R · ΔΦ_R(i) + ε
   ```
   This is the empirical version of the KKT stationarity condition.
   - Use robust regression (Huber loss) to handle outliers.
   - λ_S and λ_R should be non-negative. Use constrained regression if needed.

**Output:** `results/tables/lambda_global.csv` — global λ_S, λ_R estimates with confidence intervals

### 5.2 Per-Germline λ (L2)

**What:** Repeat L1 stratified by IGHV gene.

**How:** Same as L1, but fit separately for each IGHV gene with ≥1000 sequences.

**Output:** `results/tables/lambda_by_germline.csv`
**Figure:** Scatter plot of λ_S vs λ_R per germline, with gene labels. Look for systematic patterns:
- Are IGHV1-69 and IGHV1-2 (bnAb genes) in a distinct cluster?
- Does IGHV4-34 have anomalously high λ_R?

### 5.3 Per-Donor λ (L3)

**What:** Repeat L1 stratified by donor.

**How:** Same as L1, but fit separately per donor (only for donors with ≥5000 memory sequences).

**Output:** `results/tables/lambda_by_donor.csv`
**Figure:** Box plots of λ_S and λ_R distributions across donors. Test for inter-individual variation.

---

## Step 6: Pareto Front Mapping (Calculations P1–P4)

**File:** `src/pareto.py`
**Purpose:** Map the empirical Pareto front in the 3D objective space.

### 6.1 Compute Objective Triplets

**What:** Assign (Φ_S, Φ_A, Φ_R) to every memory sequence.

**How:**
- Φ_S(x) = composite structural penalty from Step 2 (use the per-position ω-based definition, summed across all mutated positions)
- Φ_A(x) = affinity deficit from Step 3.1
- Φ_R(x) = reactivity risk from Step 4.1

Normalize each Φ to [0, 1] range (min-max across the full memory dataset) to make them comparable.

Store: add columns `phi_S`, `phi_A`, `phi_R` to the master DataFrame.

### 6.2 Non-Dominated Sorting (P1)

**What:** Identify the Pareto front from 2.5M memory sequences.

**How:**

Use fast non-dominated sorting (Deb et al., NSGA-II):
1. Rank all sequences by Pareto dominance
2. Rank 1 = Pareto front (non-dominated set)
3. Rank 2 = non-dominated after removing rank 1, etc.

For 2.5M sequences in 3D, use an efficient implementation:
- `pygmo` library has optimized non-dominated sorting
- Or implement: for each pair, check dominance in O(k) where k=3

Store: add column `pareto_rank` to master DataFrame.

**Output:** `results/tables/pareto_front.parquet` (rank-1 sequences only)
**Figure:** 3D scatter plot of (Φ_S, Φ_A, Φ_R) for all memory sequences, with Pareto front highlighted in color. Use plotly for interactive 3D visualization.

### 6.3 Stratified Fronts (P2–P4)

**What:** Separate Pareto fronts by germline (P2), isotype (P3), and donor (P4).

**How:** Repeat 6.2 within each stratum. For each stratum, compute:
- **Front hypervolume:** volume dominated by the front (indicator of front "quality")
- **Front spread:** range of each objective on the front
- **Front centroid:** mean (Φ_S, Φ_A, Φ_R) of front members

Compare hypervolumes across germlines / isotypes / donors using statistical tests.

**Output:** `results/tables/pareto_by_germline.csv`, `results/tables/pareto_by_isotype.csv`, `results/tables/pareto_by_donor.csv`

**Key figures:**
- P2: Overlaid Pareto fronts for top 10 IGHV genes (2D projections: Φ_S vs Φ_A, Φ_S vs Φ_R, Φ_A vs Φ_R)
- P3: Pareto front shift from IgM → IgG → IgA
- P4: Donor-to-donor variation in front shape

---

## Step 7: Hamiltonian Dynamics (Calculations D1–D4)

**File:** `src/dynamics.py`
**Purpose:** Trace lineage trajectories through the 3D objective space and test Hamiltonian predictions.

### 7.1 Lineage Trajectory Mapping (D1)

**What:** For each clonal lineage, map the phylogenetic tree into (Φ_S, Φ_A, Φ_R) space.

**How:**

1. Select clonal families with ≥5 unique members (sufficient for trajectory analysis)
2. For each family:
   - The root node = germline sequence (Φ_S=0, Φ_A=baseline, Φ_R=germline_baseline)
   - Internal nodes = inferred ancestral sequences (from the existing lineage tree)
   - Leaf nodes = observed mature sequences
3. Compute (Φ_S, Φ_A, Φ_R) for every node
4. Plot the tree as a directed graph in 3D objective space

**Output:** `results/tables/lineage_trajectories.parquet` — columns: `clonal_family`, `node_id`, `parent_id`, `depth` (mutations from root), `phi_S`, `phi_A`, `phi_R`

**Figure:** 3D trajectory plots for 10 representative large lineages, with the global Pareto front as a reference surface.

### 7.2 Velocity and Acceleration (D2)

**What:** Compute the rate of objective change per mutation along each lineage edge.

**How:**

For each edge (parent → child) in the lineage tree:
```python
delta_d = mutation_count(child) - mutation_count(parent)  # typically 1-3
v_S = (phi_S(child) - phi_S(parent)) / delta_d
v_A = (phi_A(child) - phi_A(parent)) / delta_d
v_R = (phi_R(child) - phi_R(parent)) / delta_d
velocity = (v_S, v_A, v_R)
speed = norm(velocity)
```

For acceleration: compute velocity at consecutive edges and take the difference.

**Test:** Do lineages **decelerate** as they approach the Pareto front? Plot speed vs. distance-to-front. The Hamiltonian predicts that kinetic energy decreases as potential energy is minimized.

**Output:** `results/tables/lineage_velocities.parquet`
**Figure:** Speed vs. maturation depth (mutations from root), averaged across all lineages.

### 7.3 Action Comparison (D3)

**What:** Test whether observed lineage trajectories are near-optimal (minimum-action) paths.

**How:**

1. For each observed lineage trajectory, compute the discrete action:
   ```
   S_obs = Σ_t [½ ||v_t||²_M⁻¹ - U(q_t)]
   ```
   where M is the mutability mass matrix (from Step 1.1) and U = λ_S·Φ_S + λ_A·Φ_A + λ_R·Φ_R (using λ from Step 5).

2. Generate null trajectories:
   - For each observed lineage of length T, generate 1000 random walks of the same length in the mutational space M(g), using the S5F mutability model to determine step probabilities (but no selection)
   - Compute action for each random walk

3. Compare: compute p-value = fraction of random walks with action ≤ S_obs.

**Output:** `results/tables/action_comparison.csv` — per-lineage: S_obs, S_null_mean, S_null_std, p-value
**Figure:** Histogram of S_obs vs S_null distribution. If maturation is optimal, S_obs should be in the lower tail.

### 7.4 Effective Temperature (D4)

**What:** Estimate T_eff at different maturation stages.

**How:**

1. Bin all lineage edges by maturation depth (mutations from root): 1-3, 4-6, 7-9, 10-12, 13+
2. For each bin, compute the variance of velocity components: Var(v_S), Var(v_A), Var(v_R)
3. From the Langevin model: Var(v) ∝ k_B · T_eff
4. Plot T_eff vs. maturation depth

**Test:** Does T_eff decrease with depth? (i.e., does selection become more stringent over time?)

**Output:** `results/tables/effective_temperature.csv`
**Figure:** T_eff vs. maturation depth with error bars.

---

## Step 8: Biomarkers and Synthesis (α-β Analysis)

**File:** `src/biomarkers.py`

### 8.1 Strategy Vector Computation

**What:** Assign each memory sequence its maturation strategy (α, β).

**How:**
```python
alpha(x) = phi_S(x) / (phi_S(x) + phi_A(x) + phi_R(x))
beta(x)  = phi_A(x) / (phi_S(x) + phi_A(x) + phi_R(x))
gamma(x) = phi_R(x) / (phi_S(x) + phi_A(x) + phi_R(x))
# Note: alpha + beta + gamma = 1
```

**Output:** Add columns `alpha`, `beta`, `gamma` to master DataFrame.
**Figure:** Ternary plot (α, β, γ) of all memory sequences. Color by germline, isotype, or donor.

### 8.2 Strategy Profiles by Germline

**What:** Characterize each germline by its distribution on the strategy simplex.

**How:** For each IGHV gene, compute mean and covariance of (α, β, γ). Visualize as ellipses on the ternary plot.

**Output:** `results/tables/strategy_by_germline.csv`
**Figure:** Ternary plot with per-germline confidence ellipses. Key question: do bnAb germlines (IGHV1-2, IGHV1-69) cluster in a specific region?

### 8.3 Strategy Profiles by Isotype

**What:** Does isotype switching shift the strategy?

**How:** Compare (α, β, γ) distributions across IgM, IgG, IgA. Test with multivariate ANOVA.

### 8.4 λ-Signature as Germline Fingerprint

**What:** Create the composite germline-level summary.

**How:** For each germline g, compile:
- R_AID(g) — evolvability (from Step 1.2)
- λ_S(g) — structural price (from Step 5.2)
- λ_R(g) — reactivity price (from Step 5.2)
- Mean strategy (α, β, γ) (from Step 8.2)
- Pareto hypervolume (from Step 6.3)

**Output:** `results/tables/germline_fingerprints.csv`
**Figure:** Multi-panel figure showing germline fingerprints for all major IGHV genes.

---

## Execution Order Summary

```
Step 0: align.py          → data/processed/aligned_master.parquet
Step 1: mutability.py     → V_mut profiles (V1-V4)
Step 2: phi_structure.py  → Φ_S calibration (S1-S3)
Step 3: phi_affinity.py   → Φ_A calibration (A1-A3)
Step 4: phi_reactivity.py → Φ_R calibration (R1-R2)
Step 5: lagrange.py       → λ estimation (L1-L3)        [depends on Steps 2-4]
Step 6: pareto.py         → Pareto fronts (P1-P4)       [depends on Steps 2-4]
Step 7: dynamics.py       → Hamiltonian analysis (D1-D4) [depends on Steps 1-6]
Step 8: biomarkers.py     → Synthesis (α-β profiles)     [depends on Steps 5-6]
```

Steps 1–4 can run in parallel. Steps 5–8 depend on 1–4.

---

## Dependencies

```toml
[project]
name = "antibody-maturation-optimization"
requires-python = ">=3.10"
dependencies = [
    "polars>=0.20",
    "numpy>=1.24",
    "scipy>=1.11",
    "scikit-learn>=1.3",
    "matplotlib>=3.7",
    "plotly>=5.15",
    "seaborn>=0.12",
    "pygmo>=2.19",           # Pareto non-dominated sorting
    "biopython>=1.81",       # Sequence utilities
    "iggnition",             # IMGT alignment (github.com/bnemoz/iggnition)
    "statsmodels>=0.14",     # Robust regression for λ estimation
]
```

---

## Testing Strategy

Each step has a validation check (described above). Additionally:

- **Unit tests:** For each utility function (mutation calling, hotspot detection, MI computation)
- **Integration test:** Run full pipeline on a 10K-sequence subsample; verify all outputs exist and are non-empty
- **Sanity checks:**
  - Naïve sequences should have near-zero Φ_S (no mutations)
  - Memory sequences should have Φ_A < 0 on average (positive selection signal)
  - Φ_R(memory) < Φ_R(naïve) on average (negative selection for polyreactivity)
  - λ_S, λ_R should be non-negative
  - Pareto front should be a strict subset of all memory sequences
  - Lineage trajectories should generally move toward the Pareto front
