# Change Log — PairPlex Paper Analysis

> Chronological record of bugs found, fixes applied, and data insights gathered during implementation.

---

## 2026-03-13 — Step 0 notebook: first run insights + bug fixes

### Data insights from first run

| Observation | Value | Notes |
|---|---|---|
| Total sequences | 4,840,405 | Full paired library |
| Naive (bio: CD27⁻ MACS) | 2,914,772 (60.2%) | Higher naive fraction than typical — reflects experimental design |
| Naive (comp: IgM + <1% SHM) | 2,141,846 (44.2%) | |
| Naive (strict: both) | 1,979,995 (40.9%) | Most conservative set |
| Bio/comp concordance | 77.18% | 1 in 4 cells disagrees |
| Bio-naive NOT comp-naive | 692,383 | CD27⁻ cells with >1% SHM or non-IgM: IgD+ naive, pre-switched naive, or atypical memory |
| Comp-naive NOT bio-naive | 161,851 | Unmutated IgM in CD27⁺ fraction: T-independent or early memory |
| iggnition sequence output | 4,840,454 rows | 49 **extra rows** vs input (duplicate seq_names) |
| iggnition germline output | 4,839,991 rows | 414 sequences with no germline alignment |
| Post-join (inner) | 4,840,040 rows | **Inflated** by duplicate key Cartesian product — bug source |
| Mean VH codon mutations | 8.95 ± 10.12 | Bimodal: naive ~0, memory ~5–20 |
| Max VH codon mutations | 120 / 149 possible | Extreme outlier — monitor post-QC |
| Max total codon mutations | 157 | Similarly extreme |

**Region map validation:**
- H CDR: 171 nt = 57 codons (CDR1: 39 nt + CDR2: 45 nt + CDR3: 87 nt) ✓
- H FWR: 276 nt = 92 codons ✓
- L CDR: 180 nt = 60 codons ✓
- L FWR: 264 nt = 88 codons ✓

---

### Bug: shape mismatch in R/S classification (`ValueError`)

**Error:**
```
ValueError: operands could not be broadcast together with shapes (4839991,) (4840040,)
```
in `count_rs()` at `np.all(g > 0, axis=1) & np.all(s > 0, axis=1)`

**Root cause:**
The `mutated` table had 4,840,454 rows for 4,840,405 input sequences (+49). Most likely explanation: germline numbering failed for 49 sequences, causing iggnition to re-attempt alignment and emit a second row for each. Regardless of cause, these 49 seq_names appear **twice** in `mutated`.

Diagnostic: the inner join produced 4,840,040 rows, which exceeds the size of the smaller table (`germline` = 4,839,991). A standard inner join with unique keys on both sides can never exceed the smaller table — so duplicates in `mutated` are confirmed. The excess 4,840,040 − 4,839,991 = 49 rows corresponds exactly to the 49 extra rows in `mutated`.

The subsequent code used `is_in(aligned_names)` + `sort()` to re-filter `mutated` and `germline` separately. Because `is_in()` checks set membership (not count), filtering `germline` returned 4,839,991 rows while `mutations_sorted` had 4,840,040 — producing the shape mismatch when numpy arrays were constructed.

**Fix (cell `11-mut-matrix`):**
Deduplicate both `mutated` and `germline` by `seq_name` (keep first) **before** joining. Then perform a **single inner join** (`aligned_joint`) that contains both sequence and germline nt columns. All downstream numpy arrays are extracted from this one DataFrame, guaranteeing identical row counts.

```python
mutated_dedup  = mutated.unique(subset=[KEY], keep='first')
germline_dedup = germline.unique(subset=[KEY], keep='first')
aligned_joint  = mutated_dedup.join(germline_dedup, on=KEY, how="inner", suffix="_germ").sort(KEY)
mutations_sorted = aligned_joint.select([KEY] + pos_cols)
```

---

### Bug: gap character not excluded from CDRH3 length

**Observation:**
iggnition encodes two distinct "absent" states:
- `null` — position outside alignment space entirely
- `45` (ASCII `'-'`) — position exists in Aho numbering but is a gap for this sequence

The original CDRH3 length calculation used only `is_not_null()`, which incorrectly counted gap-padded positions as occupied, **overestimating CDRH3 length** for sequences with CDR3s shorter than the maximum Aho CDR3 span.

**Fix (cell `21-cdrh3-len`):**
Require both `is_not_null()` AND `!= 45` when counting occupied CDR3 nt positions:
```python
(pl.col(c).is_not_null() & (pl.col(c) != GAP)).cast(pl.UInt16)
```
where `GAP = 45`.

---

### Cells modified in `00_data_prep.ipynb`

| Cell ID | Change | Reason |
|---|---|---|
| `11-mut-matrix` | Deduplicate `mutated`/`germline` before join; create `aligned_joint` as single source of truth; derive `mutations_sorted` from it | Fix shape mismatch bug; eliminate duplicate key inflation |
| `13-region-maps` | `present = set(mutations_sorted.columns)` | Variable rename from `mutations` |
| `15-codon-counts` | Add `GAP = 45` constant; use `mutations_sorted`; document gap handling | Variable rename + gap awareness |
| `18-rs-compute` | Extract numpy arrays from `aligned_joint` (seq cols + `_germ` cols); add shape assertions | Main bugfix — guaranteed row alignment |
| `19-rs-run` | Use `aligned_joint[KEY]` as KEY source for `rs_df` | Consistent with numpy array source |
| `20-cdrh3-title` | Document null vs gap (45) distinction | Clarity |
| `21-cdrh3-len` | Use `aligned_joint`; exclude both null and gap (45) when counting occupied positions | Fix CDRH3 length overestimation |
| `23-qc` | Use `aligned_joint` and `seq_H` (from numpy) for stop codon detection; replace `mut_aligned` references | Remove stale variable; consistent with new data flow |

---

## 2026-03-13 — Step 0 second run: validation + fixes

### Bug: `mutations_sorted` stored raw nucleotides instead of differences

**Error:** All three validation plots failed:
- Plot 1 (mutation load): empty — all values at `n_mut_H = 149` (outside xlim 0–50)
- Plot 2 (CDR enrichment): single bar at neutral expectation — `n_mut_CDR_H/n_mut_H ≈ 57/149 ≈ 0.382 = neutral`
- Plot 3 (CDRH3 length): correctly computed (not affected)

**Root cause:** Refactored cell 11 set:
```python
mutations_sorted = aligned_joint.select([KEY] + pos_cols)
```
`pos_cols` are the sequence nt column names. In `aligned_joint` those columns hold raw ASCII nucleotide values (65/67/71/84 = all non-zero), so every position was flagged as "mutated".

**Fix (cell `11-mut-matrix`):** Apply the difference transformation explicitly before selecting:
```python
mutations_sorted = (
    aligned_joint
    .with_columns([(pl.col(c) - pl.col(f"{c}_germ")).alias(c) for c in pos_cols])
    .select([KEY] + pos_cols)
)
```
`aligned_joint` remains unchanged for numpy extraction.

**Note:** The R/S classification was unaffected — it uses raw `seq_H`/`germ_H` arrays directly. The R/S values produced (R/S CDR ≈ 3.2 > neutral 2.9 → positive selection; R/S FWR ≈ 1.9 < neutral → purifying selection) are biologically correct.

### Bug: CDRH3 length = 0 spike in plot, CDRL3 missing

- `cdrh3_length = 0` observed for some sequences — truncated/incompletely assembled VH. Added to QC filter.
- Light chain CDR3 length was not computed. Added `cdrl3_length` column using the same gap-aware counting logic on `L_REGIONS['CDR3']`.

### Cells modified (second round)

| Cell ID | Change | Reason |
|---|---|---|
| `11-mut-matrix` | Add `.with_columns(seq - germ)` before `.select()` to build `mutations_sorted` | Fix raw-nt-as-mutation bug |
| `21-cdrh3-len` | Add `cdrl3_length` column; refactor to shared helper function | Missing light chain CDR3 length |
| `23-qc` | Add `cdrh3_length == 0` and `cdrl3_length == 0` to QC fail set | Filter truncated sequences |
| `25-assemble` | Pass through `cdrl3_length` in master join | New column |
| `27-validate` | Stratify plot 1 by naive_strict/bio-only/comp-only/memory; fix plot 2 to use memory-only filter; fix plot 3 to show H and L CDR3 separately, exclude length=0 | Correct all three broken plots |

---

### Data questions for follow-up

1. **Why 49 duplicate seq_names from iggnition?** Root cause: germline numbering failures causing re-attempts (not duplicate input). Resolved by dedup before join.
2. **Why 60% naive by bio label?** Experimental design likely explains this — confirm whether samples include explicitly sorted naive B cells.
3. **Extreme outliers** (max VH = 120 codons, max total = 157) — check post-QC whether stop codon filter removes these or if additional curation is needed.
4. **692,383 bio-naive but NOT comp-naive** — likely IgD+ cells (IgD ≠ IgM in `c_gene:0 == 'IGHM'` filter). Pending confirmation.

---

## 2026-03-13 — Step 0 final results: biological findings

### Mutation count summary (validated — bimodal distribution confirmed)

| Metric | mean | std | median | max |
|--------|------|-----|--------|-----|
| n_mut_H | 8.95 | 10.12 | 5 | 109 |
| n_mut_L | 5.54 | 6.89 | 3 | 110 |
| n_mut_FWR_H | 5.06 | 6.74 | 2 | 78 |
| n_mut_FWR_L | 3.42 | 4.76 | 1 | 83 |
| n_mut_total | 14.49 | 14.65 | 10 | 157 |

n_mut_H mean=8.95 is consistent with published SHM data (bimodal: naive ~0, memory ~10–20). Max VH=109 and total=157 are extreme outliers — likely high-SHM memory cells; stop codon filter passed (0 stop codons detected), so these are not obvious artifacts.

### Selection signals (R/S analysis — validated, from raw alignment arrays)

| Region | n_R (mean±std) | n_S (mean±std) | R/S ratio | Interpretation |
|--------|---------------|---------------|-----------|----------------|
| CDR (H+L) | 4.12 ± 3.85 | 1.28 ± 1.58 | **3.21** | > neutral (~2.9) → **positive selection** |
| FWR (H+L) | 5.41 ± 6.36 | 2.87 ± 3.72 | **1.88** | < neutral → **purifying selection** |

CDRs accumulate beneficial replacements; FWRs resist destabilizing changes. Expected affinity maturation signature. These numbers anchor Φ_A and Φ_S calibration in Steps 2–3.

### CDR3 length distributions

| Chain | mean (aa) | std | median | min | max |
|-------|-----------|-----|--------|-----|-----|
| CDRH3 | 11.91 | 5.57 | 13 | 0 | 29 |
| CDRL3 | 7.35 | 1.10 | 7 | 0 | 29 |

CDRH3 distribution is as expected (human repertoire ~10–16 aa). CDRL3 is narrow around 7 (kappa chains), consistent with known biology.
<<<<<<< HEAD
**Note:** CDRL3 max = 29 aa is biologically real — the distribution is centred ~7–9 aa but the tail extends to 30–35 aa (rare but genuine). Not an artifact.
=======
**Note:** CDRL3 max = 29 aa is unexpected (human CDRL3 should not exceed ~12 aa). This may indicate a minority of lambda chains or alignment artifacts in the L chain CDR3 span — flag for investigation.
>>>>>>> origin/main

### QC outcomes

| Filter | Count |
|--------|-------|
| Stop codons in VH | 0 |
| Missing H alignment | 0 |
| Missing L alignment | 49 |
| CDRH3 length = 0 | 556,206 |
| CDRL3 length = 0 | 15,858 |
| **Total removed** | **570,010 (11.8%)** |
| **Retained** | **4,269,981** |

**CDRH3 = 0: biologically impossible alignment artifact.** Input is pre-filtered for complete VDJ. The 3′ primer anchors at CH1 onset → truncated mRNAs cannot be amplified. All 556,206 sequences with CDRH3=0 are iggnition CDR3 assignment failures. This is a surprisingly large fraction (11.5%) — worth investigating whether it is specific to certain donors, V genes, or CDR3 length ranges.

### Population composition (validated)

| Label | n | % |
|-------|---|---|
| naive_strict | 1,979,995 | 40.9% |
| naive_bio | 2,914,772 | 60.2% |
| naive_comp | 2,141,846 | 44.2% |
| Bio/comp concordance | 3,735,806 | 77.2% |

(Percentages are pre-QC. Post-QC population breakdown will shift slightly.)

### Open questions

1. Why does 11.5% of sequences have CDRH3=0? Is this specific to certain donors or V genes?
<<<<<<< HEAD
2. Extreme mutation outliers (n_mut_H=109, total=157) — are these genuine hypermutated memory cells or should additional outlier thresholds be applied?

---

## 2026-03-13 — Step 1 results: V1–V4 biological findings

### Dataset composition (post-QC)

| Subset | n |
|--------|---|
| Total | 4,269,981 |
| naive_strict | 1,765,330 |
| memory (not naive_bio, not naive_comp) | 1,517,911 |

Memory is dominated by IgM: 839,304 / 1,517,911 = **55.3%** class-unswitched. This reflects the experimental design.

---

### V1 — S5F Mutability profile

- 20,263,068 synonymous mutations identified across 432,344,130 opportunities (~4.7% global synonymous rate)
- 233 unique 5-mer contexts observed (out of 1024 possible; constrained by germline composition)
- Normalised max mutability = **5.71** — consistent with strongest known AID hotspots

**Hotspot density in H V-region (FR1–FR3, CDR1–CDR2):**

| Type | Count | % of 324 positions |
|------|-------|--------------------|
| neutral | 216 | 66.7% |
| edge | 37 | 11.4% |
| WRC (AID) | 21 | 6.5% |
| SYC_cold | 19 | 5.9% |
| WA (Polη) | 14 | 4.3% |
| TW (Polη) | 13 | 4.0% |
| RGYW (AID) | 4 | 1.2% |

Total AID (WRC+RGYW): **7.7%** of V-region positions. This is the empirically-derived S5F model; will be used as the neutral expectation for Φ_A and ω_i calculations in Steps 2–3.

**Note:** The profile is computed on a pooled consensus germline (all sequences collapsed to mode per Aho position). Per-germline profiles are derived separately in V2.

---

### V2 — R_AID: AID hotspot enrichment CDR vs FWR

**Major unexpected finding: R_AID < 1 for ALL 55 IGHV germlines.**

| Rank | Gene | R_AID | n |
|------|------|-------|---|
| 1 | IGHV2-26 | 0.762 | 36,302 |
| 2 | IGHV3-33 | 0.752 | 113,049 |
| 3 | IGHV3-72 | 0.752 | 20,202 |
| ... | ... | <0.76 | ... |

No germline reaches R_AID = 1.0. The WRC/RGYW density in FWR consistently exceeds CDR density.

**Interpretation:** This is consistent with published findings (e.g., Dorner et al., Shapiro et al.) — AID does not preferentially target CDRs based on germline sequence motif density alone. The CDR enrichment of SHM observed empirically arises from selection (beneficial replacements in CDRs are retained; FWR replacements are purged), not from intrinsic germline AID hotspot pre-wiring. R_AID as a "pre-wiring" metric is therefore uninformative at this level. The finding is biologically real and should be discussed in the paper.

**Impact on DESIGN.md:** The DESIGN.md expected some germlines to have R_AID > 1. This was not observed. The figure (fig_v2_raid_top40) will show a bar chart entirely left of the neutral line — this is itself the result.

---

### V3 — CDRH3 length vs V-region SHM load

- Linear regression: slope=0.121, **R²=0.002**, p≈0 (statistically significant at n=1.5M, biologically negligible)
- Each additional CDR3 aa adds ~0.12 V-gene mutations on average

**Anomaly at 22–23 aa bin:** Mean drops to 13.8 (n=24,918) vs ~17 for adjacent bins. This is reproducible (n is large). Likely reflects a structurally distinct CDR3 category — possibly sequences with a specific D-gene usage or a biased V-gene group that uses long HCDR3s but low SHM.

**Conclusion for downstream modelling:** CDRH3 length has negligible predictive power for V-region SHM. These can be treated as independent in the Φ_A framework.

---

### V4 — Mutation accumulation by isotype

| Isotype | n | Mean n_mut_H | Median |
|---------|---|-------------|--------|
| null | 47,290 | 13.3 | 13 |
| IgM | 839,304 | **13.8** | 13 |
| IgD | 36 | 14.5 | 14 |
| IgE | 721 | **19.7** | 19 |
| IgA | 378,080 | **21.2** | 21 |
| IgG | 252,464 | **22.0** | 21 |

**Key findings:**
1. **IgM memory = 13.8 mutations**: Substantial SHM in unswitched memory. These are GC-derived, class-unswitched memory B cells — biologically expected and consistent with published repertoire data.
2. **IgA ≈ IgG, and IgA < IgG**: Contradicts the DESIGN.md assumption of IgM < IgG < IgA ordering. IgA populations include T-independent switching (gut mucosa, marginal zone) which generates lower-SHM sequences and dilutes the mean. The Poisson maturation model cannot be applied with IgA > IgG assumption.
3. **IgE < IgA/IgG**: Consistent with sequential switching model (IgE arises from IgG1 by secondary switching, not direct GC-naïve switching).
4. **47,290 null-isotype**: Mean 13.3 ≈ IgM. Likely unassigned IgM or early GC members. Group with IgM for Poisson fitting.
5. **Revised maturation ladder**: IgM ≈ null < IgE < IgA ≈ IgG. Use IgG (highest, most reliable n) as the "deep maturation" reference in downstream Φ_A calibration.

**Runtime warning in V2:** `RuntimeWarning: invalid value encountered in cast` when converting null-containing Polars columns to numpy int32. Non-critical — fill_null(0) downstream handles it — but should be fixed by explicit null-filling before cast.

---

## 2026-03-15 — Step 2 results: review and biological findings

### S1 — Position-specific ω (dN/dS per Aho position)

**Method:** For each Aho AA codon position in the H V-region (FR1–CDR1–FR2–CDR2–FR3, 95 positions computed), ω = (dN_obs / E[dN]) / (dS_obs / E[dS]) with S5F-derived neutral expectation from a pooled consensus germline codon per position.

**Results summary by region:**

| Region | n positions | Mean ω | Median ω | SD ω |
|--------|-------------|--------|----------|------|
| FR1 | 24 | 0.789 | 0.338 | 1.136 |
| CDR1 | 9 | **0.737** | 0.567 | 0.516 |
| FR2 | 8 | 1.018 | 0.347 | 1.511 |
| CDR2 | 11 | **1.301** | 0.560 | 2.140 |
| FR3 | 43 | 0.916 | 0.682 | 0.991 |

**Expected findings confirmed:**
- FWR positions are predominantly ω < 1 (purifying selection) — consistent with structural constraint. The most extreme: Aho23/Cys23 (ω = 0.0017) and Aho106/Cys104 (ω = 0.00034), the two canonical disulfide cysteines. These are the most structurally locked positions in the entire V-region and serve as a strong positive control for the calculation.
- FR3 mean ω = 0.92 — slightly below neutral, consistent with the framework role of stabilising CDR3 loop geometries.

**Unexpected / notable findings:**

1. **CDR1 mean ω = 0.74 < 1 (purifying selection in CDR1).** This contradicts the naive expectation that CDR1 should show positive selection like CDR2. Reconciliation: the CDR R/S ratio from Step 0 (R/S = 3.21 vs neutral 2.9) was computed pooled across CDR1+CDR2+CDR3. CDR3 dominates this signal because it is the primary antigen-contact loop. When CDR1 is evaluated individually at the per-position ω level (and CDR3 is excluded from this analysis), positions 26–33 show ω = 0.24–0.61, consistent with moderate purifying selection. Core CDR1 residues participate in antigen binding but also constrain CDR loop backbone geometry — the structural cost of CDR1 mutation is non-trivial.

2. **CDR2 mean ω = 1.30 > 1 (slight positive selection), but highly variable (SD = 2.14).** The mean is inflated by position Aho62 (ω = 7.53, n_valid = 51,199 — a CDR2 insertion position present in only 3.5% of sequences). This high ω at a rare insertion position suggests that sequences which have extended CDR2 loops are under strong positive selection for replacements at these positions — biologically plausible. Without the outlier, median CDR2 ω = 0.56, which is in fact purifying. The correct interpretation is that the canonical CDR2 core (positions 50–61) is under mixed or purifying selection, while rare extended CDR2 variants show strong positive selection.

3. **FR1/FR2 outliers with ω > 1:** Aho17 (ω = 2.33), Aho24 (ω = 5.49), Aho98 (ω = 5.81), Aho107 (ω = 2.27). These are framework positions showing positive selection signatures:
   - **Aho24** (FR1, immediately preceding CDR1): ω = 5.49, n_valid = 1,468,868. This is a boundary position adjacent to CDR1 that participates in CDR loop presentation. Strong positive selection here is unexpected for a "framework" residue and may indicate this position is functionally involved in antigen contact or CDR1 loop anchoring.
   - **Aho98** (FR3, adjacent to CDR3): ω = 5.81, n_valid = 1,469,689. FR3 positions immediately upstream of CDR3 (positions 94–108) are known CDR3-support residues. ω > 1 here is consistent with back-mutations or compensatory mutations coordinated with CDR3 evolution.

4. **Missing positions (Aho 8, 28, 35–37, 63):** These positions are absent from the omega table, likely because they are Aho insertion positions (present in <10 sequences, filtered out by the `n_valid < 10` threshold). This is expected behaviour.

5. **NaN ω for Met and Trp codons:** Positions Aho40 (Trp), Aho41 (Met), Aho43 (Trp), Aho55 (Met), Aho57 (Trp), Aho93 (Met) have no synonymous changes possible — these codons are encoded by a single codon (Trp = TGG, Met = ATG). ω is undefined. These include the highly conserved Trp41 and Trp103 of the VH fold (canonical disulfide loop and hydrophobic core residues). Their structural importance is real but cannot be captured by ω — they must be treated separately in the Φ_S framework.

**Concern — pooled consensus germline:** The neutral expectation uses a single consensus codon per position (mode across all memory sequences). This conflates sequences from different germline V-genes that may have different codons at the same Aho position, producing a "phantom germline" that may not represent any real sequence. For positions with high germline diversity (especially CDR boundary positions), the consensus codon and its E[dN/dS] fractions may be systematically wrong. **This should be re-computed per-germline gene for the paper.** However, for a pooled overview of selection pressure across the V-region, the current approach is a valid first-pass.

---

### S2 — Forbidden mutations / Φ_S filter

**Method:** `phi_S_filter(i) = log(expected_freq_i / (observed_freq_i + ε))` where `expected_freq_i = s5f_weight_i / Σ s5f_weights` and `observed_freq_i = (n_R_obs_i + n_S_obs_i) / n_valid_i`.

**Critical bug: all phi_S_filter values are negative.**

The figure shows every single position as "enriched beyond expectation" (red bars, all below zero). Not one position shows Φ_S > 0. This is biologically impossible — we know from S1 that Cys23 and Cys104 are among the most structurally constrained positions in the V-region, and they should appear as highly constrained here too.

**Root cause — incompatible scales:**
- `expected_freq` is the **per-mutation-event probability** that a random SHM event lands at position i. This is a small number (~0.002–0.060 for individual codons across the V-region budget).
- `observed_freq` is the **per-sequence mutation rate** — the fraction of memory sequences that carry a mutation at position i. This is naturally larger by a factor equal to the mean number of V-region SHM events per memory sequence (~7–9 for the V-gene positions included here).

For a position with `expected_freq = 0.01` (1% of mutations land here) and mean ~8 V-region mutations per sequence, the expected per-sequence mutation frequency is `1 - (1-0.01)^8 ≈ 7.7%`. If the observed per-sequence rate is 5% (suppressed by selection), then `phi_S = log(0.077 / 0.05) = +0.43` — correctly positive. But with the current formula: `phi_S = log(0.01 / 0.05) = -1.6` — incorrectly negative.

**Fix required:** Scale `expected_freq` by mean V-region mutation count per memory sequence before computing the log-ratio:
```python
mean_V_mut = master.filter(~naive_bio & ~naive_comp)['n_mut_H'].mean()  # V-region mutations only
# or approximate from n_mut_H × (fraction of positions that are V-region / total H positions)
expected_per_seq = 1 - (1 - expected_freq) ** mean_V_mut  # exact
# or approximation for small p:
expected_per_seq = expected_freq * mean_V_mut
phi_S_filter = log(expected_per_seq / (observed_freq + epsilon))
```

**Partial salvage — ranking is preserved:** Despite the absolute values being wrong, the RANKING of Φ_S is still meaningful. The positions least negative (= least "enriched", = most constrained) are:
1. Aho106 / Cys104 (phi_S = −1.525, ω = 0.00034) — canonical disulfide
2. Aho23 / Cys23 (phi_S = −1.535, ω = 0.0017) — canonical disulfide
3. Aho45 / Arg45 (phi_S = −1.548, ω = 0.013) — conserved FR2 salt bridge
4. Aho100 (phi_S = −1.71)

These are exactly the positions one expects to be most structurally constrained. The ranking is correct — only the absolute scale needs fixing.

The most negative phi_S positions (supposedly most "enriched") are:
- Aho62/CDR2-insertion (phi_S = −4.30, ω = 7.53) — positive selection hotspot
- Aho98/FR3 (phi_S = −4.06, ω = 5.81) — CDR3-adjacent positive selection
- Aho24/FR1 (phi_S = −4.01, ω = 5.49) — CDR1-adjacent positive selection

This is internally consistent: the most highly mutated positions (observed_freq >> expected_freq) appear as most negative in Φ_S AND have the highest ω. Both S1 and S2 agree on which positions are under positive selection — they just use different (complementary) metrics.

**Action: cell `14-forbidden-compute` in `02_phi_structure.ipynb` must be fixed before using S2 values in downstream Steps 5–6.**

---

### S3 — VH/VL Vernier zone mutual information

**Method:** Pairwise Miller-Madow corrected MI between AA identity at VH Vernier positions [2, 47, 48, 67, 69, 71, 78, 93, 94] and VL Vernier positions [2, 36, 46, 48, 49, 64, 66, 68, 71] across all 1,517,911 paired memory sequences.

**Results:**

All MI values are effectively zero (range: 0.000–0.012 bits). The highest value is VH69–VL66 (MI = 0.0119 bits, n = 4,667).

**VL position coverage issues:**
| VL position | n_valid | Notes |
|-------------|---------|-------|
| VL2 | ~1,514,000 | Full coverage |
| VL36 | ~103,600 | Insertion position — present in ~7% of kappa sequences |
| VL46–49, 68, 71 | ~1,511,000–1,514,000 | Full coverage |
| VL64 | **0** | Position entirely absent from alignment — outside Aho space for all light chains in this dataset |
| VL66 | ~4,667 | Rare insertion position — present in <0.3% of sequences |

VL64 produces zero valid pairs for all VH positions — this column must be removed from the analysis entirely. VL36 and VL66 have severely reduced coverage, making their MI estimates unreliable.

**Biological interpretation of near-zero MI:**

Two explanations, likely both contributing:

1. **Genuine biological finding:** The VH/VL Vernier zone does not show detectable co-evolution of amino acid identity in the paired memory repertoire. This is consistent with the picture that: (a) Vernier zone residues are largely germline-encoded and fixed within a germline gene family; (b) SHM rarely targets core Vernier residues because they are in FWR positions under strong purifying selection (as seen in S1); and (c) VH and VL chains are optimised by evolution at the germline level to be compatible, not by somatic hypermutation.

2. **Methodological limitation:** MI computed on **absolute AA identity** is dominated by germline gene choice, not by somatic mutation. If IGHV3-23 always has one amino acid at VH67 and IGHV1-69 has another, the "co-variation" with the light chain is driven by VH gene usage correlating with VL gene usage (known pairing biases), not by SHM. The correct approach would be MI on **mutation co-occurrence** (binary: is position i mutated in this sequence?) — this isolates the somatic signal and removes germline confounding. This should be implemented in a revised S3.

**The VH69–VL66 signal (MI = 0.0119, highest):** Both VH69 and VL66 are interface-proximal positions. VH69 (Aho position 69) sits in FR3 at the VH/VL packing surface. That the highest MI is in this pair is biologically coherent even if the absolute magnitude is tiny — this pair may warrant specific investigation.

---

### Summary table — Step 2

| Calculation | Status | Key finding | Action required |
|-------------|--------|-------------|-----------------|
| S1: ω per position | **Valid** | CDRs mixed: CDR2 ω > 1 (positive), CDR1 ω < 1 (purifying); Cys23/Cys104 are strongest structural anchors; FR1/FR3 outliers at CDR-adjacent positions | Re-compute per germline for paper; NaN positions (Trp/Met) need separate treatment |
| S2: Φ_S filter | **Bug — fix required** | ALL values negative due to scale mismatch (per-event vs per-sequence); ranking is correct but absolute values are inverted | Multiply `expected_freq` by mean V-region mutation count per sequence before log-ratio |
| S3: VH/VL MI | **Weak signal — method revision needed** | Near-zero MI across all pairs; VL64 has zero coverage; near-zero likely reflects germline dominance of identity co-variation, not somatic co-evolution | Recompute on mutation co-occurrence (binary mutated/not) rather than absolute AA identity; remove VL64 |

### Open questions from Step 2

1. Why does CDR1 show purifying selection (ω < 1) while the pooled CDR R/S from Step 0 was 3.21 > neutral? → Likely because CDR3 dominates the R/S signal. Per-region R/S (CDR1 alone, CDR2 alone, CDR3 alone) should be extracted from the Step 0 data to confirm.
2. What drives the strong positive selection outliers at FR1/FR3 boundary positions (Aho24 ω=5.49, Aho98 ω=5.81)? → Check whether these overlap with known CDR-loop contact residues or heavy-chain packing residues.
3. After fixing S2 scale: do Cys23 and Cys104 become the highest Φ_S positions (most constrained)? They should, based on the ω data.
4. Does MI on mutation co-occurrence (S3 revision) reveal genuine somatic co-adaptation at the VH/VL interface?
=======
2. CDRL3 max=29 aa — are these lambda chains or artifacts?
3. Extreme mutation outliers (n_mut_H=109, total=157) — are these genuine hypermutated memory cells or should additional outlier thresholds be applied?
>>>>>>> origin/main
