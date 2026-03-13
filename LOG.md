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

1. **Why 49 duplicate seq_names from iggnition?** Are these truly duplicated input sequences, or a processing artefact? Worth checking if duplicates appear in the raw parquet too.
2. **Why 60% naive by bio label?** Experimental design likely explains this — confirm whether samples include explicitly sorted naive B cells.
3. **Extreme outliers** (max VH = 120 codons, max total = 157) — check post-QC whether stop codon filter removes these or if additional curation is needed.
4. **692,383 bio-naive but NOT comp-naive** — are these IgD+ cells (IgD not captured by `c_gene:0 == 'IGHM'` filter)? Clarify whether the dataset contains IgD-bearing sequences and how to classify them.
