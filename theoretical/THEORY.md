# THEORY.md — Mathematical Expression of Antibody Maturation

> A constrained multi-objective optimization framework for somatic hypermutation and germinal center selection, with static (Lagrangian/Pareto) and dynamic (Hamiltonian) formulations.

---

## 0. Notation and Definitions

### 0.1 Sequence Space

An antibody is a paired heavy + light chain. We represent it as a vector of amino acid identities at IMGT-numbered positions:

$$\mathbf{x} = (x_H, x_L) \in \Sigma^{N_H} \times \Sigma^{N_L}$$

where $\Sigma = \{A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y\}$ is the amino acid alphabet, and $N_H$, $N_L$ are the number of IMGT-aligned positions in heavy and light chains respectively. All positions are aligned and numbered using `iggnition` (github.com/bnemoz/iggnition), which provides consistent IMGT-numbered nucleotide-level alignment across all sequences.

We use the shorthand $\mathbf{x} \in \Omega$ for the full paired sequence space $\Omega = \Sigma^{N_H + N_L}$.

### 0.2 Germline Reference

Each sequence $\mathbf{x}$ has an inferred germline reference $\mathbf{g}(\mathbf{x}) \in \mathcal{G}$, determined by V(D)J gene assignment. The germline set $\mathcal{G}$ is the collection of all functional IGHV, IGHD, IGHJ, IGKV, IGKJ, IGLV, IGLJ gene segments.

The **mutation state** of $\mathbf{x}$ relative to its germline is the set:

$$\Delta(\mathbf{x}) = \{i : x_i \neq g_i(\mathbf{x})\}$$

and the **mutation count** is $d(\mathbf{x}) = |\Delta(\mathbf{x})|$.

### 0.3 Structural Regions

Each IMGT-numbered position $i$ belongs to a structural region $r(i) \in \{\text{FWR1}, \text{CDR1}, \text{FWR2}, \text{CDR2}, \text{FWR3}, \text{CDR3}, \text{FWR4}\}$. We use the indicator functions:

$$\mathbb{1}_{\text{CDR}}(i) = \begin{cases} 1 & \text{if } r(i) \in \{\text{CDR1}, \text{CDR2}, \text{CDR3}\} \\ 0 & \text{otherwise} \end{cases}$$

$$\mathbb{1}_{\text{FWR}}(i) = 1 - \mathbb{1}_{\text{CDR}}(i)$$

### 0.4 Available Metadata

Each sequence $\mathbf{x}$ in the dataset carries metadata:
- **Donor:** $\delta(\mathbf{x}) \in \mathcal{D}$ (donor identity)
- **Isotype:** $\iota(\mathbf{x}) \in \{\text{IgM}, \text{IgD}, \text{IgG1-4}, \text{IgA1-2}, \text{IgE}\}$
- **Compartment:** $c(\mathbf{x}) \in \{\text{naïve}, \text{memory}\}$
- **Clonal family:** $\mathcal{C}(\mathbf{x})$ (clonal lineage assignment, available for memory)

---

## Part I — The Static Picture

### 1. The Mutational Space $\mathcal{M}$

#### 1.1 Intrinsic Mutability

Somatic hypermutation (SHM) is not random. It is initiated by activation-induced cytidine deaminase (AID), which deaminates cytosine to uracil preferentially at WRC hotspot motifs (W = A/T, R = A/G). Downstream error-prone repair by DNA polymerase η introduces additional mutations preferentially at WA motifs. Both processes operate on the nucleotide-level sequence, but their effects are observed at the amino acid level after translation.

For each nucleotide position $i$ in germline $g$, the **intrinsic mutability** is:

$$\mu_i(g) = \sum_{b' \neq b_i} S_{5F}(c_i,\; b_i \to b')$$

where $b_i$ is the germline nucleotide at position $i$, $c_i = (b_{i-2}, b_{i-1}, b_i, b_{i+1}, b_{i+2})$ is the local 5-mer context, and $S_{5F}$ is the context-dependent mutability score (Yaari et al., *Front. Immunol.*, 2013). This captures both AID (WRC/RGYW) and Polη (WA/TW) targeting.

Since `iggnition` provides nucleotide-level aligned positions, we compute $\mu_i$ directly on the nucleotide sequence and then aggregate to the codon/amino acid level.

#### 1.2 Decomposition of Mutational Potential

The total mutational potential of germline $g$ decomposes into targeted and untargeted components:

$$V_{mut}(g) = V_{AID}(g) + V_{Pol\eta}(g) + V_{bg}(g)$$

where:

$$V_{AID}(g) = \sum_{i=1}^{L(g)} \mu_i^{AID}(g) = \sum_{i=1}^{L(g)} \sum_{b'} S_{5F}(c_i,\; b_i \to b') \cdot \mathbb{1}[c_i \in \text{WRC/RGYW}]$$

$$V_{Pol\eta}(g) = \sum_{i=1}^{L(g)} \mu_i^{Pol\eta}(g) = \sum_{i=1}^{L(g)} \sum_{b'} S_{5F}(c_i,\; b_i \to b') \cdot \mathbb{1}[c_i \in \text{WA/TW}]$$

$$V_{bg}(g) = \sum_{i=1}^{L(g)} \mu_i^{bg}(g) = \sum_{i=1}^{L(g)} \sum_{b'} S_{5F}(c_i,\; b_i \to b') \cdot \mathbb{1}[c_i \notin \text{WRC} \cup \text{WA}]$$

The regional mutational potential further decomposes by structural context:

$$V_{mut}^{CDR}(g) = \sum_{i\,:\, r(i) \in CDR} \mu_i(g), \qquad V_{mut}^{FWR}(g) = \sum_{i\,:\, r(i) \in FWR} \mu_i(g)$$

The ratio $R_{AID}(g) = V_{mut}^{CDR}(g)\, /\, V_{mut}^{FWR}(g)$ quantifies the **evolvability** of germline $g$ — the degree to which the SHM machinery is pre-targeted to the antigen-binding loops.

#### 1.3 The Feasible Mutational Space

After $t$ rounds of SHM (measured in cell divisions within the germinal center), the set of sequences reachable from germline $g$ is:

$$\mathcal{M}(g, t) = \left\{ \mathbf{x} \in \Omega : P(\mathbf{x} \mid g, t) > \epsilon \right\}$$

where $P(\mathbf{x} \mid g, t)$ is the probability of generating $\mathbf{x}$ from $g$ under neutral SHM (no selection), and $\epsilon$ is a threshold for practical reachability. The mutation count in the neutral model follows approximately:

$$d(\mathbf{x}) \mid g, t \;\sim\; \text{Poisson}\left(\lambda(g) \cdot t\right)$$

where $\lambda(g) = \sum_i \mu_i(g)$ is the total per-division mutation rate for germline $g$. Empirically, $\lambda \approx 10^{-3}$ per bp per division, yielding approximately one mutation per V-region per division.

**Crucially, $\mathcal{M}(g, t)$ defines the feasibility constraint — the "playing field" on which optimization occurs. It is not itself an objective.**

---

### 2. The Three Objective Functions

Germinal center (GC) selection acts on sequences within $\mathcal{M}$ to simultaneously satisfy three requirements. We formalize these as **penalty functions** to be minimized (lower is better):

#### 2.1 Structural Penalty $\Phi_S(\mathbf{x})$

The structural penalty measures how far $\mathbf{x}$ deviates from a foldable, stable immunoglobulin. Mutations that disrupt the hydrophobic core, the VH/VL interface, disulfide bonds, or canonical CDR loop conformations increase $\Phi_S$.

**Definition (language-model-based):**

$$\Phi_S(\mathbf{x}) = -\frac{1}{N} \sum_{i=1}^{N} \log P_{LM}(x_i \mid \mathbf{x}_{-i})$$

where $P_{LM}$ is the masked marginal probability from a protein language model (e.g., ESM-2). This pseudo-perplexity captures whether each residue is "expected" in its structural context.

**Definition (sequence-feature-based, for computation at scale):**

$$\Phi_S(\mathbf{x}) = w_1 \cdot \Phi_S^{core}(\mathbf{x}) + w_2 \cdot \Phi_S^{interface}(\mathbf{x}) + w_3 \cdot \Phi_S^{canonical}(\mathbf{x})$$

where:
- $\Phi_S^{core}(\mathbf{x})$: count of hydrophobic core positions (IMGT positions 6, 23, 41, 53, 78, 89, 104, 118 in VH; analogous in VL) mutated to non-conservative residues
- $\Phi_S^{interface}(\mathbf{x})$: deviation at VH/VL Vernier zone positions (IMGT positions 2, 47, 48, 67, 69, 71, 78, 93, 94 in VH) — requires paired data
- $\Phi_S^{canonical}(\mathbf{x})$: count of canonical CDR structure-determining positions (Chothia/Al-Lazikani classification) bearing non-canonical mutations

The weights $w_1, w_2, w_3$ are learned from data (Section 9).

#### 2.2 Affinity Deficit $\Phi_A(\mathbf{x})$

The affinity deficit measures the inverse of binding strength to cognate antigen. Since we lack antigen-binding measurements for the full repertoire, we use evolutionary proxies.

**Definition (selection-signal-based):**

$$\Phi_A(\mathbf{x}) = -\left[ \frac{R_{CDR}(\mathbf{x}) / S_{CDR}(\mathbf{x})}{(R/S)_{neutral}(g)} - 1 \right]$$

where $R_{CDR}$ and $S_{CDR}$ are the counts of replacement and silent mutations in CDR regions, and $(R/S)_{neutral}$ is the expected ratio under no selection (determined by codon structure and intrinsic mutability, typically $\approx$ 2.9 for V regions). A negative $\Phi_A$ (R/S above neutral expectation) indicates positive selection for affinity.

**Alternative proxy — CDR mutation enrichment:**

$$\Phi_A(\mathbf{x}) = -\log \frac{f_{CDR}^{obs}(\mathbf{x})}{f_{CDR}^{exp}(\mathbf{x})}$$

where $f_{CDR}^{obs}$ is the observed fraction of mutations in CDRs and $f_{CDR}^{exp}$ is the fraction expected from intrinsic mutability alone ($V_{mut}^{CDR} / V_{mut}^{total}$). Strong CDR enrichment indicates affinity-driven selection.

#### 2.3 Reactivity Risk $\Phi_R(\mathbf{x})$

The reactivity risk measures the propensity for polyreactivity or self-reactivity, which triggers negative selection in the GC and at peripheral tolerance checkpoints.

Since no functional polyreactivity data is available, we use a **sequence-based composite score** calibrated empirically (see Section 11):

$$\Phi_R(\mathbf{x}) = \alpha_R \cdot H_{H3}(\mathbf{x}) + \beta_R \cdot Q_{H3}^{+}(\mathbf{x}) + \gamma_R \cdot \ell_{H3}(\mathbf{x}) + \delta_R \cdot Y_{H3}(\mathbf{x})$$

where:
- $H_{H3}$: mean Kyte-Doolittle hydrophobicity of CDRH3 residues (Wardemann et al., *Science*, 2003: highly hydrophobic CDRH3 correlates with polyreactivity)
- $Q^+_{H3}$: net positive charge of CDRH3 at physiological pH (positively charged CDRH3 associate with nucleic acid reactivity)
- $\ell_{H3}$: CDRH3 amino acid length (longer CDRH3, especially >20 aa, correlates with polyreactivity)
- $Y_{H3}$: aromatic residue fraction in CDRH3 (aromatic residues contribute to non-specific hydrophobic contacts)

The coefficients $(\alpha_R, \beta_R, \gamma_R, \delta_R)$ are **calibrated from the data** using the naïve-to-memory transition as a natural experiment: the naïve compartment has not been filtered by GC selection, so the shift in $\Phi_R$ distribution from naïve to memory reveals the direction and magnitude of anti-polyreactivity selection.

---

### 3. The Constrained Optimization (Lagrangian Formulation)

Affinity maturation in the germinal center can now be stated as a constrained optimization. The B cell seeks to minimize the affinity deficit (maximize binding) while keeping structural disruption and self-reactivity below tolerance thresholds:

$$\min_{\mathbf{x} \in \mathcal{M}(g)} \Phi_A(\mathbf{x}) \quad \text{subject to} \quad \Phi_S(\mathbf{x}) \leq \Phi_S^{max}, \quad \Phi_R(\mathbf{x}) \leq \Phi_R^{max}$$

The **Lagrangian** is:

$$\boxed{\mathcal{L}(\mathbf{x}, \lambda_S, \lambda_R) = \Phi_A(\mathbf{x}) + \lambda_S \left[\Phi_S(\mathbf{x}) - \Phi_S^{max}\right] + \lambda_R \left[\Phi_R(\mathbf{x}) - \Phi_R^{max}\right]}$$

where $\lambda_S \geq 0$ and $\lambda_R \geq 0$ are the Lagrange multipliers — the **shadow prices** of structural integrity and non-reactivity.

**Biological interpretation of the multipliers:**
- $\lambda_S$ quantifies how much affinity improvement the system "pays" per unit of structural destabilization. Large $\lambda_S$ implies that structure is expensive — the system is conservative.
- $\lambda_R$ quantifies how much affinity improvement the system "pays" per unit of polyreactivity increase. Large $\lambda_R$ implies that tolerance checkpoints are stringent.

### 3.1 Karush-Kuhn-Tucker (KKT) Conditions

At the optimal mature antibody $\mathbf{x}^*$, the KKT conditions hold:

**Stationarity:**

$$\nabla_{\mathbf{x}} \Phi_A(\mathbf{x}^*) + \lambda_S^* \nabla_{\mathbf{x}} \Phi_S(\mathbf{x}^*) + \lambda_R^* \nabla_{\mathbf{x}} \Phi_R(\mathbf{x}^*) = \mathbf{0}$$

**Primal feasibility:**

$$\Phi_S(\mathbf{x}^*) \leq \Phi_S^{max}, \qquad \Phi_R(\mathbf{x}^*) \leq \Phi_R^{max}$$

**Dual feasibility:**

$$\lambda_S^* \geq 0, \qquad \lambda_R^* \geq 0$$

**Complementary slackness:**

$$\lambda_S^* \left[\Phi_S(\mathbf{x}^*) - \Phi_S^{max}\right] = 0, \qquad \lambda_R^* \left[\Phi_R(\mathbf{x}^*) - \Phi_R^{max}\right] = 0$$

Rearranging the **stationarity condition** yields the central equation:

$$\boxed{\nabla_{\mathbf{x}} \Phi_A = -\lambda_S \nabla_{\mathbf{x}} \Phi_S - \lambda_R \nabla_{\mathbf{x}} \Phi_R}$$

**Biological meaning:** At the mature optimum, the gradient of affinity improvement is exactly balanced by the gradients of structural and reactivity penalties, weighted by their respective prices. Every observed mutation in a mature antibody represents a local resolution of this three-way trade-off.

**Complementary slackness** yields an elegant classification:
- If $\lambda_S^* > 0$: the structural constraint is **active** (binding). The antibody has been pushed to its structural limit by affinity selection. This is the regime of highly mutated antibodies (e.g., HIV bnAbs with ~30% SHM).
- If $\lambda_S^* = 0$: structure is not a binding constraint. The antibody has structural slack. This is the regime of germline-close antibodies (e.g., m336 against MERS-CoV with a single VH mutation).
- If $\lambda_R^* > 0$: tolerance is active. The antibody has reached the polyreactivity boundary. This is the regime of antibodies like IGHV4-34-derived clones, which must "walk the line" between affinity and self-reactivity.
- If $\lambda_R^* = 0$: tolerance is not limiting. The antibody is far from autoreactive — self-reactivity is not constraining maturation.

---

### 4. Multi-Objective Optimization (Pareto Formulation)

Rather than treating affinity as the sole objective with structure and reactivity as constraints, we can equivalently cast all three as objectives in a multi-objective optimization:

$$\min_{\mathbf{x} \in \mathcal{M}(g)} \mathbf{f}(\mathbf{x}) = \left(\Phi_S(\mathbf{x}),\; \Phi_A(\mathbf{x}),\; \Phi_R(\mathbf{x})\right)$$

#### 4.1 Pareto Dominance

A sequence $\mathbf{x}$ **Pareto-dominates** $\mathbf{y}$ (written $\mathbf{x} \prec \mathbf{y}$) if:

$$\Phi_k(\mathbf{x}) \leq \Phi_k(\mathbf{y}) \;\; \forall k \in \{S, A, R\} \quad \text{and} \quad \exists \; j : \Phi_j(\mathbf{x}) < \Phi_j(\mathbf{y})$$

The **Pareto front** $\mathcal{P}$ is the set of all non-dominated sequences:

$$\mathcal{P}(g) = \left\{ \mathbf{x}^* \in \mathcal{M}(g) : \nexists\; \mathbf{x} \in \mathcal{M}(g) \text{ with } \mathbf{x} \prec \mathbf{x}^* \right\}$$

In the three-objective case, $\mathcal{P}$ is a **2-dimensional surface** (manifold) embedded in the 3D objective space $(\Phi_S, \Phi_A, \Phi_R)$.

#### 4.2 Weighted Scalarization

To trace the Pareto front, we define the family of scalarized problems parameterized by weights $\alpha, \beta \geq 0$, $\alpha + \beta \leq 1$:

$$\mathcal{L}_{\alpha,\beta}(\mathbf{x}) = \alpha\, \Phi_S(\mathbf{x}) + \beta\, \Phi_A(\mathbf{x}) + (1 - \alpha - \beta)\, \Phi_R(\mathbf{x})$$

Each $(\alpha, \beta)$ pair defines a direction on the weight simplex and selects a point on the Pareto front:

$$\mathbf{x}^*(\alpha, \beta) = \arg\min_{\mathbf{x} \in \mathcal{M}(g)} \mathcal{L}_{\alpha,\beta}(\mathbf{x})$$

The simplex $\{(\alpha, \beta) : \alpha, \beta \geq 0,\; \alpha + \beta \leq 1\}$ parameterizes the space of evolutionary strategies:

| $(\alpha, \beta)$ regime | Strategy | Biological archetype |
|---|---|---|
| $\alpha \to 1$ | Structure-first | Germline-close, broadly reactive (e.g., MERS m336) |
| $\beta \to 1$ | Affinity-first | Heavily mutated, high specificity (e.g., HIV VRC01) |
| $(1{-}\alpha{-}\beta) \to 1$ | Tolerance-first | Self-reactivity-constrained (e.g., IGHV4-34 clones) |
| $\alpha \approx \beta \approx 1/3$ | Balanced | Typical memory B cell response |

**Note on convexity:** Weighted scalarization only recovers points on the convex hull of the Pareto front. For potentially non-convex regions of the front, the $\epsilon$-constraint method is needed:

$$\min_\mathbf{x} \Phi_A(\mathbf{x}) \quad \text{s.t.} \quad \Phi_S(\mathbf{x}) \leq \epsilon_S, \;\; \Phi_R(\mathbf{x}) \leq \epsilon_R$$

Sweeping $(\epsilon_S, \epsilon_R)$ over a grid recovers the full front including non-convex portions.

---

## Part II — The Dynamic Picture

### 5. Lineage Trajectories in Objective Space

The static picture describes where mature antibodies end up. The dynamic picture describes how they get there. With clonal lineage data, each lineage is an ordered sequence of intermediate states:

$$\mathbf{x}_0 \to \mathbf{x}_1 \to \mathbf{x}_2 \to \cdots \to \mathbf{x}_T$$

where $\mathbf{x}_0 = \mathbf{g}$ is the germline (unmutated ancestor) and $\mathbf{x}_T$ is the observed mature antibody. Each step $\mathbf{x}_t \to \mathbf{x}_{t+1}$ represents one or more SHM events accepted (or drifted through) in the GC.

Each intermediate maps to a point in the 3D objective space:

$$\mathbf{f}(\mathbf{x}_t) = \left(\Phi_S(\mathbf{x}_t),\; \Phi_A(\mathbf{x}_t),\; \Phi_R(\mathbf{x}_t)\right)$$

The trajectory $\{\mathbf{f}(\mathbf{x}_t)\}_{t=0}^T$ traces a **path through the objective landscape**, from the germline origin toward the Pareto front. The geometry of these paths — their curvature, acceleration, and approach angle — encodes the dynamics of selection.

---

### 6. Hamiltonian Formulation

We embed the discrete lineage dynamics in a continuous framework. Let:

- $\mathbf{q}(t) \in \mathbb{R}^D$: position in a continuous embedding of sequence space at maturation time $t$
- $\mathbf{p}(t) = \dot{\mathbf{q}}(t)$: the "momentum" — the instantaneous velocity of evolution in the embedded space

The embedding can be constructed from ESM embeddings, from a PCA of mutation-state vectors, or directly in the 3D objective space $(\Phi_S, \Phi_A, \Phi_R)$.

#### 6.1 The Potential Energy Surface

The **potential energy** is the weighted penalty landscape:

$$U(\mathbf{q}) = \lambda_S\, \Phi_S(\mathbf{q}) + \lambda_A\, \Phi_A(\mathbf{q}) + \lambda_R\, \Phi_R(\mathbf{q})$$

This surface defines the "terrain" over which B cell lineages evolve. Germinal center selection drives lineages downhill on $U$ — toward lower total penalty.

#### 6.2 The Kinetic Energy

The **kinetic energy** represents the mutational activity — the capacity for sequence change:

$$T(\mathbf{p}) = \frac{1}{2} \mathbf{p}^T \mathbf{M}^{-1} \mathbf{p}$$

where $\mathbf{M}$ is the **mutability mass matrix** — a diagonal matrix whose entries encode the per-position resistance to mutation:

$$M_{ii} = \frac{1}{\mu_i(g) + \epsilon}$$

where $\mu_i(g)$ is the intrinsic mutability at position $i$ and $\epsilon$ is a small regularizer. Positions with high intrinsic mutability (AID hotspots in CDRs) have low "mass" (easy to move); positions with low mutability (FWR coldspots) have high "mass" (resistant to change).

#### 6.3 The Hamiltonian

$$\boxed{H(\mathbf{q}, \mathbf{p}) = \frac{1}{2}\mathbf{p}^T \mathbf{M}^{-1} \mathbf{p} + \lambda_S\, \Phi_S(\mathbf{q}) + \lambda_A\, \Phi_A(\mathbf{q}) + \lambda_R\, \Phi_R(\mathbf{q})}$$

**Hamilton's equations of motion:**

$$\dot{\mathbf{q}} = \frac{\partial H}{\partial \mathbf{p}} = \mathbf{M}^{-1} \mathbf{p}$$

$$\dot{\mathbf{p}} = -\frac{\partial H}{\partial \mathbf{q}} = -\lambda_S \nabla \Phi_S - \lambda_A \nabla \Phi_A - \lambda_R \nabla \Phi_R$$

**Biological meaning:**
- **First equation:** The rate of sequence evolution at each position is proportional to the momentum divided by the mutational resistance. Highly mutable positions evolve faster.
- **Second equation:** The force on the system (the acceleration of evolution) is the negative gradient of the total penalty. Selection pushes the trajectory toward regions of lower penalty.

#### 6.4 Stochastic Extension: Langevin Dynamics

SHM is inherently stochastic — not all generated mutations are beneficial, and drift operates alongside selection. Adding noise yields the **Langevin equation**:

$$\dot{\mathbf{q}} = -\mathbf{M}^{-1} \nabla U(\mathbf{q}) + \sigma(\mathbf{q})\, \boldsymbol{\eta}(t)$$

where $\boldsymbol{\eta}(t)$ is Gaussian white noise and $\sigma(\mathbf{q})$ is the position-dependent diffusion coefficient:

$$\sigma_i(\mathbf{q}) = \sqrt{2\, k_B T_{eff} \cdot \mu_i(g)}$$

Here, $T_{eff}$ is an effective "temperature" controlling the exploration-exploitation balance:
- **High $T_{eff}$**: broad exploration, low selection stringency (early GC, dark zone proliferation)
- **Low $T_{eff}$**: focused refinement, high selection stringency (late GC, tight Tfh competition)

The time-dependence of $T_{eff}(t)$ models the well-documented increase in selection stringency over the course of a GC reaction. An exponential cooling schedule $T_{eff}(t) = T_0 \cdot e^{-\kappa t}$ is a natural first model (analogous to simulated annealing).

#### 6.5 Conservation Laws and the Action Principle

In the deterministic limit ($T_{eff} \to 0$), the Hamiltonian $H$ is conserved along trajectories. This yields a testable prediction: **the total "energy" (kinetic + potential) should be approximately constant along observed lineage trajectories.** Departures from conservation indicate external perturbations — changing antigen, viral escape, re-entry into GCs.

The **action** along a trajectory is:

$$S[\mathbf{q}] = \int_0^T \left[ T(\dot{\mathbf{q}}) - U(\mathbf{q}) \right] dt$$

The principle of least action states that observed trajectories minimize $S$. By comparing the action of observed lineage trajectories to random walks of the same length in $\mathcal{M}$, we can test whether maturation follows near-optimal paths.

---

### 7. Connecting Static and Dynamic

The two formulations are related by the following correspondence:

| Static (Lagrangian/Pareto) | Dynamic (Hamiltonian) |
|---|---|
| Pareto front $\mathcal{P}$ | Terminal states of optimal trajectories |
| Lagrange multipliers $\lambda_S, \lambda_R$ | Weights in the potential $U(\mathbf{q})$ |
| Scalarization weights $(\alpha, \beta)$ | Direction of approach to the front |
| Feasible set $\mathcal{M}(g)$ | Initial condition $\mathbf{q}(0) = \mathbf{g}$ + mass matrix $\mathbf{M}$ |
| Complementary slackness | Active constraints at terminal time |

The Lagrangian tells us **where** mature antibodies end up. The Hamiltonian tells us **how** they got there. Together, they provide a complete description of antibody maturation as a physical process.

---

## Part III — Empirical Parameterization (Calculation Roadmap)

### 8. Calibrating $V_{mut}$ from Data

| ID | Calculation | Method | Output |
|---|---|---|---|
| V1 | Per-germline mutability profiles | Compute S5F mutability at every `iggnition`-aligned nucleotide position per germline | $\mu_i(g)$ matrix |
| V2 | AID hotspot density (CDR vs FWR) | Count WRC/RGYW and WA/TW motifs per region per germline | $R_{AID}(g)$ evolvability index |
| V3 | CDR3 length effect | Stratify memory sequences by CDRH3 length, compute mean $d(\mathbf{x})$ | Regression: $d$ vs $\ell_{H3}$ |
| V4 | Mutation accumulation by isotype | Plot $d(\mathbf{x})$ vs isotype class; fit Poisson model | Effective time $t(\iota)$ per isotype |

### 9. Calibrating $\Phi_S$ from Data

| ID | Calculation | Method | Output |
|---|---|---|---|
| S1 | Position-specific selection ($\omega_i$) | $\omega_i = (dN_i / E[dN_i]) / (dS_i / E[dS_i])$ with S5F as neutral baseline | Per-position $\omega$ map |
| S2 | Forbidden mutations | $\Phi_S(i) = \log[\mu_i(g) / (f_i^{obs}(g) + \epsilon)]$ | Structural filter profile |
| S3 | VH/VL interface co-variation | Mutual information at Vernier zone position pairs (paired data) | $MI(i_H, j_L)$ matrix |

### 10. Calibrating $\Phi_A$ from Data

| ID | Calculation | Method | Output |
|---|---|---|---|
| A1 | CDR replacement enrichment | $\Delta(R/S) = [R_{CDR}/S_{CDR}] / (R/S)_{neutral}$ per sequence | Selection intensity distribution |
| A2 | Convergent/public clonotypes | Find CDRH3 shared across $\geq 3$ donors; compare mutation profiles to private clonotypes | Public vs private mutation spectra |
| A3 | Isotype-stratified affinity proxy | Within clonal families with IgM + IgG members: identify IgG-enriched mutations | Candidate affinity mutations by position |

### 11. Calibrating $\Phi_R$ from Data

| ID | Calculation | Method | Output |
|---|---|---|---|
| R1 | Naïve-to-memory shift | Logistic regression: $P(\text{memory} \mid H, Q^+, \ell, Y)$; fitted coefficients are empirical weights | $(\hat{\alpha}_R, \hat{\beta}_R, \hat{\gamma}_R, \hat{\delta}_R)$ |
| R2 | Polyreactivity risk by germline | Compute $\langle \Phi_R \rangle$ per IGHV gene | Germline-level reactivity propensity |

### 12. Estimating the Lagrange Multipliers

| ID | Calculation | Method | Output |
|---|---|---|---|
| L1 | Gradient decomposition | For each accepted mutation: compute $\Delta\Phi_S$, $\Delta\Phi_A$, $\Delta\Phi_R$; fit $\Delta\Phi_A \approx -\hat{\lambda}_S \Delta\Phi_S - \hat{\lambda}_R \Delta\Phi_R$ | $\hat{\lambda}_S$, $\hat{\lambda}_R$ (global) |
| L2 | Per-germline $\lambda$ | Repeat L1 per IGHV gene | $\lambda_S(g)$, $\lambda_R(g)$ per germline |
| L3 | Per-donor $\lambda$ | Repeat L1 per donor | $\lambda_S(\delta)$, $\lambda_R(\delta)$ per individual |

### 13. Mapping the Pareto Front

| ID | Calculation | Method | Output |
|---|---|---|---|
| P1 | Empirical 3D front | Compute $(\Phi_S, \Phi_A, \Phi_R)$ for all 2.5M memory sequences; non-dominated sorting | Pareto surface $\mathcal{P}$ |
| P2 | Front by germline | Separate fronts per IGHV gene | Germline-specific front shapes |
| P3 | Front by isotype | Separate fronts per isotype class | Isotype-specific front locations |
| P4 | Front by donor | Separate fronts per individual | Inter-individual variation |

### 14. Dynamic Analysis: Lineage Trajectories

| ID | Calculation | Method | Output |
|---|---|---|---|
| D1 | Trajectory mapping | Map each node in lineage trees to $(\Phi_S, \Phi_A, \Phi_R)$ | 3D lineage trajectories |
| D2 | Velocity and acceleration | $\mathbf{v} = \Delta\mathbf{f}/\Delta d$; $\mathbf{a} = \Delta\mathbf{v}/\Delta d$ per edge | Deceleration toward front? |
| D3 | Action comparison | Discrete action $S$ for observed vs random trajectories | Optimality of maturation paths |
| D4 | Effective temperature | Variance of trajectory velocities at different maturation depths | $T_{eff}$ vs maturation stage |

---

## Part IV — Synthesis and Biomarkers

### 15. The $(\alpha, \beta)$ Strategy Biomarker

Each mature antibody's position on the Pareto front can be characterized by its **strategy vector**:

$$\hat{\alpha}(\mathbf{x}) = \frac{\Phi_S(\mathbf{x})}{\Phi_S(\mathbf{x}) + \Phi_A(\mathbf{x}) + \Phi_R(\mathbf{x})}, \quad \hat{\beta}(\mathbf{x}) = \frac{\Phi_A(\mathbf{x})}{\Phi_S(\mathbf{x}) + \Phi_A(\mathbf{x}) + \Phi_R(\mathbf{x})}$$

The distribution of $(\hat{\alpha}, \hat{\beta})$ across a repertoire defines the **maturation strategy profile**.

**Testable hypotheses:**
1. Effective vaccine responses produce a characteristic $(\alpha, \beta)$ distribution.
2. Autoimmune disease is characterized by anomalous $\lambda_R$ (insufficient reactivity filtering).
3. Aging shifts/contracts the Pareto front (loss of strategy diversity).
4. Broadly neutralizing responses require extreme $\alpha$ tolerance (high structural slack for extensive SHM).

### 16. The $\lambda$-Signature

The per-germline Lagrange multiplier profile $(\lambda_S(g), \lambda_R(g))$ constitutes a **functional fingerprint** of each germline gene — encoding how much structural and reactivity "tax" it imposes on maturation. This signature:
- Predicts which germlines are most amenable to deep maturation (low $\lambda_S$, low $\lambda_R$)
- Identifies germlines requiring compensatory mutations (high $\lambda_S$)
- Flags germlines with inherent autoreactivity risk (high baseline $\Phi_R$, requiring high $\lambda_R$)

---

## References (Key)

- Yaari et al. (2013). *Models of somatic hypermutation targeting and substitution based on synonymous mutations from high-throughput immunoglobulin sequencing data.* Front. Immunol. 4:358.
- Wardemann et al. (2003). *Predominant autoantibody production by early human B cell precursors.* Science 301:1374-1377.
- Julian et al. (2017). *Efficient affinity maturation of antibody variable domains requires co-selection of compensatory mutations to maintain thermodynamic stability.* Sci. Rep. 7:45259.
- Fernández-Quintero et al. (2020). *Local and global rigidification upon antibody affinity maturation.* Front. Mol. Biosci. 7:182.
- Wei et al. (2015). *Overlapping hotspots in CDRs are critical sites for V region diversification.* PNAS 112:E728-E737.
- Hoehn et al. (2022). *Deep learning model of somatic hypermutation reveals importance of sequence context beyond hotspot targeting.* iScience 25:103668.
- Amitai et al. (2022). *Affinity maturation for an optimal balance between long-term immune coverage and short-term resource constraints.* PNAS 119:e2113512119.
- Mishra & Mariuzza (2018). *Insights into the structural basis of antibody affinity maturation from next-generation sequencing.* Front. Immunol. 9:117.
