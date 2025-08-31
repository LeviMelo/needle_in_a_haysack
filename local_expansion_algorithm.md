Alright—here’s the **one-stop, everything-in-one-place** write-up of our framework, with a clean skeleton, where each concept fits, what’s optional vs. core, and how the three conductances (internal / external / global) change the picture.

I’ll structure this as:

1. Vocabulary & objects (U, C, F, H, degrees)
2. Core signals we compute (per node & per state)
3. Gates (what admits a node) and why
4. Evaluation mode (synchronous vs. beam/batch) and synergy
5. Conductance—local deltas, sweep curve (diagnostic), and the 3 conductances (internal/external/global)
6. Megahubs & gateways (definitions and handling)
7. Diagnostics & stopping (objective, parameter-light)
8. Optional choices (edge-strength families, queueing, H\_core vs H\_full, etc.)
9. Edge cases & robustness
10. Skeleton (pseudocode-level plan you can implement)
11. What changed after introducing the 3 conductances (how it refines/automates the earlier plan)
12. Ultra-brief recap (TL;DR)

No new ideas here—just synthesis and clear placement of everything we’ve already discussed.

---

# 1) Vocabulary & Objects

* **U**: the **universe** graph (e.g., all of PubMed as a citation network). We do **not** materialize U, but for **any node** we *can* fetch the **full 1-hop ID lists** (references + citers), hence its **true total degree** $\deg_{\text{tot}}(v)$.

* **C**: current accepted cluster (starts at seeds $S$).

* **F**: **frontier** = all 1-hop neighbors (refs or citers) of any node in $C$ that are **not** in $C$, after filters (e.g., recency, hub containment).

* **H**: induced graph on $C \cup F$. We use **H\_full by default**:
  includes **all edges among $C\cup F$** (C–C, C–F, F–F). (H\_core—only edges that touch C—is a fallback for scale; with your data, H\_full is feasible and better.)

* **Undirected vs. directed**: for **cohesion/community** metrics (cut/condutância, blocks, articulation), we use an **undirected** view; direction is only used to build ref/citer sets that feed **structural similarity** features (e.g., triangles, AA/RA).

* **Degrees**:

  * $\deg_H(v)$: degree of $v$ inside H.
  * $\deg_{\text{tot}}(v)$: degree of $v$ in U (exact: from global 1-hop ID lists).
  * **Edges leaving H** from $v$: $\deg_{\text{outU}}(v)=\deg_{\text{tot}}(v)-\deg_H(v)$ (exact).

* **Internal/outside split relative to $C$** for a frontier node $v$:

  * $N_C(v) = N_H(v)\cap C$.
  * **Qualified internal degree** $d^{(\star)}_{\text{in}}(v;C)$: sum of **strong** ties $v$–$u$ for $u\in C$ (see §2.1).
  * **Raw out** inside H: $d_{\text{out}}(v;C)=\deg_H(v)-d^{(\star)}_{\text{in}}(v;C)$.
  * **Qualified out**: $d^{(\text{qual})}_{\text{out}}(v;C)$ (downweights out-edges into fringe nodes already attracted to C; §2.2).

---

# 2) Core Signals We Compute

## 2.1 Strong ties (edge strength $w^\star$) and $d^{(\star)}_{\text{in}}$

You choose **one** (or a mix) to define the **strength** of $v$–$u$ when $u\in C$:

* **Triangles in C** (closure):
  $T_{C\triangle}(v)$ = number of **C-edges** $(u,w)\in E(C)$ such that $v$ connects to both $u$ and $w$.
  Triangle **density** $\text{tri\_dens}(v)=\frac{T_{C\triangle}(v)}{1+|N_C(v)|}$.

* **CN$_C$** / **AA$_C$** / **RA$_C$** (all restricted to C and H):
  Common-neighbors counts for $v$ with each $u\in C\cap N(v)$; AA/RA **downweight** locally popular C-nodes (IDF-like, but purely structural).

* **EigAff-edge** (optional leash to C’s nucleus):
  Compute centrality (or simply degree) **on C** only; set $w^\star(v,u)=\xi(u)$ or $\deg_C(u)$.

Then assemble **qualified internal degree**:

$$
d^{(\star)}_{\text{in}}(v;C) = \sum_{u\in C\cap N(v)} \mathbb{1}\{w^\star(v,u)\ge\tau\}\cdot w^\star(v,u),
$$

where $\tau$ is a **local threshold** (e.g., wave-wise quantile 0.6–0.7 of the positive $w^\star$ values).

> **When to use which**:
> • If your data are dense and cliquish, **triangles** or **CN/AA/RA** shine.
> • If C has a clear “core”, **EigAff-edge** helps anchor to it (or use **degree-in-C** as a cheap proxy).
> • A **mix** (e.g., $0.7\,\xi(u)+0.3\,\text{AA}_C$) often works best.

## 2.2 Qualified out-degree $d^{(\text{qual})}_{\text{out}}$

Not all out-edges are equally “leaky”: penalize less the out-edges that point to F-nodes already gravitating to C.

* **Anchor score** for any $u\notin C$:
  $a(u)=\min\!\Big\{1,\ \kappa\cdot\dfrac{d^{(\star)}_{\text{in}}(u;C)}{\deg_{\text{tot}}(u)}\Big\}$, with $\kappa\in[1,3]$; cap at 1.

* **Qualified out** for $v$:

$$
d^{(\text{qual})}_{\text{out}}(v;C) = \sum_{u\in N(v)\setminus C} \Big(1-a(u)\Big)\cdot w_{\text{edge}}(v,u),
$$

with $w_{\text{edge}}$ = 1 (unweighted) or a simple structural weight.

This feeds $\Delta\phi$ (next) and reduces punishment for edges into the **same local basin**.

## 2.3 Per-candidate conductance delta (admission score)

Given current $\text{cut}(C)$ and $\mathrm{vol}(C)$ in H (undirected):

$$
\text{cut}' = \text{cut}(C) - d^{(\star)}_{\text{in}}(v;C) + d^{(\text{qual})}_{\text{out}}(v;C),
\quad
\mathrm{vol}' = \mathrm{vol}(C) + \deg_H(v),
$$

$$
\Delta\phi(v\mid C)= \frac{\text{cut}'}{\mathrm{vol}'} - \frac{\text{cut}(C)}{\mathrm{vol}(C)}.
$$

We require $\Delta\phi \le \epsilon_\phi$ (non-worsening; $\epsilon_\phi$ tiny).

## 2.4 “Hidden outside” pressure (U-aware)

Because you know global 1-hop IDs:

* **Exact** $\deg_{\text{outU}}(v)=\deg_{\text{tot}}(v)-\deg_H(v)$.
* **ECR (Excess-Core Ratio)** gate input:

  $$
  \text{ECR}(v) \;=\; \frac{d^{(\star)}_{\text{in}}(v;C)}{d_{\text{out}}(v;C) + \lambda\,\deg_{\text{outU}}(v) + \epsilon},
  $$

  with $\lambda\in[0.7,1]$. We gate on $\text{ECR}(v)\ge \tau_{\text{ECR}}$.

This kills “small local out, huge global out” boundary-nodes.

---

# 3) Gates (the admission checklist)

A candidate $v\in F$ is **admitted** if it passes **all**:

1. **k-in online**: $d^{(\star)}_{\text{in}}(v;C)\ge k_{\text{in}}$ (typically 2; raise to 3 if drift risk).
2. **Specificity/closure**: e.g., $\text{tri\_dens}(v)\ge \tau_\triangle$ (or AA/RA gate).
3. **ECR (U-aware)**: $\text{ECR}(v)\ge \tau_{\text{ECR}}$.
4. **Boundary tightening** (optional single-step check): $-d^{(\star)}_{\text{in}} + d_{\text{out}} + \deg_{\text{outU}}(v) \le 0$.
5. **Conductance non-worsening**: $\Delta\phi(v\mid C)\le \epsilon_\phi$, using **qualified** out and including $\deg_{\text{outU}}$ in the cut/volume update if you want a pessimistic guard.
6. **Recency** (if used): year within $\text{median}(C)\pm \Delta$.

> These gates are **orthogonal**: (1)/(2) assert **tight local anchoring**; (3)/(4) block **global leakage**; (5) enforces **community shape**; (6) keeps temporal coherence.

---

# 4) Evaluation Mode and Synergy

Two modes:

* **Synchronous (queue-less, recommended)**: in each wave, evaluate **all** $v\in F$ *against the same fixed $C$*, admit all that pass, commit them together. Order-free; misses only “same-round synergies,” which are tiny in practice.

* **Beam/batch (order-aware with synergy)**: build a beam $Q$ (size $b$) by a priority (e.g., $\alpha d^{(\star)}_{\text{in}} + \beta \text{tri\_dens} - \gamma d_{\text{out}} - \delta \deg_{\text{outU}}$).
  Select a subset $B\subseteq Q$ maximizing a score **with pairwise bonus** $\mu\sum_{v<w\in B}\mathbf{1}\{(v,w)\in E(F)\}$, subject to the **per-node gates**. Then commit $B$. This captures the fact that admitting both ends of an F–F edge shrinks the boundary by an extra **2**.

Either way, we **commit at end** of the wave (keeps decisions reproducible and unbiased within the wave).

---

# 5) Conductance: deltas, sweep curve (diagnostic), and the 3 conductances

## 5.1 Per-node (or per-batch) conductance **delta**

Already covered: $\Delta\phi(v\mid C)\le \epsilon_\phi$ using H\_full and qualified out. For batches $B$, add the pairwise $-2$ corrections for F–F edges admitted together.

## 5.2 **Sweep curve** (diagnostic only)

Pick **one** ranking (e.g., PPR on H or your priority). Form **prefixes** $C_k=S\cup\{v_1,\ldots,v_k\}$, plot $\phi(C_k)$ vs $k$.

* **One deep minimum** → one basin; stopping around that size is natural.
* **Two minima** (dip–rise–dip) → two basins linked by a neck; tighten gates so you stop at the **first** minimum.

## 5.3 The **three conductances** (U-aware state metrics)

* **Internal** $\phi_{C\to F} = \dfrac{|E(C,F)|}{\mathrm{vol}_H(C)}$
  “How leaky is C into F (inside H)?”

* **External** $\phi_{F\to U} = \dfrac{|E(F,U\setminus H)|}{\mathrm{vol}_U(F)} = \dfrac{\sum_{v\in F}\deg_{\text{outU}}(v)}{\sum_{v\in F}\deg_{\text{tot}}(v)}$
  “How strongly is F glued to the outside universe we haven’t materialized?”

* **Global** $\phi_{H\to U} = \dfrac{|E(H,U\setminus H)|}{\mathrm{vol}_U(H)}$
  “How sealed is the *whole* neighborhood we’re working in?”

A handy **risk index**:

$$
\mathcal{R} = \frac{\phi_{C\to F}}{\phi_{F\to U} + \varepsilon}.
$$

If $\mathcal{R}$ is **small**, the frontier is tightly glued outside → **auto-tighten** gates (raise $k_{\text{in}}$, $\tau$, $\tau_{\text{ECR}}$, shrink recency). If **large**, you can be a bit more permissive.

> **Did these change our approach?**
> Yes—**they reduce hand-tuned constants**. Gates become **data-driven**: their thresholds adapt from these conductances and observed distributions each wave. The decision logic (k-in, ECR, $\Delta\phi$, etc.) stays the same; how strict it is becomes **self-calibrated** by U-aware signals.

---

# 6) Megahubs & Gateways

* **Megahub (degree/tilt phenomenon)**:
  $\deg_{\text{tot}}(v)$ extremely high (e.g., top 1%) **or** tilt $t(v)=\dfrac{\deg_{\text{tot}}(v)-d^{(\star)}_{\text{in}}(v;C)}{d^{(\star)}_{\text{in}}(v;C)+1}$ extreme.
  Policy: never auto-expand; must pass **all** gates; ECR and qualified $\Delta\phi$ usually reject it early.

* **Gateway (boundary/neck phenomenon)**:
  (i) Inside $C$: nodes that are **articulation points** whose removal splits $C$ into ≥2 **large** pieces (or belong to multiple large biconnected blocks).
  (ii) On frontier: **gateway propensity**
  $g(v)=\dfrac{d^{(\text{qual})}_{\text{out}}(v;C)-d^{(\star)}_{\text{in}}(v;C)}{\deg_{\text{tot}}(v)}$.
  If rejected with large $g(v)$ or large positive $\Delta\phi$, **log as exit** (catalog where expansion would leave the theme).

* **Block-cut (Tarjan)** on $C$ each wave: count articulation points; track sizes of the two largest blocks $B_1, B_2$ and ratio $|B_2|/|B_1|$. Rising $|B_2|/|B_1|$ → two-nuclei risk → **tighten** gates next wave.

---

# 7) Diagnostics & Stopping

Per wave, record:

* $\phi_{C\to F},\, \phi_{F\to U},\, \phi_{H\to U},\, \mathcal{R}$.
* Acceptance rate; best $\Delta\phi$.
* **Internalization ratio** $\rho$ = fraction of edges incident to newly accepted nodes that land inside the **previous** $C$ (qualified variant: downweight out-edges to well-anchored fringe).
* Block health: # articulation points; $|B_2|/|B_1|$.
* Sweep minima positions (optional diagnostic).

**Stop** after $w$ consecutive waves (e.g., 2–3) where all hold:

1. No candidate/batch passes the pessimistic $\Delta\phi$ gate.
2. Boundary won’t tighten: $\min_v\{-d^{(\star)}_{\text{in}} + d_{\text{out}} + \deg_{\text{outU}}(v)\}\ge 0$.
3. Specificity spent: top quantiles of $\text{tri\_dens}$ (or AA/RA) in the beam fall below $\tau$.
4. External pull remains strong (low $\phi_{F\to U}$, small $\mathcal{R}$).
5. Block health not improving.

These are **objective**, no “frontier entanglement” heuristics.

---

# 8) Optional Choices (and how they swap)

* **Edge strength $w^\star$** (choose one or mix): triangles / CN / AA / RA / EigAff-edge (or degree-in-C).
* **Evaluation**: synchronous (queue-less) **or** beam/batch with synergy.
* **Priority (if using beams)**: $\alpha d^{(\star)}_{\text{in}} + \beta\text{tri\_dens} - \gamma d_{\text{out}} - \delta \deg_{\text{outU}}$ **or** PPR on H (ranking only).
* **H variant**: H\_full (default) **or** H\_core (fallback) with U-aware penalties to compensate missing F–F edges.
* **Recency**: on/off; window $\Delta$ auto-tuned from acceptance history.
* **Parameterized vs. data-driven thresholds**: you can start with fixed $k_{\text{in}}$, $\tau$, $\tau_{\text{ECR}}$, $\epsilon_\phi$, then switch to **auto-tuning** from quantiles and $\mathcal{R}$.

---

# 9) Edge Cases & Robustness

* **Tiny seed set S**: raise $k_{\text{in}}$ to 2–3; use stricter $\tau$; rely more on triangles/CN (or degree-in-C).
* **Nodes with giant ref lists**: internal-ratio effect is captured by ECR; also $\Delta\phi$ punishes large $\deg_H(v)$ with small $d^{(\star)}_{\text{in}}$.
* **Sparse areas (few triangles)**: prioritize EigAff-edge (or degree-in-C) and the ECR/Δφ gates.
* **“Two nuclei” emerging**: detected by sweep minima & blocks; raise $k_{\text{in}}$, $\tau$, shrink $\Delta$; optionally run a **sweep-cut** inside $C$ to trim back to the first minimum.

---

# 10) Skeleton (pseudocode-level)

**Config / knobs** (can be auto-tuned per wave):

* `k_in`, `tau_strong`, `tau_ECR`, `epsilon_phi`, `recency_window`, `kappa_anchor`, `lambda_outU`, `beam_size`, `mu_synergy`.

**Per wave:**

1. **Frontier & H**

   ```
   F = neighbors_1hop(C)  // refs OR citers
   F = apply_recency(F, recency_window)
   H = induced_full_graph(C ∪ F)  // undirected for cohesion
   ```

2. **Per-candidate features (for all v ∈ F)**

   ```
   deg_tot[v]      // from global 1-hop IDs
   deg_H[v]        // in H
   deg_outU[v] = deg_tot[v] - deg_H[v]

   // Strong ties
   for u in N(v) ∩ C:
       w_star[v,u] = edge_strength(v,u)  // triangles/CN/AA/RA/EigAff-edge or mix
   tau_strong = quantile({w_star>0}, q≈0.6–0.7)
   d_in_star[v] = sum_{u∈C∩N(v)} 1{w_star≥tau_strong} * w_star[v,u]

   // Qualified out
   for u in N(v) \ C:
       a(u) = min(1, kappa_anchor * d_in_star[u]/deg_tot[u])   // if u∈F; else a(u)=0
   d_out_qual[v] = sum_{u∈N(v)\C} (1 - a(u)) * w_edge(v,u)     // w_edge=1 is fine

   // Local out inside H (raw)
   d_out[v] = deg_H[v] - d_in_star[v]

   // ECR
   ECR[v] = d_in_star[v] / ( d_out[v] + lambda_outU * deg_outU[v] + eps )
   ```

3. **State metrics (U-aware & topology health)**

   ```
   phi_CtoF = |E(C,F)| / vol_H(C)
   phi_FtoU = (sum_{x∈F} deg_outU[x]) / (sum_{x∈F} deg_tot[x])
   phi_HtoU = (sum_{x∈H} deg_outU[x]) / (sum_{x∈H} deg_tot[x])
   R = phi_CtoF / (phi_FtoU + eps)

   // Block-cut diagnostics on C (undirected)
   (articulation_pts, blocks) = tarjan_biconnected(C)
   B1,B2 = sizes_of_two_largest_blocks(blocks)
   ```

4. **Auto-tune gates (optional, data-driven)**

   ```
   k_in      = smallest k∈{2,3,4} that stabilizes phi_CtoF and reduces B2/B1 on beam
   tau_ECR   = median(ECR) among top quartile by tri_dens (bump if R small or B2/B1 big)
   epsilon_phi = tiny or 5th percentile of past improvements; clamp to 0 if R small
   recency_window = adjust to keep ≥70% of last accepted inside window
   ```

5. **Evaluation**

   * **Synchronous**: for all v∈F, compute Δφ(v|C) with `d_out_qual` (and include `deg_outU` if pessimistic).
     Admit v iff all gates pass:

     ```
     (d_in_star[v] ≥ k_in) AND
     (tri_dens[v] ≥ tau_triangle OR AA/RA gate) AND
     (ECR[v] ≥ tau_ECR) AND
     (Δφ[v] ≤ epsilon_phi) AND
     (recency ok)
     ```

     Commit all admitted at once.

   * **Or Beam/Batch with synergy**:
     Build priority; pick subset B ⊆ Q maximizing node score + μ \* (# F–F edges inside B) subject to gates; commit B.

6. **Stopping**
   If no admissions or the objective stop conditions (in §7) hold for w consecutive waves, **stop**.

7. **Logging & exits**
   For rejected nodes with high gateway propensity $g(v)=(d_out_qual - d_in_star)/deg_tot$ or large positive Δφ: log as **exits/gateways** with neighbor IDs.

---

# 11) What the 3 Conductances Changed (vs. initial plan)

**Before**: thresholds like $k_{\text{in}}$, $ \tau$ (strong-edge), $ \tau_{\text{ECR}}$, $\epsilon_\phi$ were **fixed** or tuned by hand; stopping relied partly on local notions (e.g., frontier patterns).

**After introducing $ \phi_{C\to F}, \phi_{F\to U}, \phi_{H\to U}$**:

* We **measure** how sealed C is inside H (**internal**), how glued F is to the universe (**external**), and how sealed H is (**global**).
* A simple risk index $ \mathcal{R} = \phi_{C\to F}/(\phi_{F\to U}+\varepsilon)$ now **drives strictness**: when external pull is strong (small $ \mathcal{R}$), the algorithm **auto-tightens**—raising $k_{\text{in}}$, strong-edge threshold, ECR threshold, shrinking recency, and clamping $\epsilon_\phi$ to zero.
* Stopping is **objective**: no candidate/batch can pass pessimistic Δφ, boundary won’t tighten, specificity spent, external pull high, blocks unhealthy.
* Net: **same decision logic**, **far fewer hand-picked constants**—thresholds and strictness **adapt** to the measured geometry of U around H.

---

# 12) TL;DR

* Work on **H\_full = induced(C ∪ F)** with all edges among them; keep decisions **1-hop** and **structural**.
* Admit a frontier node only if it (i) has **enough strong ties** to C (k-in, triangles/AA/RA/EigAff), (ii) passes **ECR** (anchors outweigh its local+global out), (iii) **doesn’t worsen conductance** (Δφ ≤ ε), and (iv) respects **recency**.
* Use **qualified d\_out** (downweight out-edges into already anchored fringe) and **U-aware penalties** ($\deg_{\text{outU}}$).
* **Synchronous** evaluation is clean and order-free; batch/beam adds a **synergy** bonus for F–F edges.
* Track **block-cut** (articulation points, biconnected components) in C and the **sweep curve** as diagnostics; tighten gates if a second basin appears.
* Compute **three conductances** (C→F, F→U, H→U) and a **risk index** $ \mathcal{R}$ to **auto-tune thresholds** and **decide when to stop**.
* Log **gateways** (exits) systematically.

That’s the whole system, compact and complete. If you want next, I can turn the skeleton into concrete function signatures with exact inputs/outputs and incremental updates for each metric.

[Original thread](https://chatgpt.com/share/68b36fbb-05d0-800a-b416-c0d503d4cb75)