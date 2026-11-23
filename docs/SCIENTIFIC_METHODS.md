# ITME v23.0 (`v2.3.0`) Scientific Methods & Equation Reference

This document details the mathematical formalism, information-theoretic metrics, thermodynamic energetics, toroidal spatial ecology, and planetary biogeochemical couplings implemented in `itme_evolution_sim` (`v2.3.0`), licensed under **GPL-3.0-or-later**.

All equations are typeset in plain-text ASCII and Unicode notation (no LaTeX macros) so every symbol maps directly to executable Rust code in [`src/model.rs`](../src/model.rs), [`src/organisms.rs`](../src/organisms.rs), [`src/spatial.rs`](../src/spatial.rs), [`src/regimes.rs`](../src/regimes.rs), and [`src/engine.rs`](../src/engine.rs).

---

## 1. Macro-Evolutionary State & Order Velocity

ITME separates raw sequence bloat from functional adaptation by normalizing emergent work-and-information complexity against conserved non-neutral structure and downside persistence stability:

```text
1. Biospheric / Ecosystem Order:
   O_pop(t) = O(t) * ln(1 + N_pop(t)) * M_div(t) * M_web(t) * M_ETI(t)

2. Lineage Order & Order Velocity:
   O(t)     = S_T_minus(t) * L_I_bar(t)
   v_O(t)   = sum_k((k - k_bar) * (O_k - O_bar)) / sum_k((k - k_bar)^2)

3. Fitness-Weighted Informational Leverage:
   L_I_bar(t) = sum_i(w_i * L_I_i) / sum_i(w_i)
   L_I_i      = C_i / max(S_min, S_eff_i)       [S_min = 2.0 bits]

4. Asymmetric Downside Stability:
   S_T_minus(t) = clamp((1.0 - H_comb(t)) * (1.0 - P_downside(t)), 0.02, 1.0)
   H_comb(t)    = clamp(H_thermo(t) + H_meltdown(t), 0.0, 0.94)
   P_downside   = (0.25 + N_losses) / (4.25 + N_pairs)
```

- **`M_div(t)`**: `1.0 + 0.16 * max(0, N_species - 1)` rewards speciation across reproductively isolated clades.
- **`M_web(t)`**: `1.0 + 0.18 * H_web` scales with Shannon trophic guild entropy `H_web` across Producers, Syntrophs, Predators/Parasites, and Holobionts.
- **`M_ETI(t)`**: `1.0 + sum_j(beta_j * n_j)` aggregates individuality increments across holobionts, serial plastids, reduced MROs, clonal germ-soma organisms, eusocial superorganisms, agro-colonies, cultural demes, verified archives, and formal AI symbionts.
- **`P_downside(t)`**: Laplace-smoothed fraction of adjacent tick pairs in `TelemetryWindow` where `L_I(t) - L_I(t - 1) < -epsilon` (`epsilon = 0.08` default).

---

## 2. Structural Compression & Adami Physical Complexity

To penalize tandem sequence repetition while isolating functional information conserved by selection (`src/model.rs`):

```text
5. Lempel-Ziv Raw Structure & Multi-Genome Aggregation:
   S_seq     = c_LZ76(seq) * (log2(A) + 0.35 * log2(L)) + E_wiring * (0.25 * log2(L))
   S_raw_i   = S_host + 0.65 * S_endo + 0.55 * S_nm + S_neural

6. Adami (2000) Physical Complexity & Effective Conserved Structure:
   H_pos     = -sum_{a in A} p_a(pos) * log2(p_a(pos))
   rho_Adami = clamp(1.0 - (1.0 / L_eval) * sum_pos(H_pos / log2(A)), 0.08, 1.0)
   lambda_N  = N_rep / (N_rep + 8.0)          [0.0 when N_rep <= 1]
   S_eff_i   = S_min + max(0.0, S_raw_i - S_min) * ((1.0 - lambda_N) * rho_phys + lambda_N * rho_Adami)
```

- **`compute_adami_conservation`** (`src/model.rs`): Directly implements Adami, Ofria, & Collier (2000) by normalizing per-locus Shannon site entropy `H_pos` against the full alphabet capacity `log2(A)` (`h_max = (alphabet_size as f64).log2()`), rather than clamping to `log2(N_present)`. When `sequences.len() < 2`, `alphabet_size < 2`, `max_len == 0`, or no locus has `present >= 2`, the function safely returns finite `1.0`.
- **`rho_phys`**: Single-molecule physical folding constraint active prior to Darwinian ensemble replication.

---

## 3. Eigen Error Threshold & Price Selection Decomposition

```text
7. Eigen Critical Fidelity Threshold (Eq. 11, src/model.rs::compute_critical_fidelity):
   D_eff   = max(4.5, 0.22 * S_raw)
   theta_f = clamp(1.0 - ln(max(1.25, sigma)) / D_eff, 0.72, 0.999)

8. Error-Catastrophe Hazard:
   H_meltdown = clamp(2.4 * (theta_f - f_blend), 0.18, 0.58)   [if f_blend < theta_f & (1 - rho) > 0.05]
              = 0.0                                            [otherwise]

9. Price Equation Selection-Transmission Identity:
   Delta_z_bar = Cov(w, z) / w_bar + E[w * Delta_z] / w_bar
```

- Raising the upper clamp of `compute_critical_fidelity` to `0.999` ensures large megabase/symbolic structures (`S_raw > ~214 bits` at `sigma ~ 10`) accurately require near-unity per-symbol fidelity (`theta_f -> 0.999`), aligning directly with Eigen's quasispecies error threshold `1 - ln(sigma) / D_eff`.
- `PriceDecomposition::from_parent_offspring_records` separates selective covariance (`Cov(w, z) / w_bar`) from transmission/mutational drift (`E[w * Delta_z] / w_bar`), separating Darwinian regimes `R(3.0-DAR)` from quasispecies error catastrophe `R(2.0-QS)`.

---

## 4. Landauer Maintenance Dissipation & Net Exergy Transduction

Living and computational agents operate as open dissipative structures (`src/organisms.rs`, `src/model.rs`):

```text
10. Exergy Harvest, Landauer Heat & Starvation Hazard:
    W_harvest  = W_primary + W_syntrophic + J_trophic
    Q_Landauer = k_Landauer * (S_host + S_endo + S_nm + S_klepto + S_coord_import + 2.8 * norm1(m))
    Phi_net    = W_harvest - Q_Landauer
    Phi_exergy = clamp(0.76 + 0.48 * min(1.6, Phi_net), 0.65, 1.72)   [if Phi_net >= -0.02]
               = clamp(0.28 * I_int, 0.008, 0.32)                     [if Phi_net < -0.02]
    H_thermo   = clamp(0.42 + 0.55 * (-Phi_net), 0.35, 0.92)          [if Phi_net < -0.02]
               = clamp((0.042 / (0.55 + Phi_net)) / (1 + 0.65 * N), 0.008, 0.065)
```

- When `Q_Landauer` exceeds `W_harvest` (`Phi_net < -0.02`), thermodynamic starvation drives emergent transition into `R(0.5) Entropic Decay`.
- Functional complexity combines structural scaffold scaling, causal epistasis, exergy margin, 8-epoch transition synergy `Omega_synergy`, and lifetime plasticity:

```text
11. Emergent Functional Complexity:
    C_i = max(0.08, 4.25 * (S_eff_i ^ 0.62) * (0.38 + 1.45 * E_epistatic) * Phi_exergy * M_eng * Omega_synergy * I_int * (1.0 + 0.72 * m_bar * pi_mem))
```

---

## 5. Multi-Channel Heritable Mutual Information

ITME tracks up to five cumulative inheritance channels (`src/model.rs::compute_heritable_mutual_info`, `src/engine.rs`):

```text
12. Channel Mutual Information:
    H_cond(f, A)  = -(f * log2(f) + (1.0 - f) * log2((1.0 - f) / (A - 1.0)))
    I_chan(f,A,B) = B * clamp(1.0 - H_cond(f, A) / log2(A), 0.0, 1.0)
    I_herit_total = I_chan(f_gen, A, S_raw) + I_chan(f_cult, A, B_cult) + I_chan(f_arch, A, B_arch) + I_chan(f_proof, A, B_proof)
```

---

## 6. Toroidal 2D Spatial Reaction-Diffusion (`src/spatial.rs`)

Population ecology runs on a `16 x 16` (`DEFAULT_GRID_DIM = 16`) toroidal lattice tracking three coupled scalar fields: primary exergy `E(x, y)`, diffusible cross-feeding byproduct `B(x, y)`, and persistent niche capital `N(x, y)`.

- **Unified Periodic Indexing (`toroidal_stencil`)**: Both `SpatialGrid::idx(x, y)` and `SpatialGrid::step_reaction_diffusion(...)` (as well as `best_guild_neighbor`) wrap indices modulo `(width, height)`:
  - `left  = ((x + width - 1) % width, y)`
  - `right = ((x + 1) % width, y)`
  - `up    = (x, (y + height - 1) % height)`
  - `down  = (x, (y + 1) % height)`
- **Exact Discrete Laplacian Conservation**: Because the periodic 5-point stencil `lap(c) = left + right + up + down - 4 * center` sums strictly to zero across the torus, pure diffusion + linear decay conserves total mass (`sum_after = sum_before * (1 - decay)`) across boundaries without artificial edge reflections.
- **Non-Finite Input Guard**: `step_reaction_diffusion` checks `bath_exergy_flux.is_finite()`, `exergy_diffusion.is_finite()`, `niche_decay.is_finite()`, and `niche_diffusion.is_finite()` up front so `NaN` or `±Inf` inputs cannot corrupt grid fields.
- **Wright's `F_ST` & Speciation (`analyze_spatial_speciation`)**: Clusters lineages by Hamming sequence divergence (`<= theta_spec`) across `BiomeZone::Alpha` and `BiomeZone::Beta`, evaluating population differentiation `F_ST = (H_total - H_sub) / H_total`.

---

## 7. Prognostic Planetary Biogeochemistry (`O2 PAL Index`, `CO2`, `N/P`)

In coupled `PlanetaryBiosphere` mode (`src/model.rs::step_biogeochemical_feedbacks`), `atmospheric_o2` models the normalized **`O2 (PAL Index)`** (Present Atmospheric Level index bounded to `[0.02, 0.94]`), where `O2 (PAL) >= 0.25` marks the Great Oxidation Event (`GOE`) threshold:

```text
13. Coupled Atmosphere & Nutrient Equations:
    d(NP)/dt  = gain * (0.115 * F_syntroph + 0.065 * F_superorg + 0.080 * s_gaia2 * F_steward - 0.055 * F_autotroph)
    d(O2)/dt  = gain * (0.085 * F_auto_eff + 0.035 * s_gaia2 * F_steward - 0.028 - 0.045 * O2 - 0.095 * b_burn * F_tech * (1.0 - 0.75 * s_gaia2))
    d(CO2)/dt = gain * (14.0 + 135.0 * b_burn * (0.35 + F_tech) * (1.0 - 0.82 * s_gaia2) - 38.0 * F_auto_eff - 26.0 * F_superorg - 95.0 * s_gaia2 * (0.25 + F_steward))
    E_planet  = clamp(J_bath * (0.72 + 0.48 * NP * (0.4 + 0.6 * O2) + 0.52 * s_gaia2 * min(1.4, F_steward)) * (1.0 - 0.72 * B_load), 0.22, 2.35)
```

- **`PlanetaryState`**: Transitions emergently among `AnoxicArchean` (`O2 < 0.25`), `OxygenatedPhanerozoic`, `AnthropoceneDisequilibrium` (`boundary_transgression_load >= 0.42`), and `HomeostaticAwakening` (`Gaia-2.0`, requiring `verified_stewardship_flux >= 0.25`, `s_gaia2 >= 0.65`, `O2 >= 0.55`, and `CO2 <= 380 ppm`).

---

## 8. The 8 Continuous Epochs, 13 Emergent Regimes & 24 Milestones

| Epoch (`GrandEpochStage`) | Age (`Ga`) | Active Channels | Core Physical & Evolutionary Transitions |
| --- | --- | --- | --- |
| **E1** `E1_Homochirality` | `4.250` | `1` | Magnetite `Fe3O4` `CISS` `D`-RNA symmetry breaking + `3'-NCCA` `L`-aminoacylation |
| **E2** `E2_SGC_RAF` | `4.180` | `1` | Citrate-chelated `Mg2+` `RAF` protocells + collective `HGT` `20-AA` `SGC` convergence |
| **E3** `E3_TDNA_LUCA` | `4.050` | `2` | Anaerobic Class III `RNR` `RNA -> U-DNA -> T-DNA`, `ThyX/UDG` repair, virosphere & `LUCA` |
| **E4** `E4_GOE_LECA` | `1.850` | `2` | Oxygenic photosynthesis (`GOE`), Asgard-alphaproteobacterial `LECA`, `1°/2°/3°` plastids & `MRO` reduction |
| **E5** `E5_GermSoma` | `0.800` | `3` | Clonal `r = 1.0` zygotic bottleneck, `OXPHOS` germ-soma sequestration & mitochondrial apoptosis |
| **E6** `E6_Superorg_Agri` | `0.060` | `3` | Monogamy-window eusocial superorganisms, worker egg policing & `8n` clonal crops + *Pseudonocardia* |
| **E7** `E7_Archives_FormalAI` | `0.005` | `4` | Symbolic cumulative culture, peer-verified written archives & formal-verifier anti-collapse AI |
| **E8** `E8_Gaia2_Homeostasis` | `0.000` | `5` | Closed-loop `O2 (PAL)`, `CO2`, `N/P`, Gibbs disequilibrium & renewable solar/geothermal exergy homeostasis |

The emergent regime classifier (`src/regimes.rs::classify_emergent_state`) assigns a 5-component `DiagnosticVector(V_oc, M_c, R_c, V_c, P_c)` solely from live observables (`Phi_net`, `v_O`, plasticity gain, birth activity, `rho_Adami`, `f_blend` vs `theta_f`, Price covariance, crowding, assimilation, and recombination) across 13 regimes (`R(0.0)` through `R(4.0-AUT)`).

---

## 9. Software License

`itme_evolution_sim` is free software licensed under the **GNU General Public License v3.0 or later (`GPL-3.0-or-later`)** — see [`LICENSE`](../LICENSE).
