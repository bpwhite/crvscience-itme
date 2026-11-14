# ITME v23.0 (`v2.3.0`) — Continuous 4.25 Ga Multi-Epoch Hadean-to-Gaia Grand Trajectory Engine

Native Rust simulation unifying **24 evolutionary milestones** into a single unbroken **8-Epoch Continuous Trajectory (`4.25 Ga -> Gaia 2.0`)** on one persistent `SimulationWorld` and one reticulate `PhylogenyTracker`:
1. **Epoch 1 (`~4.25 Ga`) — Magnetite `CISS` Homochirality & Minihelix `L-Aminoacylation`** ([Ozturk et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10246896/); [Tamura & Schimmel, 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1564265/))
2. **Epoch 2 (`~4.18 Ga`) — Collective `20-AA SGC` Crystallization & `RAF` Protocells** ([Vetsigian et al., 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1502294/); [Damer & Deamer, 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7133448/); [Hordijk et al., 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6315399/))
3. **Epoch 3 (`~4.05 Ga`) — `RNR` `RNA -> T-DNA`, `ThyX/UDG` Repair, Virosphere & `LUCA`** ([Lundin et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4390871/); [Sousa et al., 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2660913/); [Forterre, 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1450140/))
4. **Epoch 4 (`~2.40–1.85 Ga`) — `GOE` Oxygenation & Asgard-Alphaproteobacterial Symbiogenesis (`LECA`)** ([Krissansen-Totton et al., 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC5787383/); [Eme / Vosseberg et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12872458/))
5. **Epoch 5 (`~0.80 Ga`) — Clonal Germ-Soma Multicellularity & Mitochondrial Apoptosis** ([Ratcliff et al., 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3277146/); [Szathmáry, 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4547294/))
6. **Epoch 6 (`~0.06 Ga`) — Eusocial Superorganisms & Clonal Agricultural Ectosymbioses** ([Boomsma, 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2781870/); [Schultz & Brady, 2008](https://pmc.ncbi.nlm.nih.gov/articles/PMC2291119/))
7. **Epoch 7 (`~0.01 Ga–Present`) — Symbolic Culture, Peer-Verified Archives & Machine-Verified Synthesis** ([Lewis & Laland, 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3385684/); [Smaldino & McElreath, 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5043322/); [Shumailov et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11269175/))
8. **Epoch 8 (`Gaia-2.0`) — Planetary Biogeochemical Closure & Homeostatic Awakening** ([Kleidon, 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3261436/); [Kempes et al., 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5686401/); [Wong & Bartlett, 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9065981/))

## Quick Start

```bash
cargo run --release                                                                                       # Interactive Live Multi-Epoch GUI
cargo run --release -- --headless --grand-trajectory --export-newick grand.nhx --export-phylo-svg grand.svg # Continuous 8-Epoch Grand Benchmark
cargo test                                                                                                # Run Full 24-Suite Verification
```

## Core Evolution Equation System (ITME v2.3.0)

All mathematical expressions in the manuscript are typeset in plain-text ASCII notation so undergraduate readers can map every variable directly to executable lines in `src/model.rs`, `src/organisms.rs`, and `src/engine.rs` without LaTeX macros:

```text
1. Biospheric / Ecosystem Order (Macro-Evolutionary State Equation):
   O_pop(t) = O(t) * ln(1 + N_pop(t)) * M_div(t) * M_web(t) * M_ETI(t)

2. Lineage Order and Order Velocity (Phase-Transition Criterion):
   O(t)     = S_T_minus(t) * L_I_bar(t)
   v_O(t)   = sum_k((k - k_bar) * (O_k - O_bar)) / sum_k((k - k_bar)^2)

3. Fitness-Weighted Population Informational Leverage:
   L_I_bar(t) = sum_i(w_i * L_I_i) / sum_i(w_i)
   L_I_i      = C_i / max(S_min, S_eff_i)

4. Asymmetric Downside Stability (Thermodynamic + Error-Cliff Survival):
   S_T_minus(t) = clamp((1.0 - H_comb(t)) * (1.0 - P_downside(t)), 0.02, 1.0)
   H_comb(t)    = clamp(H_thermo(t) + H_meltdown(t), 0.0, 0.94)
   P_downside   = (0.25 + N_losses) / (4.25 + N_pairs)

5. Price Selection-Transmission Decomposition of Evolutionary Change:
   Delta_z_bar = Cov(w, z) / w_bar + E[w * Delta_z] / w_bar

6. Emergent Functional Complexity (Sublinear Structural Scaffold x Epistasis x 8-Epoch Synergy):
   C_i = max(0.08, 4.25 * (S_eff_i ^ 0.62) * (0.38 + 1.45 * E_epistatic) * Phi_exergy * M_eng * Omega_synergy * I_int * (1.0 + 0.72 * m_bar * pi_mem))

7. Effective Conserved Structure (LZ76 Compression Blended with Adami Ensemble Conservation):
   S_raw_i   = S_host + 0.65 * S_endo + 0.55 * S_nm
   S_seq     = c_LZ76(seq) * (log2(A) + 0.35 * log2(L)) + E_wiring * (0.25 * log2(L))
   S_eff_i   = S_min + max(0.0, S_raw_i - S_min) * ((1.0 - lambda_N) * rho_phys + lambda_N * rho_Adami)
   lambda_N  = N_rep / (N_rep + 8.0)
   rho_Adami = 1.0 - (1.0 / L_eval) * sum_pos(H_pos / H_max)

8. Net Exergy Transduction and Landauer Maintenance Cost:
   Phi_net    = W_harvest - Q_Landauer
   W_harvest  = (W_primary + W_syntrophic) * M_origin + J_trophic
   Q_Landauer = k_Landauer * (S_host + S_endo + S_nm + S_klepto + S_coord_import + 2.8 * norm1(m))
   Phi_exergy = clamp(0.76 + 0.48 * min(1.6, Phi_net), 0.65, 1.72)   [if Phi_net >= -0.02]
              = clamp(0.28 * I_int, 0.008, 0.32)                     [if Phi_net <  -0.02]

9. Prebiotic-to-Cellular Origin Gate Product (Epochs 1 to 3) and Post-LUCA Synergy Product (Epochs 4 to 8):
   M_origin      = m_chiral_SGC * m_abio_RAF * m_DNA_virosphere * m_Gaia_atm
   Omega_synergy = Omega_holo_MRO * Omega_multicell * Omega_superorg * Omega_agro * Omega_culture * Omega_archive * Omega_formal_AI * Omega_Gaia

10. Multi-Channel Heritable Mutual Information and Prognostic Biosphere Feedbacks:
   I_herit_total = I_chan(f_gen, A, S_raw) + I_chan(f_cult, A, B_cult) + I_chan(f_arch, A, B_arch) + I_chan(f_proof, A, B_proof)
   I_chan(f,A,B) = B * (1.0 - (-(f * log2(f) + (1.0 - f) * log2((1.0 - f) / (A - 1.0)))) / log2(A))
   E_planet      = clamp(J_bath * (0.72 + 0.48 * NP * (0.4 + 0.6 * O2) + 0.52 * s_gaia2 * min(1.4, F_steward)) * (1.0 - 0.72 * B_load), 0.22, 2.35)

11. Phase-Transition Ratchet Sub-Equations (Homochirality, SGC, Eigen Cliff, Deamination, Kin Selection, Autophagy & Atmosphere):
   P_cross_inhib = 0.82 * (1.0 - ee_D)^2
   MP            = (phi_random_mean - phi_code) / (phi_random_mean - phi_theoretical_min)
   L_eigen_max   = floor(ln(sigma_superiority) / mu_per_base)
   U_lesion(t+1) = clamp(0.62 * U_lesion(t) + J_deamination - 0.56 * e_udg, 0.008, 0.95)
   Kin gate      : r_relatedness * B_recipient_benefit > C_actor_cost
   A_autophagy   = clamp(A_autophagy + r_synth * (1.0 - v_formal) - 1.35 * gamma_empirical, 0.0, 1.0)
   d(NP)/dt      = gain * (0.115 * F_syntroph + 0.065 * F_superorg + 0.080 * s_gaia2 * F_steward - 0.055 * F_autotroph)
   d(O2)/dt      = gain * (0.085 * F_auto_eff + 0.035 * s_gaia2 * F_steward - 0.028 - 0.045 * O2 - 0.095 * b_burn * F_tech * (1.0 - 0.75 * s_gaia2))
   d(CO2)/dt     = gain * (14.0 + 135.0 * b_burn * (0.35 + F_tech) * (1.0 - 0.82 * s_gaia2) - 38.0 * F_auto_eff - 26.0 * F_superorg - 95.0 * s_gaia2 * (0.25 + F_steward))
```

### Table 3A: Macro-Evolutionary Order, Stability, and Selection Variables (Equations 1 to 5)

| Symbol | Variable Name | Domain / Units | Definition and Theoretical Grounding |
| ----- | ----- | ----- | ----- |
| `O_pop(t)` | Biospheric Ecosystem Order | `[0.0, +inf)` dimensionless | Total stability-weighted informational order scaled by population size, clade richness, trophic entropy, and major transition tiers ([Szathmáry, 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4547294/)). |
| `O(t)` | Lineage Order | `[0.0, +inf)` dimensionless | Product of downside persistence stability `S_T_minus(t)` and mean informational leverage `L_I_bar(t)`. |
| `v_O(t)` | Order Velocity | `(-inf, +inf)` per tick | Ordinary-least-squares linear slope of `O(t)` across telemetry window `W_dt` (`k` indexes window ticks, `k_bar` and `O_bar` are window means); positive `v_O` diagnoses active adaptive expansion. |
| `L_I_i`, `L_I_bar` | Informational Leverage | `[0.04, +inf)` ratio | Emergent functional complexity `C_i` sustained per bit of effective conserved structure `S_eff_i`, averaged across individuals weighted by fitness `w_i`. |
| `N_pop(t)` | Extant Population Size | Integer `[0, K_cap]` | Total viable replicating units on the `16x16` spatial reaction-diffusion lattice; enters `O_pop` logarithmically via `ln(1 + N_pop)`. |
| `M_div(t)` | Clade Diversity Multiplier | `[1.0, +inf)` factor | `1.0 + 0.16 * max(0, N_species - 1)`, rewarding reproductive and ecological speciation (`N_species` distinct clades). |
| `M_web(t)` | Trophic Web Multiplier | `[1.0, 1.36]` factor | `1.0 + 0.18 * H_web`, where `H_web` (`0.0` to `2.0` bits) is Shannon guild entropy across Producers, Syntrophs, Predators/Parasites, and Holobionts ([Boyle et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC11848429/)). |
| `M_ETI(t)` | Population Transition Bonus | `[1.0, +inf)` factor | `1.0 + sum_j(beta_j * n_j)`, summing per-individual individuality increments across holobionts, serial plastids, MROs, germ-soma organisms, superorganisms, agro-colonies, cultural demes, verified archives, and formal AI symbionts ([Szathmáry, 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4547294/)). |
| `S_T_minus(t)` | Asymmetric Downside Stability | `[0.02, 1.0]` probability | Survival complement against combined thermodynamic/mutational hazards `H_comb` and empirical leverage drawdowns `P_downside`. |
| `H_thermo(t)` | Thermodynamic Starvation Hazard | `[0.0, 0.92]` probability | Risk of dissipative collapse when net exergy margin `Phi_net` falls negative, buffered by local niche construction capital ([Kleidon, 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3261436/)). |
| `H_meltdown(t)` | Error-Catastrophe Hazard | `[0.0, 0.58]` probability | `clamp(2.4 * (theta_f - f_blend), 0.18, 0.58)` whenever blended replication fidelity `f_blend` drops below Eigen's critical threshold `theta_f = 1.0 - ln(sigma) / D_eff` ([Wolf & Koonin, 2007](https://pmc.ncbi.nlm.nih.gov/articles/PMC1894784/)). |
| `P_downside(t)` | Downside Drawdown Frequency | `[0.0, 0.96]` probability | Laplace-smoothed fraction of adjacent tick pairs (`N_losses` out of `N_pairs`) where `L_I(t) - L_I(t - 1)` falls below `-epsilon` (`epsilon = 0.08` default noise floor). |
| `w_i`, `w_bar` | Darwinian Fitness | `[0.02, 4.60]` relative rate | Composite reproductive viability `(0.34 * E_epistatic + 0.46 * (L_I_i / 2.4) + 0.20 * V_energy) * I_int` and its population mean `w_bar`. |
| `Delta_z_bar` | Price Trait Shift | `(-inf, +inf)` trait/tick | Decomposes per-generation trait change `z` into natural selection covariance `Cov(w, z) / w_bar` plus transmission/plasticity bias `E[w * Delta_z] / w_bar`. |

### Table 3B: Structural Compression, Adami Conservation, and Landauer Energetics (Equations 6 to 8)

| Symbol | Variable Name | Domain / Units | Definition and Theoretical Grounding |
| ----- | ----- | ----- | ----- |
| `C_i` | Emergent Functional Complexity | `[0.08, +inf)` bits-equiv | Realized adaptive work-and-information capacity combining sublinear structural scaling (`S_eff ^ 0.62`), epistasis, exergy margin, transition synergies, and lifetime plasticity. |
| `S_raw_i` | Raw Multi-Genome Structure | `[2.0, +inf)` bits | Weighted description length across host chromosome (`S_host`), endosymbiont/apicoplast (`0.65 * S_endo`), and periplastidial nucleomorph (`0.55 * S_nm`). |
| `c_LZ76(seq)` | Lempel-Ziv Phrase Count | Integer `[1, L]` phrases | Number of distinct dictionary phrases parsed by Lempel-Ziv 1976 compression of integer sequence `seq`, penalizing tandem repetitive bloat. |
| `A`, `L` | Alphabet Size and Sequence Length | Integers (`A` at least `2`, `L` at least `1`) | Symbol cardinality `A` (`12` molecular, `10` hominid/AI, `20` `SGC` amino acids) and polymer/locus count `L`. |
| `E_wiring` | Epistatic Wiring Edge Count | Integer `[0, L]` edges | Number of active intramolecular stem-loops, regulatory links, or algorithmic pipeline edges in the unit's interaction graph. |
| `S_eff_i`, `S_min` | Effective Conserved Structure | `[2.0, S_raw_i]` bits | Non-neutral structural bits conserved against neutral drift above minimal floor `S_min = 2.0` bits ([Kempes et al., 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5686401/)). |
| `lambda_N` | Ensemble Blending Weight | `[0.0, 1.0)` weight | `N_rep / (N_rep + 8.0)` (`0` when `N_rep` is `1` or `0`), smoothly shifting structural valuation from single-molecule physical folding `rho_phys` to population Adami conservation `rho_Adami`. |
| `rho_phys` | Physical Constraint Fraction | `[0.05, 1.0]` fraction | Intrinsic thermodynamic/structural folding constraint active even in pre-Darwinian (`N_rep` at most `1`) dissipative assemblies. |
| `rho_Adami` | Adami Ensemble Conservation | `[0.08, 1.0]` fraction | `1.0 - mean_pos(H_pos / H_max)`, measuring per-locus Shannon entropy reduction maintained by purifying selection across `N_rep` sequences ([Kempes et al., 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5686401/)). |
| `E_epistatic` | Causal Epistatic Score | `[0.08, 2.10]` score | Context-specific functional synergy (catalytic triad + RNP cooperativity in molecules; ordered social-causal chains in hominids; framed memory-verification pipelines in AI). |
| `W_harvest` | Harvested Useful Work | `[0.0, +inf)` exergy/tick | Sum of primary solar/geochemical work `W_primary`, syntrophic byproduct cross-feeding `W_syntrophic`, and predatory/parasitic flux `J_trophic` ([Kleidon, 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3261436/)). |
| `Q_Landauer` | Landauer Maintenance Dissipation | `[0.0, +inf)` exergy/tick | Thermodynamic cost `k_Landauer` per bit required to replicate, proofread, import proteins, and maintain somatic/neural/archival states ([Kempes et al., 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5686401/)). |
| `S_coord_import` | Organelle & Social Overhead Bits | `[0.0, +inf)` bits | Sum of `EGT`/nucleomorph/apicoplast protein-import costs, `MRO` `ISC` maintenance, multicellular soma coordination, eusocial caste signaling, agro-crypt upkeep, and offloaded neural load `cap_neural * (1.0 - 0.58 * f_offload)` ([Pegado, 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9226919/)). |
| `m_bar`, `norm1(m)`, `pi_mem` | Plastic Memory Terms | `[0.0, 1.0]` state/gain | Mean lifetime memory state `m_bar`, L1 memory magnitude `norm1(m)` contributing to Landauer upkeep, and environmental learning plasticity gain `pi_mem`. |
| `I_int` | Structural Integrity | `[0.02, 1.0]` health | Multiplicative viability health eroded by unchelated `Mg2+` vesicle precipitation, Eigen error cliffs, `C->U` mutator meltdowns, neoplastic cheaters, model autophagy, or planetary burnout. |

### Table 3C: 8-Epoch Origin Gates, Transition Synergies, and Biosphere Couplings (Equations 9 to 11)

| Symbol | Variable Name | Domain / Units | Definition and Theoretical Grounding |
| ----- | ----- | ----- | ----- |
| `m_chiral_SGC` | Homochirality & `SGC` Gate | `{0.22 .. 1.00}` factor | Scales work from racemic `3'` cross-inhibited baths (`0.22`) through magnetite `CISS` `D`-RNA (`0.54`, [Ozturk et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10246896/)), `3'-NCCA` `L`-aminoacylation (`0.76`, [Tamura & Schimmel, 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1564265/)), and 20-AA `SGC` (`1.00`, [Vetsigian et al., 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1502294/)). |
| `m_abio_RAF` | Hadean Abiogenesis Gate | `{0.28 .. 1.00}` factor | Progresses across monomer soup (`0.28`), citrate-`Mg2+` `RAF` protocells (`0.52`, [Adamala & Szostak, 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC4104020/); [Hordijk et al., 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6315399/)), ribozyme hypercycles (`0.74`, [Papastavrou et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10990123/)), progenotes (`0.88`), and `LUCA` (`1.00`, [Koonin, 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4008548/)). |
| `m_DNA_virosphere` | `RNR` `T-DNA` & Virosphere Gate | `{0.46 .. 1.00}` factor | Scales chromosome stability from labile `2'-OH` RNA (`0.46`) through radical-SAM `RNR` `U-DNA` (`0.68`, [Lundin et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4390871/)), `ThyX/UDG`-repaired `T-DNA` (`0.88`, [Sousa et al., 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2660913/)), and viral polymerase symbiogenesis (`1.00`, [Forterre, 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1450140/)). |
| `Omega_holo_MRO` | Symbiogenetic & `MRO` Synergy | `[0.32, 1.85]` factor | Couples Asgard-alphaproteobacterial `LECA` fusion ([Dominant Asgard Study, 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12872458/)), nuclear `EGT` lock-in, `2°/3°` plastids, apicoplast streamlining ([McFadden & Yeh, 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5406208/)), and anaerobic hydrogenosome/mitosome `ISC` retention ([Stairs et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4571565/)). |
| `Omega_multicell` | Clonal Germ-Soma Synergy | `[0.18, 1.34]` factor | `(1.0 + 0.34 * f_soma) * (1.0 - 0.85 * L_cheat)`, rewarding `OXPHOS` soma-germline specialization when unicellular bottlenecks and mitochondrial apoptosis suppress somatic cheater load `L_cheat` ([Ratcliff et al., 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3277146/); [Radzvilavicius et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5172535/)). |
| `Omega_superorg` | Eusocial Superorganism Synergy | `[0.16, 1.42]` factor | `(1.0 + 0.42 * f_caste) * (1.0 - 0.86 * L_rebel)`, rewarding sterile worker/soldier caste division `f_caste` stabilized by monogamy and worker egg policing ([Boomsma, 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2781870/)). |
| `Omega_agro` | Clonal Agriculture Synergy | `[0.12, 1.38]` factor | `(1.0 + 0.38 * Y_gongyl) * (1.0 - 0.88 * L_pathogen)`, coupling `8n` clonal cultivar yield `Y_gongyl` to *Pseudonocardia* antimicrobial suppression of garden blight `L_pathogen` ([Schultz & Brady, 2008](https://pmc.ncbi.nlm.nih.gov/articles/PMC2291119/)). |
| `Omega_culture` | Cumulative Culture Synergy | `[0.12, 1.76]` factor | `(1.0 + 0.36 * (T_ratchet / 4.0)) * c_factor * (1.0 - 0.86 * L_drift)`, peaking under high-fidelity symbolic teaching (`f_cult` at least `0.82`) and intermediate inter-deme connectivity (`0.12` to `0.75`) ([Lewis & Laland, 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3385684/); [Derex & Boyd, 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4801235/)). |
| `Omega_archive` | Exosomatic Institutional Synergy | `[0.12, 1.48]` factor | `1.48 * (1.0 - 0.88 * L_FDR)` in peer-verified written archives, penalizing unvetted false-discovery flood `L_FDR` ([Smaldino & McElreath, 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5043322/)). |
| `Omega_formal_AI` | Neuro-Symbolic Proof Synergy | `[0.08, 1.54]` factor | `1.54 * D_tail * (1.0 - 0.90 * L_autophagy)`, sustaining theorem discovery when formal verifiers (`v_formal` at least `0.72`) and empirical data anchors (`gamma_emp` at least `0.12`) preserve tail diversity `D_tail` ([Shumailov et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11269175/); [Romera-Paredes et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10794145/)). |
| `Omega_Gaia`, `E_planet` | Planetary Homeostasis Coupling | `{0.74 .. 1.36}`, `[0.22, 2.35]` | Couples guild `O2`, `CO2`, and `NP` nutrient recycling (`F_steward`, `s_gaia2`, `b_burn`) to planetary exergy budget `E_planet` and boundary transgression load `B_load` ([Richardson et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10499318/); [Wong & Bartlett, 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9065981/)). |
| `I_herit_total` | Multi-Channel Mutual Information | `[0.0, +inf)` bits | Sum of non-erased Shannon mutual information transmitted across genetic/organellar (`f_gen`, `S_raw`), cultural (`f_cult`, `B_cult`), archival (`f_arch`, `B_arch`), and formal theorem (`f_proof`, `B_proof`) channels. |

## Citation & CRV Science Apps Hub

- **CRV Science Apps Hub**: <https://crvscience.com/apps>
- **Companion Preprint (`APP-PREPRINT-07`)**: [ITME: Continuous Four-Billion-Year Information-Theoretic Macroevolution and Planetary Biogeochemical Feedback Simulation](https://crvscience.com/apps/preprints/itme-planetary-biosphere-evolution)
- **Concept DOI**: [`10.5281/zenodo.24222386`](https://doi.org/10.5281/zenodo.24222386)
- **Source Repository**: <https://github.com/bpwhite/crvscience-itme>

See [`CITATION.cff`](CITATION.cff) for machine-readable BibTeX/CFF metadata and [`validation/README.md`](validation/README.md) for the headless validation benchmark outputs (`cargo test --test headless_validation`).

