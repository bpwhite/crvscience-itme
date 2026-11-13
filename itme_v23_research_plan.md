# Research Plan: Unifying Evolutionary Phase Transitions from Prebiotic Homochirality to Planetary Homeostasis in the Information-Theoretic Model of Evolution (ITME v2.3.0)

**Document Type:** Stage-1 Preliminary Research Plan and Quantitative Benchmark Blueprint

**Target Venue Style:** Undergraduate Peer-Reviewed Interdisciplinary Biology and Systems Science Journal

**Target Full-Manuscript Length:** 4,600 words (range: 4,400 to 4,800 words, excluding tables and reference list)

**Simulation Codebase:** `itme_evolution_sim` v2.3.0 (12 Rust modules, 24 integration verification suites, deterministic Xoshiro256 PRNG)

**Formatting Standard:** Plain-text ASCII equations, standard tabular summaries, non-hyperbolic empirical prose

## 1. Executive Editorial Scope and Research Thesis

Living systems persist far from thermodynamic equilibrium by coupling free-energy transduction to heritable information storage ([Kleidon, 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3261436/); [Kempes et al., 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5686401/)). Across Earth's 4.25-billion-year history, increases in organismal and ecological complexity did not accumulate through smooth, unconstrained drift alone; instead, lineage persistence repeatedly reorganized across discrete Major Evolutionary Transitions in Individuality (`ETI`), each introducing new mechanisms of error suppression, conflict policing, division of labor, and multi-channel inheritance ([Szathmáry, 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4547294/)).

This preliminary research plan outlines a 4,600-word undergraduate-journal article presenting **ITME v2.3.0**, a context-free computational framework implemented in Rust that chains eight sequential transitions inside a single continuous 240-tick simulation (`--grand-trajectory`):

1. **Epoch 1 (4.250 Ga):** Prebiotic `D/L` homochiral symmetry breaking on magnetized magnetite (`Fe3O4`) surfaces ([Ozturk et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10246896/)) and stereoselective `3'-NCCA` RNA minihelix `L`-aminoacylation ([Tamura & Schimmel, 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1564265/)).

2. **Epoch 2 (4.180 Ga):** Collective horizontal gene transfer (`HGT`) crystallization of the 20-amino-acid Standard Genetic Code (`SGC`) ([Vetsigian, Woese & Goldenfeld, 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1502294/)) inside citrate-chelated `Mg2+` Reflexively Autocatalytic Food-Generated (`RAF`) protocells ([Adamala & Szostak, 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC4104020/); [Hordijk et al., 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6315399/)).

3. **Epoch 3 (4.050 Ga):** Anaerobic Class III radical-SAM ribonucleotide reductase (`RNR`) `2'-OH` to `2'-H` deoxyribose synthesis ([Lundin et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4390871/)), thymidylate synthase plus uracil-DNA glycosylase (`UDG`) cytosine-deamination repair ([Sousa et al., 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2660913/)), primordial virosphere polymerase domestication ([Forterre, 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1450140/)), and crossing of Carl Woese's Darwinian threshold into the Last Universal Common Ancestor (`LUCA`) ([Koonin, 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4008548/)).

4. **Epoch 4 (1.850 Ga):** Cyanobacterial oxygenic photosynthesis driving the Great Oxidation Event (`GOE`) ([Krissansen-Totton et al., 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC5787383/)) and endosymbiotic integration of an alphaproteobacterial respiratory symbiont within an Asgard archaeal host to form the Last Eukaryotic Common Ancestor (`LECA`) ([Dominant Asgard Contribution Study, 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12872458/)).

5. **Epoch 5 (0.800 Ga):** Clonal multicellularity via incomplete cytokinesis, single-cell zygotic bottlenecks (`r = 1.0`), oxidative-phosphorylation (`OXPHOS`) soma-germline division of labor, and mitochondrial apoptosis policing somatic defectors ([Ratcliff et al., 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3277146/); [Radzvilavicius et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5172535/); [Aktipis et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4581024/)).

6. **Epoch 6 (0.060 Ga):** Monogamy-window eusocial superorganisms with sterile worker castes, mutual egg policing ([Boomsma, 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2781870/); [West et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4547252/)), polyploid clonal crop domestication, and cuticular antimicrobial actinobacterial symbionts ([Schultz & Brady, 2008](https://pmc.ncbi.nlm.nih.gov/articles/PMC2291119/); [Shik et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11451569/)).

7. **Epoch 7 (0.005 Ga):** Partially connected symbolic cumulative culture ([Lewis & Laland, 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3385684/); [Derex & Boyd, 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4801235/)), written exosomatic archives decoupling cultural capacity from individual brain volume ([Pegado, 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9226919/)), peer-verified epistemic institutions ([Smaldino & McElreath, 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5043322/)), and neuro-symbolic proof verification preventing recursive synthetic autophagy ([Shumailov et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11269175/); [Romera-Paredes et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10794145/)).

8. **Epoch 8 (0.000 Ga / Gaia-2.0):** Planetary-scale biogeochemical co-regulation across autotrophic, syntrophic, superorganismal, and machine-verified guilds stabilizing atmospheric `O2`, `CO2`, `N/P` nutrient closure, and solar exergy capture ([Boyle et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC11848429/); [Richardson et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10499318/); [Wong & Bartlett, 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9065981/)).

## 2. Manuscript Section Architecture and Word-Count Budget

| Section Number and Title | Target Words | Primary Pedagogical Objective and Content Focus | Visual and Tabular Anchors | 
| ----- | ----- | ----- | ----- | 
| **Abstract** | 250 | Summarize the unified thermodynamic-informational framework, the 8-epoch Rust simulation design, monotonic expansion of Informational Leverage (`L_I`: `3.016` to `124.792`), and counterfactual knockout findings ([Szathmáry, 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4547294/); [Wong & Bartlett, 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9065981/)). | Standalone structured summary | 
| **1. Introduction: Free Energy, Heredity, and Phase Transitions** | 550 | Introduce open thermodynamic systems ([Kleidon, 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3261436/)), Landauer-bounded biological computation ([Kempes et al., 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5686401/)), and why major transitions require coupled error-control and conflict-policing ratchets ([Szathmáry, 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4547294/)). | Conceptual overview diagram | 
| **2. Mathematical Formalism and Context-Free Metrics** | 750 | Define Informational Leverage (`L_I`), Order (`O`), Order Velocity (`v_O`), Adami sequence conservation (`rho`), Price selection-transmission decomposition, and Extended Newick (`eNewick`) reticulate network encoding ([Cardona et al., 2008](https://pmc.ncbi.nlm.nih.gov/articles/PMC2621367/)). | Plain-Text Equations 1 to 10 | 
| **3. Prebiotic Geochemistry to Megabase LUCA (Epochs 1 to 3)** | 850 | Trace `Fe3O4` spin-selective `RAO` crystallization ([Ozturk et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10246896/)), `3'-NCCA` minihelix acylation ([Tamura & Schimmel, 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1564265/)), hot-spring wet-dry plus alkaline vent `RAF` protocells ([Damer & Deamer, 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7133448/); [Herschy et al., 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4247476/)), `SGC` convergence ([Vetsigian et al., 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1502294/)), and `RNR + ThyX/UDG` `T-DNA` synthesis ([Lundin et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4390871/); [Sousa et al., 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2660913/)). | Table 2 (Epochs 1 to 3) | 
| **4. Eukaryogenesis, Multicellularity, and Agro-Superorganisms (Epochs 4 to 6)** | 800 | Analyze cyanobacterial `GOE` oxygenation ([Krissansen-Totton et al., 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC5787383/)), Asgard-alphaproteobacterial symbiogenesis ([Dominant Asgard Contribution Study, 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12872458/)), reductive MRO/apicoplast pathways ([Stairs et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4571565/); [McFadden & Yeh, 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5406208/)), clonal germ-soma apoptosis ([Ratcliff et al., 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3277146/)), and attine-style polyploid agriculture ([Schultz & Brady, 2008](https://pmc.ncbi.nlm.nih.gov/articles/PMC2291119/)). | Table 2 (Epochs 4 to 6); Reticulate SVG Phylogram | 
| **5. Cultural Ratchets, Exosomatic Archives, and Gaia-2.0 (Epochs 7 to 8)** | 700 | Examine symbolic teaching thresholds ([Lewis & Laland, 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3385684/)), cognitive offloading into peer-verified archives ([Pegado, 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9226919/); [Smaldino & McElreath, 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5043322/)), formal proof filters against recursive model collapse ([Shumailov et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11269175/)), and closed-loop planetary homeostasis ([Wong & Bartlett, 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9065981/)). | Table 2 (Epochs 7 to 8) | 
| **6. Counterfactual Knockout Synthesis and Methodological Limits** | 450 | Contrast the unbroken 8-epoch chain against single-gate knockouts (`c_Mg`, `b_ciss`, `h_sgc`, `e_udg`, symbiogenesis, apoptosis, formal verification, and `s_gaia2`), and state idealizations of discrete-tick lattice dynamics ([Boyle et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC11848429/)). | Table 3 (Knockout Matrix) | 
| **7. Undergraduate Keyword Glossary** | 250 | Provide concise, jargon-free operational definitions for 16 foundational thermodynamic, geochemical, phylogenetic, and information-theoretic terms used across the article. | Section 6 Glossary List | 
| **Total Planned Manuscript Text** | **4,600** | Complete Stage-2 full-draft target (within 2,500 to 5,000 word window) | 3 Data Tables + 10 Plain-Text Equations | 

## 3. Core Plain-Text Evolution Equation System and Complete Variable Definitions

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

## 4. Quantitative Synthesis Tables from ITME v2.3.0 Benchmarks

### Table 1: Continuous 240-Tick Multi-Epoch Grand Trajectory (`4.250 Ga` to `Gaia-2.0`)

Values below are taken directly from the deterministic 240-tick `--grand-trajectory` benchmark (`30` ticks per epoch, `seed = 2301`, exported in `itme_v23_grand_trajectory_benchmark.csv`).

| Epoch & Code | Age (Ga) | End Tick | Pop | Active Channels | Heritable MI (bits) | Leverage `L_I` | Atm `O2` | Atm `CO2` (ppm) | Planetary Exergy | Key Transition Mechanism and Primary Reference | 
| ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | 
| **E1** `E1_Homochirality` | 4.250 | 30 | 24 | 1 | 1.58 | 3.016 | 0.06 | 1240.1 | 0.48 | Magnetite `CISS` (`99% D-ee`) + `3'-NCCA` `4:1` `L`-aminoacylation ([Ozturk et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10246896/); [Tamura & Schimmel, 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1564265/)) | 
| **E2** `E2_SGC_RAF` | 4.180 | 60 | 24 | 1 | 4.32 | 11.547 | 0.07 | 1180.2 | 0.62 | Citrate-`Mg2+` `RAF` vesicles + collective `HGT` 20-AA `SGC` (`MP = 81.1%`) ([Adamala & Szostak, 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC4104020/); [Vetsigian et al., 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1502294/)) | 
| **E3** `E3_TDNA_LUCA` | 4.050 | 90 | 29 | 2 | 18.67 | 31.261 | 0.08 | 1120.2 | 0.74 | Class III `RNR` + `ThyX/UDG` megabase `T-DNA` (`1,873 kb`) + `LUCA` ([Lundin et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4390871/); [Sousa et al., 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2660913/); [Forterre, 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1450140/)) | 
| **E4** `E4_GOE_LECA` | 1.850 | 120 | 34 | 2 | 34.85 | 44.583 | 0.78 | 342.0 | 1.18 | Cyanobacterial `GOE` + Asgard-alphaproteobacterial `LECA` fusion ([Krissansen-Totton et al., 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC5787383/); [Dominant Asgard Study, 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12872458/)) | 
| **E5** `E5_GermSoma` | 0.800 | 150 | 39 | 3 | 52.41 | 56.913 | 0.80 | 320.1 | 1.24 | Clonal `r = 1.0` zygote + `OXPHOS` germ-soma split + apoptosis ([Ratcliff et al., 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3277146/); [Radzvilavicius et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5172535/)) | 
| **E6** `E6_Superorg_Agri` | 0.060 | 180 | 44 | 3 | 78.87 | 74.251 | 0.82 | 290.1 | 1.32 | Kin-policed eusocial castes + `8n` clonal crops + *Pseudonocardia* ([Boomsma, 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2781870/); [Schultz & Brady, 2008](https://pmc.ncbi.nlm.nih.gov/articles/PMC2291119/)) | 
| **E7** `E7_Archives_FormalAI` | 0.005 | 210 | 49 | 4 | 332.08 | 96.413 | 0.80 | 354.8 | 1.48 | Symbolic culture + peer-verified archives + 64-locus proof kernel ([Lewis & Laland, 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3385684/); [Smaldino & McElreath, 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5043322/); [Shumailov et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11269175/)) | 
| **E8** `E8_Gaia2_Homeostasis` | 0.000 | 240 | 52 | 5 | 412.28 | 124.792 | 0.88 | 242.3 | 1.90 | Closed-loop `O2/CO2/NP` co-regulation + renewable solar exergy ([Kleidon, 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3261436/); [Boyle et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC11848429/); [Wong & Bartlett, 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9065981/)) | 

### Table 2: Counterfactual Knockout and Phase-Boundary Failure Matrix Across ITME Suites

Each row summarizes a deterministic ablation verified across the 24 integration tests and multi-regime Python benchmarks (`v1.1.0` through `v2.3.0`).

| Ablated Gate or Parameter | Critical Threshold Violated | Benchmark Outcome (`140` to `240` Ticks) | Mechanistic Failure or Stagnation Mode | Grounding Literature | 
| ----- | ----- | ----- | ----- | ----- | 
| **Wet-Dry & Vent Concentration** | `A_wd = 0.03`, `delta_pH = 0.15` | `Pop = 24`, stalled at `L = 10 nt` (`L_I = 0.626`) | Water-paradox hydrolysis outpaces phosphodiester condensation | [Damer & Deamer, 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7133448/); [Herschy et al., 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4247476/) | 
| **Citrate `Mg2+` Chelation** | `c_Mg = 0.10` at `50 mM Mg2+` | `Pop = 0` (`25` vesicle collapses) | Free `Mg2+` precipitates fatty-acid bilayers and cleaves RNA (`0.074 h^-1`) | [Adamala & Szostak, 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC4104020/) | 
| **Magnetite `CISS` Symmetry Bias** | `b_ciss = 0.04` (`ee_D = 0.03`) | `Pop = 0` (`25` chain terminations) | Opposite-handed `L`-ribonucleotides cap `3'` primers via cross-inhibition | [Ozturk et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10246896/); [Chen & Ma, 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC6974302/) | 
| **Ribozyme Replication Fidelity** | `mu = 0.052` (`L > L_eigen_max`) | `Pop = 0` (`24` error catastrophes) | Quasispecies master sequence melts across Eigen's error cliff | [Wolf & Koonin, 2007](https://pmc.ncbi.nlm.nih.gov/articles/PMC1894784/); [Papastavrou et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10990123/) | 
| **Communal Progenote `HGT`** | `h_sgc = 0.04` (`MP = 35.2%`) | `Pop = 0` (`24` Babel collapses) | Isolated lineages freeze into incompatible suboptimal genetic codes | [Vetsigian et al., 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1502294/); [Koonin & Novozhilov, 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC3293468/) | 
| **Thymine + `UDG` Excision** | `e_udg = 0.06` (`C->U = 0.774`) | `Pop = 0` (`24` mutator meltdowns) | Unrecognized `C -> U` cytosine deamination erodes `U-DNA` coding loci | [Sousa et al., 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2660913/); [Vertessy & Toth, 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2732909/) | 
| **Asgard-Mito Symbiogenesis** | Endosymbiosis disabled (`E4`) | `Pop = 24`, capped at `L_I = 31.907`, `2` channels | Viable prokaryotes persist without bioenergetic scaling to `LECA` or tissues | [Lane, 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC3996473/); [Dominant Asgard Study, 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12872458/) | 
| **Zygotic Bottleneck & Apoptosis** | Aggregative chimera (`r < 0.5`), no PCD | `Pop = 0` (somatic neoplastic overrun) | Unpoliced somatic defectors out-replicate cooperative tissue lineages | [Ratcliff et al., 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3277146/); [Aktipis et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4581024/) | 
| **Cuticular Antimicrobial Shield** | `8n` clonal crop without *Pseudonocardia* | `Pop = 0` (`26` *Escovopsis*-type blights) | Genetically uniform polyploid monoculture succumbs to specialized pathogens | [Schultz & Brady, 2008](https://pmc.ncbi.nlm.nih.gov/articles/PMC2291119/); [Shik et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11451569/) | 
| **Deme Network Connectivity** | `c_deme = 0.00` without written archives | `Pop = 0` (`24` skill-loss collapses) | Demographic drift erodes cumulative cultural repertoires (Tasmanian effect) | [Lewis & Laland, 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3385684/); [Derex & Boyd, 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4801235/) | 
| **Formal Proof & Empirical Anchor** | `v_formal < 0.72` or `gamma_emp < 0.12` | `Pop = 0` (`6.0%` tail diversity, autophagy) | Recursive training on unverified synthetic outputs erases distribution tails | [Shumailov et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11269175/); [Trinh et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10794143/) | 
| **Syntrophic `N/P` Remineralization** | `F_syntroph = 0.00` | `NP = 0.08`, stalled at `L_I = 3.005` | Refractory detritus locks up nutrients and throttles planetary `O2` and exergy | [Boyle et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC11848429/); [Component Gaia Study, 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC11892576/) | 
| **Planetary Stewardship Coupling** | `s_gaia2 = 0.06`, `b_burn = 0.92` | `Pop = 0` (`CO2 = 1800 ppm`, `24` burnouts) | Superlinear carbon extraction breaches planetary boundaries | [Richardson et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10499318/); [Wong & Bartlett, 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9065981/) | 

## 5. Detailed Section-by-Section Drafting Instructions for Stage 2

### Section 1: Introduction (Target: 550 words)

* Open with the thermodynamic baseline: Earth intercepts low-entropy solar radiation (`~5800 K`) and re-radiates infrared heat (`~255 K`), generating free energy (`exergy`) that abiotic and biotic dissipative structures convert into chemical disequilibrium ([Kleidon, 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3261436/)).

* Frame the central research question without dramatic metaphors: how do evolving lineages overcome physical degradation limits (hydrolysis, racemic chain termination, mutational error cliffs, somatic cheating, epistemic drift, and biogeochemical nutrient lock-up) across 4.25 billion years?

* State the contribution of ITME v2.3.0: a unified Markov-like discrete-step state and population engine where each Major Evolutionary Transition ([Szathmáry, 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4547294/)) acts as an error-correcting or conflict-policing ratchet that raises the ceiling on `I_herit` and `L_I`.

### Section 2: Mathematical Formalism and Context-Free Metrics (Target: 750 words)

* Derive Equations 1 and 2 in plain prose so an upper-division biology or computer science undergraduate can compute `C_emergent`, `S_eff`, `L_I`, `O`, and `v_O` on a hand calculator for a toy 16-locus sequence.

* Explain why `L_I` normalizes useful work plus heritable mutual information against Landauer erasure and maintenance cost `Q_Landauer`, connecting ribosomal translation efficiency near the Landauer bound ([Kempes et al., 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5686401/)) to organismal and institutional architectures.

* Introduce the Price equation (`Delta_z_bar = Cov(w, z)/w_bar + E(w * Delta_z)/w_bar`) as the observable separator between Darwinian selection covariance and non-Darwinian transmission bias, alongside Extended Newick (`eNewick`) `#H` hybrid tags for reticulate symbiogenesis, horizontal gene transfer, cultivar inheritance, and archival synthesis ([Cardona et al., 2008](https://pmc.ncbi.nlm.nih.gov/articles/PMC2621367/)).

### Section 3: Prebiotic Geochemistry to Megabase LUCA — Epochs 1 to 3 (Target: 850 words)

* Compare subaerial anoxic geothermal hot-spring pools (`K+/Na+ > 1.0`, hydration-dehydration-gel lipid encapsulation; [Damer & Deamer, 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7133448/); [Mulkidjanian et al., 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3325685/)) with submarine alkaline hydrothermal vents (`Fe(Ni)S` mineral pores, `delta_pH ~ 3.5`, Wood-Ljungdahl `Acetyl-CoA` carbon fixation; [Herschy et al., 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4247476/)). Explain why ITME v2.3.0's volcanic-island hydrological coupling combines wet-dry phosphodiester condensation with hydrothermal redox cofactors.

* Walk sequentially through the four molecular gates: (a) `Fe3O4` spin-selective `RAO` crystallization (`CISS`) relieving `3'` enantiomeric cross-inhibition ([Ozturk et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10246896/); [Joyce & Szostak, 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6120706/)); (b) Tamura-Schimmel `3'-NCCA` minihelix steric selection (`105 deg` Burgi-Dunitz angle, `4:1` `L:D` preference; [Tamura, 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3155382/)); (c) citrate `Mg2+` chelation protecting fatty-acid vesicles ([Adamala & Szostak, 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC4104020/)) plus `RAF` autocatalytic closure ([Xavier et al., 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7126077/)) and communal `HGT` convergence onto the `SGC` (`MP = 81.1%`; [Vetsigian et al., 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1502294/)); and (d) anaerobic Class III `RNR` (`2'-OH -> 2'-H`), `ThyX + UDG` cytosine-deamination repair, and viral polymerase domestication unlocking `1,873 kb` `T-DNA` genomes across Woese's Darwinian threshold ([Lundin et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4390871/); [Sousa et al., 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2660913/); [Forterre, 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1450140/)).

### Section 4: Eukaryogenesis, Multicellularity, and Agro-Superorganisms — Epochs 4 to 6 (Target: 800 words)

* Synthesize Epoch 4: cyanobacterial oxygenic photosynthesis drives atmospheric `O2` across the `GOE` (`O2 = 0.08 -> 0.78`), creating redox disequilibrium ([Krissansen-Totton et al., 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC5787383/)) while Asgard archaeal hosts integrate alphaproteobacterial respiratory symbionts into `LECA` ([Dominant Asgard Contribution Study, 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12872458/)). Briefly note ITME's reversible reductive modules in hypoxic or nutrient-dense niches (cryptic apicoplasts retaining `MEP/IPP` and `Fe-S` synthesis, [McFadden & Yeh, 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5406208/); anaerobic hydrogenosomes and `ISC`-retaining mitosomes, [Stairs et al., 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4571565/)).

* Synthesize Epochs 5 and 6: show how the exact same multilevel-selection logic (`r * B > C` plus active cheater policing) recurs first at the cellular level via unicellular zygotic bottlenecks (`r = 1.0`), `OXPHOS` soma-germline sequestration, and mitochondrial apoptosis ([Ratcliff et al., 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3277146/); [Radzvilavicius et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5172535/)), and subsequently at the colony level via lifetime queen monogamy (`r >= 0.50`), worker egg policing, polyploid `8n` fungal cultivar domestication, and *Pseudonocardia* antimicrobial defense ([Boomsma, 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2781870/); [Schultz & Brady, 2008](https://pmc.ncbi.nlm.nih.gov/articles/PMC2291119/)).

### Section 5: Cultural Ratchets, Exosomatic Archives, and Gaia-2.0 — Epochs 7 to 8 (Target: 700 words)

* Detail Epoch 7's three informational decoupling steps: (a) high-fidelity symbolic teaching (`f >= 0.82`) across partially connected deme networks (`c = 0.36`) preventing isolated skill loss while preserving cross-deme recombination ([Lewis & Laland, 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3385684/); [Derex & Boyd, 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC4801235/)); (b) written exosomatic archives and peer-verification ratchets (`r_peer >= 0.65`) reducing per-individual neural storage load from `90%` to `68%` ([Pegado, 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9226919/); [Smaldino & McElreath, 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5043322/)); and (c) neuro-symbolic formal proof kernels (`v_formal = 0.92`, `gamma_emp = 0.35`) preventing recursive synthetic variance collapse ([Shumailov et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11269175/); [Romera-Paredes et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC10794145/)).

* Detail Epoch 8: contrast unguided extractive energy scaling (`s_gaia2 = 0.06`, `CO2 = 1800 ppm`, planetary boundary breach) against machine-verified closed-loop stewardship (`s_gaia2 = 0.92`, `O2 = 0.88`, `CO2 = 242.3 ppm`, `Exergy = 1.90`, `L_I = 124.792`), connecting biospheric persistence selection ([Boyle et al., 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC11848429/)) to homeostatic awakening ([Wong & Bartlett, 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9065981/)).

### Section 6: Counterfactual Knockouts and Methodological Limitations (Target: 450 words)

* Walk through Table 2's 13 ablation rows to show that every transition gate in ITME v2.3.0 is mechanistically load-bearing rather than decorative: removing prebiotic error or repair filters causes population extinction (`Pop = 0`), whereas removing eukaryogenesis (`no_mito_symbio`) leaves a viable prokaryotic biosphere permanently capped at `2` inheritance channels and `L_I = 31.907`.

* Clearly state three modeling idealizations for undergraduate transparency: (1) coarse-grained discrete ticks compress unequal geological intervals onto uniform 30-tick relaxation windows; (2) spatial transport uses a 16x16 toroidal lattice rather than coupled 3D general-circulation ocean-atmosphere fluid dynamics; and (3) locus sequences (`16` to `64` integer slots) represent effective functional blocks rather than full nucleotide-level molecular dynamics.

## 6. Undergraduate Keyword Glossary (16 Core Terms)

* **Adami Complexity (`C_emergent`):** The amount of sequence information conserved across a population within a specific environment because purifying selection removes deleterious variants ([Kempes et al., 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5686401/)).

* **Chiral-Induced Spin Selectivity (`CISS`):** Interaction between electron spin polarization on magnetic minerals such as magnetite (`Fe3O4`) and chiral molecules that selectively crystallizes a single enantiomer (`D`-ribose precursors) from a racemic solution ([Ozturk et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10246896/)).

* **Cytosine Deamination (`C -> U`):** Spontaneous hydrolytic loss of an amine group converting cytosine into uracil; repairing it without corrupting `A:T` pairs requires methylating native uracil into thymine (`T`) and excising stray `U` via uracil-DNA glycosylase (`UDG`) ([Sousa et al., 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2660913/)).

* **Darwinian Threshold:** Carl Woese's transition point where rampant communal horizontal gene exchange among progenotes gives way to cohesive, vertically inherited cellular lineages (`LUCA`) ([Koonin, 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4008548/)).

* **Eigen Error Catastrophe:** Loss of sequence heredity when the per-base mutation rate `mu` exceeds `ln(sigma) / L`, causing a replicating polymer population to diffuse into random sequences ([Wolf & Koonin, 2007](https://pmc.ncbi.nlm.nih.gov/articles/PMC1894784/)).

* **Enantiomeric Cross-Inhibition:** Premature termination of non-enzymatic RNA template extension when a nucleotide of opposite handedness (`L`-ribose) incorporates at the `3'` primer end of a `D`-RNA strand ([Joyce & Szostak, 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6120706/)).

* **Endosymbiotic Gene Transfer (`EGT`):** Migration of genes from an endosymbiont or organelle genome into the host nuclear chromosome, locking a formerly autonomous partner into a unified holobiont ([McFadden & Yeh, 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5406208/)).

* **Exergy:** The maximum useful work obtainable as an open system or biosphere relaxes toward thermodynamic equilibrium with its surrounding heat bath ([Kleidon, 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3261436/)).

* **Exosomatic Memory:** Heritable symbolic information stored outside biological nervous systems in material notation, written archives, and verified digital repositories ([Pegado, 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9226919/)).

* **Extended Newick (`eNewick`):** A text standard (`#H` hybrid nodes plus `&&NHX` key-value annotations) representing phylogenetic networks with reticulate symbiogenesis and horizontal transfers rather than strictly branching trees ([Cardona et al., 2008](https://pmc.ncbi.nlm.nih.gov/articles/PMC2621367/)).

* **Homeostatic Awakening (`Gaia-2.0`):** A planetary transition in which cognitive and computational institutions couple energy expansion to closed-loop biogeochemical regulation inside safe planetary boundaries ([Wong & Bartlett, 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9065981/); [Richardson et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10499318/)).

* **Informational Leverage (`L_I`):** Context-free ITME ratio `(W_dot_useful + alpha_I * I_herit) / Q_Landauer` quantifying free-energy capture and heritable mutual information sustained per unit thermodynamic maintenance cost.

* **Major Evolutionary Transition in Individuality (`ETI`):** The integration of formerly independent replicating entities into a higher-level unit of selection via cooperative division of labor and conflict suppression ([Szathmáry, 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4547294/)).

* **Order Velocity (`v_O`):** Discrete time derivative `dO/dt` of stability-weighted informational leverage (`O = L_I * S_T_minus`), measuring the rate at which a lineage or biosphere expands organized thermodynamic capacity.

* **Reflexively Autocatalytic Food-Generated Set (`RAF`):** A chemical reaction network in which every reaction is catalyzed by at least one molecule produced within the network or supplied by an ambient geochemical food set ([Hordijk et al., 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6315399/)).

* **Recursive Model Collapse (Synthetic Autophagy):** Progressive loss of distribution tail variance when generative models repeatedly train on unverified synthetic outputs without empirical data anchoring or formal verification ([Shumailov et al., 2024](https://pmc.ncbi.nlm.nih.gov/articles/PMC11269175/)).