# Headless Validation Artifacts (`crvscience-itme`)

Generated directly by running the compiled Rust integration test suite:

```bash
cargo test --test headless_validation
```

All outputs are kept under `384 KB` (`393,216` bytes) per file and capture genuine runtime execution of `itme_evolution_sim::engine::SimulationWorld` (`v2.3.0`):

1. **`planetary_epochs_trajectory.csv`** (`T1`): Continuous 4.25 Ga 8-epoch Hadean-to-Gaia trajectory (`run_continuous_grand_trajectory(12)`) across `4.2500 -> 0.0000 Ga`, tracking `population_size`, `heritable_channels_active`, `heritable_mutual_info_bits`, `informational_leverage_li`, `order_velocity_vo`, atmospheric `planetary_o2`, `planetary_co2_ppm`, and `planetary_exergy_budget`.
2. **`emergent_biosphere_tree.nwk`** (`T2`): Balanced Newick tree (`export_newick(false)`) exported after 40 baseline ticks, a combined environmental shock + apex predator cull (`trigger_environmental_shock` + `trigger_apex_predator_cull`), and 50 recovery ticks.
3. **`itme_validation_summary.json`** (`T1 + T2 + T3`): Machine-readable telemetry record including T1 epoch telemetry, T2 shock-recovery & Sackin imbalance metrics, and T3 multi-seed (`[101, 202, 303]`) geothermal flux (`bath_exergy_flux`) and `heritable_variation` ensemble summaries with bit-exact determinism verification.
