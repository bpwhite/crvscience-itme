// Copyright (C) 2026 CRV Science Network (https://crvscience.com)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://gnu.org>.

//! Information-Theoretic Model of Evolution (ITME v23.0 / `v2.3.0`)
//! **Continuous 4.25-Billion-Year Multi-Epoch Grand Trajectory Engine**:
//! Unbroken simulation chaining `Epoch 1: Homochirality & Minihelix L-AA` -> `Epoch 2: 20-AA SGC & RAF Protocells` ->
//! `Epoch 3: RNR RNA-to-DNA, ThyX/UDG, Virosphere & LUCA` -> `Epoch 4: GOE & Asgard-Mito Symbiogenesis (LECA)` ->
//! `Epoch 5: Clonal Germ-Soma Multicellularity` -> `Epoch 6: Eusocial Superorganisms & Clonal Agriculture` ->
//! `Epoch 7: Symbolic Culture, Peer-Verified Archives & Machine-Verified Synthesis` ->
//! `Epoch 8: Gaia-2.0 Planetary Biogeochemical Homeostasis`!

pub mod console;
pub mod engine;
pub mod geometry;
#[cfg(feature = "gui")]
pub mod gui;
pub mod model;
pub mod organisms;
pub mod phylogeny;
pub mod regimes;
pub mod spatial;

#[cfg(test)]
mod tests {
    use crate::model::{
        compute_adami_conservation, compute_critical_fidelity, lz76_complexity,
        EnvironmentControls, PlanetaryBiosphere, PriceDecomposition,
    };
    use crate::regimes::{classify_emergent_state, RegimeId};
    use crate::spatial::SpatialGrid;

    #[test]
    fn unit_test_adami_conservation_clonal_vs_polymorphic() {
        let a = [1u16, 4, 7, 2, 9, 3, 5, 8];
        let b = [1u16, 4, 7, 2, 9, 3, 5, 8];
        let c = [0u16, 2, 3, 5, 1, 6, 4, 0];
        let clonal = compute_adami_conservation(&[&a, &b], 12);
        let poly = compute_adami_conservation(&[&a, &c], 12);
        assert!((clonal - 1.0).abs() < 1e-12);
        assert!(poly < 1.0);
    }

    #[test]
    fn unit_test_adami_conservation_log2_alphabet_normalization_and_edge_cases() {
        let rho = compute_adami_conservation(&[&[0, 1], &[0, 2]], 20);
        assert!(
            rho > 0.75 && rho < 1.0,
            "Adami log2(20) site normalization must yield > 0.75, got {rho}"
        );
        let expected = 1.0 - (0.5 / 20.0_f64.log2());
        assert!((rho - expected).abs() < 1e-12);

        // Edge cases must remain completely safe and finite
        assert!((compute_adami_conservation(&[], 20) - 1.0).abs() < 1e-12);
        assert!((compute_adami_conservation(&[&[0, 1]], 20) - 1.0).abs() < 1e-12);
        assert!((compute_adami_conservation(&[&[], &[]], 20) - 1.0).abs() < 1e-12);
        assert!((compute_adami_conservation(&[&[0], &[]], 20) - 1.0).abs() < 1e-12);
        assert!((compute_adami_conservation(&[&[0, 1], &[0, 2]], 1) - 1.0).abs() < 1e-12);
        assert!((compute_adami_conservation(&[&[0, 1], &[0, 2]], 0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn unit_test_critical_fidelity_upper_clamp_0999_and_o2_pal_index_bounds() {
        let f_high = compute_critical_fidelity(2500.0, 2.0);
        assert!(
            f_high > 0.95 && f_high <= 0.999,
            "Expected Eq. 11 threshold in (0.95, 0.999], got {f_high}"
        );
        let f_sat = compute_critical_fidelity(1.0e7, 2.0);
        assert!((f_sat - 0.999).abs() < 1e-12, "Expected upper clamp 0.999, got {f_sat}");
        let f_low = compute_critical_fidelity(2.0, 20.0);
        assert!((f_low - 0.72).abs() < 1e-12, "Expected lower clamp 0.72, got {f_low}");
        assert!(compute_critical_fidelity(f64::NAN, f64::NAN).is_finite());

        let mut ctrl = EnvironmentControls::default();
        ctrl.biogeochemical_feedback_gain = 1.0;
        ctrl.planetary_stewardship_coupling = 1.0;
        let mut bio = PlanetaryBiosphere::new_for_controls(&ctrl);
        for _ in 0..80 {
            bio.step_biogeochemical_feedbacks(3.0, 2.5, 2.5, 0.0, 2.5, &ctrl);
            assert!(
                (0.02..=0.94).contains(&bio.atmospheric_o2),
                "O2 (PAL Index) escaped [0.02, 0.94]: {}",
                bio.atmospheric_o2
            );
        }
        ctrl.extractive_burn_intensity = 1.0;
        ctrl.planetary_stewardship_coupling = 0.0;
        for _ in 0..80 {
            bio.step_biogeochemical_feedbacks(0.0, 0.0, 0.0, 3.0, 0.0, &ctrl);
            assert!(
                (0.02..=0.94).contains(&bio.atmospheric_o2),
                "O2 (PAL Index) escaped [0.02, 0.94]: {}",
                bio.atmospheric_o2
            );
        }
    }

    #[test]
    fn unit_test_spatial_toroidal_mass_conservation_and_nonfinite_flux_guard() {
        let mut grid = SpatialGrid::new(10, 8, 1.0);
        assert_eq!(grid.idx(10, 8), grid.idx(0, 0));
        assert_eq!(grid.idx(19, 15), grid.idx(9, 7));

        grid.byproduct.fill(0.0);
        grid.niche.fill(0.0);
        let corner = grid.idx(0, 0);
        grid.byproduct[corner] = 1.6;
        grid.niche[corner] = 1.8;

        grid.step_reaction_diffusion(1.0, 0.12, 0.02, 0.10);

        // Verify seamless periodic wrap-around across (0, 0) <-> (width - 1, 0) and (0, height - 1)
        assert!(grid.local_byproduct(9, 0) > 0.0);
        assert!(grid.local_byproduct(0, 7) > 0.0);
        assert!((grid.local_byproduct(9, 0) - grid.local_byproduct(1, 0)).abs() < 1e-12);
        assert!((grid.local_byproduct(0, 7) - grid.local_byproduct(0, 1)).abs() < 1e-12);
        assert!((grid.local_niche(9, 0) - grid.local_niche(1, 0)).abs() < 1e-12);
        assert!((grid.local_niche(0, 7) - grid.local_niche(0, 1)).abs() < 1e-12);

        // Verify exact mass conservation under toroidal diffusion + linear decay
        let total_b: f64 = grid.byproduct.iter().sum();
        let total_n: f64 = grid.niche.iter().sum();
        assert!(
            (total_b - 1.6 * (1.0 - 0.085)).abs() < 1e-12,
            "Toroidal byproduct mass mismatch: {total_b}"
        );
        assert!(
            (total_n - 1.8 * (1.0 - 0.02)).abs() < 1e-12,
            "Toroidal niche mass mismatch: {total_n}"
        );

        // Verify NaN / Inf / -Inf bath_exergy_flux cannot corrupt fields
        let prev_exergy = grid.exergy.clone();
        let prev_byproduct = grid.byproduct.clone();
        let prev_niche = grid.niche.clone();
        for bad_flux in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            grid.step_reaction_diffusion(bad_flux, 0.12, 0.02, 0.10);
            assert!(grid.exergy.iter().all(|v| v.is_finite()));
            assert!(grid.byproduct.iter().all(|v| v.is_finite()));
            assert!(grid.niche.iter().all(|v| v.is_finite()));
            assert_eq!(grid.exergy, prev_exergy);
            assert_eq!(grid.byproduct, prev_byproduct);
            assert_eq!(grid.niche, prev_niche);
        }
    }

    #[test]
    fn unit_test_lz76_nonperiodic_archival_vs_periodic_wrap() {
        let wrapped16: Vec<u16> = (0..38u16).map(|k| (11 + k * 7) % 16).collect();
        let nonperiodic64: Vec<u16> = (0..38u16).map(|k| (11 + k * 7) % 64).collect();
        let c16 = lz76_complexity(&wrapped16);
        let c64 = lz76_complexity(&nonperiodic64);
        assert!(c16 <= 18, "expected periodic wrap c16 <= 18, got {c16}");
        assert!(c64 > 30, "expected nonperiodic archival c64 > 30, got {c64}");
        assert!(c64 > c16);
    }

    #[test]
    fn unit_test_clonal_replicator_ignores_spurious_price_noise() {
        let price = PriceDecomposition {
            selection_covariance: 0.015,
            transmission_bias: 0.0,
            selection_share: 1.0,
        };
        let (_, _, _, regime) = classify_emergent_state(
            0.82,
            0.50,
            0.25,
            0.05,
            0.30,
            0.0,
            1.0,
            0.995,
            0.75,
            0.0,
            0.0,
            &price,
            false,
        );
        assert_eq!(regime, RegimeId::R15Replicative);
    }
}
