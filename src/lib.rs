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
    use crate::model::{compute_adami_conservation, lz76_complexity, PriceDecomposition};
    use crate::regimes::{classify_emergent_state, RegimeId};

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
