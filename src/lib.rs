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
