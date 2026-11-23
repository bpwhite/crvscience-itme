# ITME v23.0 (`v2.3.0`) User Guide

**Information-Theoretic Model of Evolution (`itme_evolution_sim` v2.3.0)** is an interactive native Rust desktop simulator and headless benchmark engine distributed under the **GNU General Public License v3.0 or later (`GPL-3.0-or-later`)**.

---

## 1. Launching the Application

### Native Desktop GUI

```bash
cargo run --release
```

Force the GUI explicitly when passing other options:

```bash
cargo run --release -- --gui --seed 20250813
```

### Headless & Feature-Stripped Builds

If you are running on a server or CI node without graphics libraries, compile or run without the `gui` feature:

```bash
cargo check --no-default-features --all-targets
cargo run --no-default-features --release -- --headless --grand-trajectory
```

---

## 2. Four-Card Main Menu Launcher

Upon opening the GUI, ITME starts at the **4-Card Main Menu**:

| Card | Domain Scope | Sandbox Isolation & Key Mechanics |
| --- | --- | --- |
| **Domain 1 · Molecular Evolution** | `RNA/DNA/Protein` | Isolated sandbox covering Hadean `FeNiS` / wet-dry protocells, `Fe3O4` CISS homochirality, 20-AA `SGC`, Class III `RNR` `T-DNA`, Eigen's fidelity cliff, `1°/2°/3°` endosymbiosis, and `Mito -> Hydro -> MitoS -> MRO->0`. |
| **Domain 2 · Hominid Populations** | `HominidPopulation` | Isolated sandbox covering foraging band fissioning, territorial catchment competition, symbolic teaching ratchets, agro-ectosymbiosis, and external written archives with peer verification. |
| **Domain 3 · Artificial Intelligences** | `ArtificialIntelligence` | Isolated sandbox covering compute cluster load-balancing, container forking, online test-time RL, synthetic model collapse (Shumailov cliff), teacher distillation, and formal proof verifiers. |
| **Domain 4 · Coupled Planetary Biosphere** | `4.25 Ga PlanetaryBiosphere` | Couples autotrophic, syntrophic, superorganismal, and machine-verified guilds to live atmospheric `O2 (PAL Index)`, `CO2 (ppm)`, `N/P` nutrient pool, Gibbs chemical disequilibrium, and the 8-Epoch Grand Trajectory. |

### Launch Footer Options

- **Mode**: `13-Regime Guided Tour` (cycles automatically across all 13 emergent regimes) or `Autonomous Sandbox` (stays in manual/emergent exploration).
- **Initial Regime**: `R(0.0) Basal`, `R(3.0) Darwinian`, `R(3.0) Soft/Distill`, or `R(3.5) Integrated`.
- **Seed**: Deterministic Xoshiro256** PRNG seed accepting decimal integers (e.g., `20250813`) or `0x`-prefixed hexadecimal (`u64`). Invalid inputs display an inline warning and safely fall back to the active valid seed.
- **Ticks/Regime**: Slider (`12..=80` on the launcher, `12..=120` in-session) controlling stage dwell time during Guided Tour mode.

---

## 3. Global Keyboard Shortcuts & Mouse Controls

Global hotkeys are active whenever no text input field has keyboard focus (`!ctx.wants_keyboard_input()`):

| Key / Input | Action | Notes |
| --- | --- | --- |
| `Space` | Play / Pause simulation | Resets step timer cleanly on resume |
| `.` (Period) | Single Step | Advances simulation by 1 discrete tick |
| `R` | Reset 3D Camera | Restores `Camera3D::default()` (`yaw=0.58`, `pitch=0.34`, `zoom=1.0`) |
| `1` | Spatial + Network View | 2D reaction-diffusion spatial lattice + reticulate phylogram |
| `2` | Spatial + 3D View | 2D spatial lattice + interactive 3D structure |
| `3` | Full Lineage View | Full-width anastomosing phylogeny / distillation graph |
| `4` | 3D Structure View | Depth-interleaved 3D nodes, glow halos, dashed bonds/links |
| `5` | 2D Blueprint View | 2D schematic blueprint with directed edge arrowheads |
| `6` | 4.25 Ga Timeline View | 8-Epoch Grand Trajectory charts & table (`PlanetaryBiosphere` mode) |
| `F1` | Toggle About Modal | Displays version, license (`GPL-3.0-or-later`), and CRV Science links |
| `Esc` | Close Modal / Main Menu | Closes the About modal if open; otherwise returns to Main Menu |
| **Mouse Drag** | Orbit 3D Camera | Drag inside the central viewport to rotate yaw and pitch |
| **Scroll Wheel** | Zoom 3D Camera | Smoothly zooms between `0.45x` and `2.8x` |

---

## 4. Top Toolbar, Export Dropdown & Collapsible Milestones

- **Simulation Controls**: `<- Main Menu`, `Play`/`Pause` (`Space`), `Step` (`.`), and the `ticks/s` speed slider (`1.0..=30.0`).
- **Compact `Export...` Dropdown**: Consolidates all four exporters into a single menu:
  - **Telemetry CSV (`.csv`)**: Exports per-tick `ItmeMetrics` (or the 8-epoch Grand Trajectory table when viewing `4.25 Ga Timeline`).
  - **Stage Atlas SVG (`.svg`)**: Renders a 13-regime `3D/2D` vector overview across the active domain.
  - **Extended Newick (`.nhx`)**: Exports `#H` reticulate hybrid tags and `MRO/APICO/KLEPTO` NHX annotations.
  - **Phylogram SVG (`.svg`)**: Exports the anastomosing evolutionary phylogram as vector SVG.
- **View Switcher (`1`–`6`) & Render Toggles**: Switch layouts, toggle `Auto-Orbit`, toggle `Labels` (pill backdrops with overlap culling), and click `Milestones [-]` / `Milestones [+]` to collapse or expand the quick pulse toolbar so the viewport and telemetry plots receive maximum vertical space.
- **Collapsible Milestone Pulse Row**:
  - *PlanetaryBiosphere*: `4.25 Ga Grand Trajectory`, `Gaia-2.0 Awakening`, `Extractive CO2 Overshoot`.
  - *Molecular*: `Hadean Abiogenesis -> LUCA`, `Homochiral D/L + SGC`, `RNR RNA->T-DNA`, `1°/2° Symbio + MRO`.
  - *Hominid*: `Kin Division & Lithic Guilds`, `Symbolic Teaching + Ratchet`, `Agro-Ectosymbiosis & Crops`, `Written Archive + Peer Audit`.
  - *Artificial Intelligence*: `Container Checkpoint Scale`, `Adversarial Arena (PBT)`, `Formal Proof Verifier`, `Anti-Collapse Safeguard`.

---

## 5. Viewport Rendering & Bottom Telemetry Plots

- **3D Structure View (`4`, `2`)**: Draws depth-interleaved `nodes_3d` and `edges_3d`, sphere shading with atmospheric depth fog, `Node3D::glow` halos, true dashed segments (`Edge3D::dashed`), and hover node inspection.
- **2D Blueprint View (`5`)**: Renders directed graph edges (`Edge2D::directed`) with clean arrowheads terminating at target node boundaries plus hover tooltips.
- **Bottom Telemetry Plots**:
  - When fewer than `2` ticks exist (`history.len() < 2`), an onboarding prompt instructs you to press `Play` (`Space`) or `Step` (`.`).
  - When `history.len() >= 2`, both panels display horizontal and vertical gridlines, Y-axis tick labels, X-axis `Ga` time labels (`4.25 Ga -> 0.00 Ga`), and synchronized hover crosshairs mapping pointer `(x, y)` to `Ga`, tick `step`, cursor Y value, and exact series values (`Order O` on the left; `O2 (PAL Index)`, `N/P Pool`, and `Gibbs/26` or domain-specific series on the right).

---

## 6. Headless CLI & Benchmark Flags

Run `cargo run --release -- --help` to see all supported flags:

| Flag | Description | Default |
| --- | --- | --- |
| `--headless`, `--console` | Run terminal console simulation without launching a window | Off |
| `--grand-trajectory` | Run the continuous 4.25 Ga 8-epoch grand trajectory benchmark (implies `--headless`) | Off |
| `--gui` | Force native desktop GUI mode | On (when `gui` feature enabled) |
| `--screenshot <FILE>` | Render a `1600x1000` offscreen `wgpu` PNG screenshot and exit | None |
| `--screenshot-view <1\|2\|3\|4>` | Select preset (`1=Main Menu`, `2=Molecular`, `3=AI`, `4=Biosphere`) | `1` |
| `--organism <rna\|hominid\|ai\|all>` | Filter target organism substrate(s) | `all` |
| `--regime <CODE\|all>` | Run `all` 13 regimes or one code (`0.0`, `0.5`, `1.0`, `1.5-AA`, `1.5-RP`, `2.0-EC`, `2.0-BAL`, `2.0-QS`, `3.0-DAR`, `3.0-CUL`, `3.5-SYM`, `4.0-ENG`, `4.0-AUT`) | `all` |
| `--steps <N>` | Discrete simulation ticks per regime (or per epoch in `--grand-trajectory`) | `24` |
| `--seed <U64>` | Deterministic Xoshiro256** PRNG seed | `20250813` |
| `--no-ascii-3d` | Omit inline 3D ASCII & spatial map projections from console output | Show ASCII |
| `--export-csv <FILE>` | Export step-by-step telemetry (or grand trajectory summary) CSV | None |
| `--export-svg <FILE>` | Export 13-regime 3D/2D structural stage atlas SVG | None |
| `--export-newick <FILE>` | Export Extended Newick (`eNewick #H` + `MRO/APICO/KLEPTO` NHX) | None |
| `--export-phylo-svg <FILE>` | Export anastomosing SVG phylogram | None |
| `-V`, `--version` | Print package version (`itme_evolution_sim 2.3.0`) and exit | — |
| `-h`, `--help` | Print CLI help summary and exit | — |

### Common Headless Workflows

```bash
# Continuous 8-epoch grand trajectory with eNewick + phylogram + CSV export
cargo run --release -- --headless --grand-trajectory --steps 12 --export-csv grand.csv --export-newick grand.nhx --export-phylo-svg grand.svg

# Single-domain high-speed sweep (Molecular R(3.5-SYM) for 40 ticks)
cargo run --release -- --headless --organism rna --regime 3.5-SYM --steps 40 --no-ascii-3d

# Offscreen GPU screenshot of the Coupled Planetary Biosphere view
cargo run --release -- --screenshot biosphere.png --screenshot-view 4
```
