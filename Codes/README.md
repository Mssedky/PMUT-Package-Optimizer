# PMUT Acoustic Horn Simulation & Optimization Pipeline

This repository contains four MATLAB scripts that together form a design pipeline
for an acoustic horn/housing that boosts the on-axis pressure of a 0.7 mm PMUT
operating at 180 kHz. They use [k-Wave](http://www.k-wave.org/) for the acoustic
physics and a custom genetic algorithm (GA) for shape optimization.

| File | Role |
|---|---|
| `Solver_Acoustic.m` | Single-design validation/visualization script. Also exports STL and DXF geometry files. |
| `acoustic_cost.m` | Cost function used by the GA — runs one k-Wave simulation for a given design vector and returns a scalar cost. |
| `optimize_acoustic_housing.m` | Genetic algorithm driver that repeatedly calls `acoustic_cost.m` to search for the best horn design. |
| `PMUT_Horn_GUI.m` | Interactive slider-based 3D viewer for a horn shape, with STL export. Uses a **different geometry parameterization** than the other three files (see [Important note](#important-note-two-geometry-parameterizations) below). |

All four scripts model the same physical setup: a bare PMUT membrane sits at the
bottom of a rigid cylindrical housing, and a horn-shaped hole cut into that housing
guides/focuses the emitted acoustic wave. The scripts differ in what they *do* with
that geometry (simulate, optimize, or visualize).

---

## 1. Dependencies

- **MATLAB** (R2018b or newer recommended, for built-in `stlwrite`/`triangulation` support).
- **k-Wave toolbox** (`kWaveGrid`, `kspaceFirstOrder3D`, `kspaceFirstOrder3DG`) — must be on your MATLAB path.
- **Parallel Computing Toolbox** — required for GPU acceleration (`gpuArray-single`, `gpuDevice`). An NVIDIA CUDA-capable GPU is used automatically if detected; otherwise the scripts fall back to CPU (much slower).
- No toolboxes are required for `PMUT_Horn_GUI.m` beyond base MATLAB graphics/UI.

Each script is self-contained (no shared `.m` helper files besides k-Wave itself),
so you can run any one of them independently as long as k-Wave is on the path.

---

## 2. Shared Physical Model

All simulation scripts (`Solver_Acoustic.m`, `acoustic_cost.m`) share the same setup:

- **PMUT**: 0.7 mm diameter circular membrane, driven at `f0 = 180 kHz` with a
  clamped-membrane (0,1)-mode Bessel `J0` radial velocity profile.
- **Housing**: an 8 mm outer diameter cylinder. A short solid cylinder section sits
  directly above the PMUT (`cylinder_height ≈ 0.475 mm`), followed by the horn —
  a generalized ellipsoidal hole whose radius varies with height as:

  ```
  r(z) = r_base*(1 - t)^(1/p) + r_top * t^(1/q),   t = z / horn_height
  ```

  with optional azimuthal/x/y/z sinusoidal modulations (`A_theta`, `A_x`, `A_y`,
  `A_z` and their mode numbers) layered on top for non-axisymmetric shapes.
- **Medium**: air (`c = 343 m/s`, `ρ = 1.21 kg/m³`) with the horn walls modeled as
  a dense, faster material (`c = 2730 m/s`, `ρ = 1190 kg/m³`, e.g. aluminum) acting
  as an effectively rigid boundary.
- **Grid**: a 5.5 × 5.5 × 5 mm domain discretized at 30 points per wavelength
  (~0.064 mm resolution at 180 kHz), with a 5-cell PML absorbing boundary.
- **Metric of interest**: peak pressure captured on an observation plane 4 mm above
  the horn exit, compared against a baseline (bare PMUT, no horn) pressure.

---

## 3. `Solver_Acoustic.m` — Single-Design Validation & Export

**Purpose:** Run one fixed horn design through k-Wave, visualize the results in
detail, and export manufacturable geometry files (STL solid model + DXF profile
curve). This is also how you generate the "no-horn" baseline pressure value that
feeds into `acoustic_cost.m` (`base_amplitude`) — just set `makeHorn = 0` and read
off the reported pressure.

**How to run:**
1. Open the script and edit the parameters under `%% HORN OPTIMIZATION PARAMETERS`
   (`horn_height`, `horn_shape.r_top/p/q/r`, and the azimuthal/x/y/z modulation
   terms) to the design you want to inspect — typically the best design output by
   `optimize_acoustic_housing.m`.
2. Run the script (no function call needed — it's a plain script, not a function).
3. It will, in order:
   - Print grid size, resolution, and CFL number sanity checks.
   - Build the horn geometry and voxelize it into `horn_mask`.
   - Export **`figure_1_acoustic_horn_mm.stl`** — a triangulated 3D surface of the
     horn interior surface, in millimeters, suitable for CAD/3D printing.
   - Export **`figure_1_horn_profile_origin_at_base.dxf`** — the axisymmetric 2D
     radius-vs-height profile curve (origin placed at the horn base), useful for
     lathe/CNC or as a CAD sketch reference.
   - Build the PMUT velocity source and run the k-Wave simulation (GPU if available).
   - Print numeric results: max pressure at the top and bottom observation planes,
     dB re 20 µPa, on-axis pressure, and the peak location.
   - Produce four figure windows:
     - **Figure 1**: pressure maps at both observation heights, the analytical horn
       profile, and an xz cross-section of the pressure field.
     - **Figure 2**: voxelized material map (air vs. horn) in the xz plane, with the
       PMUT and observation plane locations marked.
     - **Figure 3**: the PMUT's Bessel J0 velocity profile (2D map + radial slice).
     - **Figure 4**: a polished xz pressure plot with the horn/cylinder outline
       overlaid, for use directly in a paper or presentation figure.

**Key outputs to look for:** the printed `Maximum pressure up` value (Pa) and the
`Amplitude improvement (top/bottom)` ratio — this is the same figure of merit the
optimizer maximizes.

---

## 4. `acoustic_cost.m` — GA Cost Function

**Purpose:** Given an 11-element design vector `lambda`, build the horn geometry,
run one k-Wave simulation, and return a scalar cost (negative pressure-improvement
ratio, since the GA driver minimizes cost). This file is not meant to be run
directly — it's called by `optimize_acoustic_housing.m` inside the GA loop, but you
can also call it manually to spot-check a single design:

```matlab
lambda = [1.5e-3, 1.8e-3, 4.9, 6.5, 2.4, 0, 0, 0, 0, 0, 0]; % example
cost = acoustic_cost(lambda);
```

**Design vector `lambda` (11 elements, in order):**

| Index | Variable | Meaning |
|---|---|---|
| 1 | `horn_height` | Height of the horn (m) |
| 2 | `horn_shape.r_top` | Radius at the horn exit (m) |
| 3 | `horn_shape.p` | Ellipsoid exponent controlling the base-side curvature |
| 4 | `horn_shape.q` | Ellipsoid exponent controlling the top-side curvature |
| 5 | `horn_shape.r` | Exponent applied to the normalized height in the modulation scale function |
| 6 | `horn_shape.A_theta` | Azimuthal modulation amplitude |
| 7 | `horn_shape.m_theta` | Azimuthal modulation mode number |
| 8 | `horn_shape.A_x` | X-direction modulation amplitude |
| 9 | `horn_shape.A_y` | Y-direction modulation amplitude |
| 10 | `horn_shape.A_z` | Z-direction modulation amplitude |
| 11 | `horn_shape.m_z` | Z-direction modulation mode number |

**How the cost is computed:**
1. Build the horn/cylinder mask exactly as in `Solver_Acoustic.m`, from the given
   `lambda`.
2. Position the PMUT source, define the pressure sensor plane 4 mm above the horn,
   and run `kspaceFirstOrder3D` on GPU.
3. Compute `improvement_ratio = max_pressure / base_amplitude`, where
   `base_amplitude = 106.77 µPa` is the pre-measured no-horn baseline at 4 mm
   (update this constant if your baseline geometry/frequency changes — regenerate
   it with `Solver_Acoustic.m` and `makeHorn = 0`).
4. Return `cost = -improvement_ratio` (so the GA, which minimizes, effectively
   maximizes pressure improvement).
5. **Stability guard:** if `improvement_ratio` exceeds 5000, the design is treated
   as numerically unstable/unphysical and assigned a large penalty cost
   (`1e7`) instead, to keep the GA from chasing simulation artifacts.

**Tuning notes:**
- `Nsteps = 5000` and `deltaT = 4e-9 s` control simulation length/resolution; both
  affect runtime and are shared with `Solver_Acoustic.m`.
- The function auto-computes an even-sized grid from `ppw = 30` points per
  wavelength at 180 kHz — don't change `f0` without also reconsidering `ppw`,
  `Nsteps`, and `deltaT` for CFL stability.
- `input_args` currently hardcodes `'DataCast', 'gpuArray-single'` — if you don't
  have a CUDA GPU, remove that option and expect a large slowdown (unlike
  `Solver_Acoustic.m`, this file does **not** auto-detect GPU availability).

---

## 5. `optimize_acoustic_housing.m` — Genetic Algorithm Driver

**Purpose:** Search the 11-parameter design space for the horn geometry that
maximizes pressure improvement, by repeatedly calling `acoustic_cost.m`.

**How to run:**
```matlab
S1 = optimize_acoustic_housing();
```
No input arguments — all GA settings and parameter bounds are hardcoded at the top
of the function under `%% Givens`. It clears the workspace (`clear all`) on entry,
so save anything you need first.

**GA settings:**
- `S = 30` — population size (designs per generation).
- `G = 100` — number of generations.
- `K = 10` — number of elite parents carried forward each generation.
- `dv = 11` — dimensionality of the design vector (matches `acoustic_cost.m`).
- Parameter bounds are set via `*_min`/`*_max` pairs for each of the 11 variables.
  **Note:** in the current configuration, `A_theta`, `m_theta`, `A_x`, `A_y`,
  `A_z`, and `m_z` all have `min == max == 0`, meaning the search is effectively
  over just the first 5 parameters (`horn_height`, `r_top`, `p`, `q`, `r`) —
  i.e., an axisymmetric horn. Widen those bounds if you want the GA to explore
  non-axisymmetric shapes.

**Algorithm flow:**
1. Randomly initialize a population of 30 design strings within the given bounds.
2. Evaluate `acoustic_cost.m` for every string in generation 1; sort by cost.
3. For each subsequent generation:
   - Take the top `K = 10` elite parents (costs are reused/cached, not
     re-simulated, since they haven't changed).
   - Breed `K` children via arithmetic crossover (`kids = φ·parent_i + (1-φ)·parent_j`
     pairing parents 1↔2, 3↔4, …) with a random blend factor `φ` per pair.
   - Fill the remaining `S - 2K = 10` slots with fresh random strings (mutation via
     resampling, not perturbation).
   - Evaluate cost only for the new children + random strings (10 new simulations
     per generation, not 30 — this is what keeps runtime manageable).
   - Sort, log, and repeat.
4. Track best/average cost per generation and stop after `G` generations (there is
   also a commented-out early-stopping tolerance `TOL`, currently disabled).

**Outputs:**
- **`ga_convergence_log.txt`** — tab-separated `Generation \t MinCost` log, written
  incrementally as the GA runs (useful for monitoring long runs or resuming later
  with your own logic).
- **Figure 1** — semilog plot of best cost vs. generation (convergence curve).
- **Console output** — for every simulated design, the full parameter vector and
  resulting cost are printed; at the end, the top 4 performing design vectors and
  their costs are printed and returned.
- **`S1`** (function return value) — the single best design's parameter vector
  (`horn_height`, `r_top`, `p`, `q`, `r`, and the modulation terms). Feed this
  directly into `Solver_Acoustic.m`'s parameter block to validate, visualize, and
  export the winning design.

**Runtime:** with `S = 30`, `G = 100`, ~10 new simulations per generation after
gen 1, expect roughly `1 + 99×10 ≈ 991` total k-Wave simulations. On an NVIDIA RTX
2060, this pipeline has taken on the order of ~15 hours end-to-end — budget
accordingly, and consider reducing `G` or `S` for faster iteration during early
exploration.

---

## 6. `PMUT_Horn_GUI.m` — Interactive Horn Viewer & STL Export

**Purpose:** A standalone, non-simulation tool for interactively sculpting and
inspecting a horn's 3D shape in real time, then exporting it as an STL. This is
useful for building physical intuition about how each parameter changes the
geometry before committing to a full k-Wave run, or for quickly exporting a
manufacturable shape once you already know the parameters you want.

**How to run:**
```matlab
pMUT_Horn_GUI
```
This opens a figure window with a 3D rendered horn on the left and a parameter
panel with sliders (plus editable numeric text boxes) on the right.

**Adjustable parameters:**
- `r_base`, `r_top`, `height`, `p`, `q` — same base-shape controls as the ellipsoid
  radius profile used elsewhere in the pipeline.
- `dx_top`, `dy_top` — lateral (x/y) offset of the horn's *exit* centerline
  relative to its base, letting you tilt/shift the horn off-axis.
- `exit_tilt_x`, `exit_tilt_y` — tilt angles (radians) of the exit-plane normal
  relative to the base normal.
- `tilt_pow`, `shift_pow` — exponents controlling how quickly the tilt and lateral
  shift ramp in from base to exit (higher = more of the change concentrated near
  the top).

Dragging a slider (or typing a value into its adjacent edit box) regenerates the
surface live via `generateHorn`, which builds the shape by sweeping a circular
cross-section of radius `r(z)` along a smoothly tilting/shifting centerline and
orthonormal frame — a fundamentally different construction from the
azimuthal-modulation approach in `acoustic_cost.m`/`Solver_Acoustic.m` (see below).

**STL export:** click **"Export STL"** to save the currently displayed shape
(scaled from meters to millimeters) via a standard MATLAB file-save dialog.

---

## Important note: two geometry parameterizations

`PMUT_Horn_GUI.m` uses a **tilt/shift-based** horn parameterization
(`dx_top`, `dy_top`, `exit_tilt_x`, `exit_tilt_y`, `shift_pow`, `tilt_pow`) that
bends and offsets the horn's centerline from base to exit.

`acoustic_cost.m`, `optimize_acoustic_housing.m`, and `Solver_Acoustic.m` instead
use an **azimuthal/axial modulation** parameterization (`A_theta`, `m_theta`,
`A_x`, `A_y`, `A_z`, `m_z`) that keeps the centerline fixed on-axis and instead
perturbs the cross-sectional radius as a function of angle and height.

These two representations are **not currently interchangeable** — a design vector
from the GA (`optimize_acoustic_housing.m`) cannot be typed directly into the GUI's
parameter panel, and vice versa. If you want to preview a GA-optimized design in
the GUI's live 3D viewer, or run a GUI-tilted design through k-Wave, you'll need to
either (a) restrict yourself to the shared parameters (`r_base`/`r_top`/`height`/
`p`/`q`, with all tilt/shift/azimuthal terms zeroed out — this is effectively what
the current GA bounds already do), or (b) port one geometry generator's math into
the other script.

---

## Suggested end-to-end workflow

1. **Baseline:** Run `Solver_Acoustic.m` with `makeHorn = 0` to get the bare-PMUT
   pressure at 4 mm; update `base_amplitude` in `acoustic_cost.m` if it differs
   from the current hardcoded value.
2. **Optimize:** Run `optimize_acoustic_housing.m` (budget several hours on GPU).
   Monitor `ga_convergence_log.txt` and the printed per-generation best cost.
3. **Validate:** Take the returned `S1` (or the top design printed at the end) and
   plug its values into the `HORN OPTIMIZATION PARAMETERS` block in
   `Solver_Acoustic.m`. Run it to get detailed pressure maps, the improvement
   ratio, and exported STL/DXF geometry files.
4. **Inspect/refine geometry:** Optionally, translate the winning shape parameters
   into `PMUT_Horn_GUI.m` (base/top radius, height, p, q — leave tilt/shift/azimuthal
   terms at 0) to interactively inspect the 3D form and re-export a clean STL for
   fabrication.
