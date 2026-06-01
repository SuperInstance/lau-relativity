# lau-relativity

A pure-Rust library for special and general relativity — spacetime geometry, Lorentz transformations, relativistic kinematics, energy-momentum, tensor formulation, geodesic integration, gravitational redshift, and cosmology.

**87 tests** · `nalgebra` + `serde` + `num-complex` · MIT license

---

## What This Does

`lau-relativity` implements the core mathematical structures of Einstein's relativity as composable Rust modules. From the Minkowski metric and Lorentz boosts of special relativity through the Schwarzschild metric, Christoffel symbols, and geodesic integration of general relativity, to FLRW cosmology and ΛCDM parameters.

Every function is pure (no side effects, no I/O), uses SI units throughout, and is backed by `nalgebra` for the linear algebra. Test coverage includes interval invariance under boosts, energy-momentum relation verification, and Christoffel symbol symmetry.

---

## Key Idea

Relativity is fundamentally about **coordinate transformations** and **invariant quantities**. This crate mirrors that structure:

| Module | Domain | Key Structures |
|---|---|---|
| `minkowski` | Flat spacetime | 4-vectors, metric η_μν, spacetime intervals |
| `lorentz` | Boosts & rotations | Boost matrices (x/y/z), velocity addition, rapidity |
| `kinematics` | SR kinematics | Time dilation, length contraction, Doppler effect |
| `energy_momentum` | Relativistic mechanics | 4-momentum, E² = p²c² + m²c⁴, kinetic energy |
| `tensor` | Curved spacetime | Schwarzschild metric, Christoffel symbols, FLRW metric |
| `geodesic` | Geodesic motion | Effective potential, RK4 integration, ISCO, photon sphere |
| `redshift` | Gravitational redshift | Wavelength/frequency shift, Shapiro delay |
| `cosmology` | Cosmology | Friedmann equations, ΛCDM parameters, Hubble law |
| `agent_spacetime` | Multi-agent analogy | Causal structure for distributed agent systems |

---

## Install

```toml
[dependencies]
lau-relativity = "0.1"
```

Or:

```sh
cargo add lau-relativity
```

---

## Quick Start

### Lorentz Boost a 4-Vector

```rust
use lau_relativity::minkowski::FourVector;
use lau_relativity::lorentz::{boost, Axis, gamma};

let event = FourVector::from_spatial(1e-6, 100.0, 0.0, 0.0); // t=1μs, x=100m
println!("Interval s² = {:.3}", event.interval());

let boosted = boost(&event, 0.6, Axis::X);  // β = 0.6c along x
println!("Boosted: ct'={:.3}, x'={:.3}", boosted.ct, boosted.x);
println!("Interval preserved: {:.3}", boosted.interval());  // same as original
```

### Time Dilation

```rust
use lau_relativity::kinematics::{time_dilated, length_contracted};

let proper_time = 1.0;  // 1 second on the ship's clock
let beta = 0.99;        // travelling at 0.99c

let earth_time = time_dilated(proper_time, beta);
println!("1s ship time = {:.3}s Earth time", earth_time);  // ~7.09s

let ship_length = length_contracted(100.0, beta);
println!("100m ship appears {:.1}m from Earth", ship_length);  // ~14.1m
```

### Energy-Momentum

```rust
use lau_relativity::energy_momentum::{FourMomentum, relativistic_kinetic_energy};

let p = FourMomentum::from_rest_mass_and_beta(1.0, 0.8);  // 1 kg at 0.8c
println!("Total energy E = {:.3e} J", p.energy());
println!("Rest mass from invariant: {:.6} kg", p.rest_mass());  // 1.0
println!("E² = (pc)² + (mc²)² holds: {}", p.verify_energy_momentum_relation());

let ke = relativistic_kinetic_energy(1.0, 0.8);
println!("Kinetic energy = {:.3e} J", ke);
```

### Schwarzschild Geodesic

```rust
use lau_relativity::geodesic::{schwarzschild_christoffel, geodesic_step_rk4, isco_radius};
use lau_relativity::tensor::schwarzschild_radius;

let sun_mass = 1.989e30;
let r_s = schwarzschild_radius(sun_mass);  // ~2953 m
println!("Schwarzschild radius of Sun: {:.0} m", r_s);
println!("ISCO: {:.0} m", isco_radius(r_s));           // ~8859 m

// Integrate a geodesic step near a black hole
let pos = [0.0, 1e7, std::f64::consts::FRAC_PI_2, 0.0];
let vel = [1.0, 0.0, 0.0, 0.001];
let (new_pos, new_vel) = geodesic_step_rk4(&pos, &vel, r_s, 100.0);
```

### Cosmology (ΛCDM)

```rust
use lau_relativity::cosmology::LCDMParams;

let params = LCDMParams::planck2018();
println!("Ω_m = {:.3}, Ω_Λ = {:.3}", params.omega_m, params.omega_lambda);
println!("Ω_k = {:.6} (≈ flat)", params.omega_k());
println!("Age of universe ≈ {:.2e} s ({:.1} Gyr)",
    params.age_estimate(),
    params.age_estimate() / 3.156e16);

let h_at_z1 = params.h_at(0.5);  // H at scale factor 0.5 (z=1)
println!("H(z=1) = {:.3e} s⁻¹", h_at_z1);
```

---

## API Reference

### `minkowski` — Minkowski Spacetime

| Type / Function | Description |
|---|---|
| `C: f64` | Speed of light = 299,792,458 m/s |
| `minkowski_metric()` | Returns η_μν = diag(+1, −1, −1, −1) as `Matrix4` |
| `FourVector` | Spacetime 4-vector (ct, x, y, z) |
| `.interval()` | s² = (ct)² − x² − y² − z² |
| `.is_timelike()` / `.is_spacelike()` / `.is_lightlike()` | Interval classification |
| `.proper_time()` | τ = √(s²)/c for timelike intervals |
| `.dot(&other)` | Minkowski inner product |
| `.spatial_norm()` | √(x² + y² + z²) |
| `.time()` | Coordinate time in seconds |

### `lorentz` — Lorentz Transformations

| Function | Description |
|---|---|
| `gamma(β)` | Lorentz factor γ = 1/√(1 − β²) |
| `boost_x(β)` / `boost_y(β)` / `boost_z(β)` | Boost matrices along each axis |
| `boost(&v, β, axis)` | Apply boost to a `FourVector` |
| `rotation_xy(θ)` / `rotation_xz(θ)` | Spatial rotations (4×4 matrix) |
| `velocity_addition(β₁, β₂)` | Relativistic velocity addition |
| `rapidity(β)` / `from_rapidity(φ)` | Rapidity parameter: β = tanh(φ) |

### `kinematics` — Relativistic Kinematics

| Function | Description |
|---|---|
| `time_dilated(τ, β)` | Δt = γτ |
| `proper_time_from_dilated(Δt, β)` | τ = Δt/γ |
| `length_contracted(L₀, β)` | L = L₀/γ |
| `proper_length_from_contracted(L, β)` | L₀ = γL |
| `doppler_factor(β, approaching)` | Relativistic Doppler: √((1±β)/(1∓β)) |
| `ReferenceFrame` | Encapsulates velocity + direction; provides γ, time-dilation & length-contraction factors |

### `energy_momentum` — Relativistic Mechanics

| Type / Function | Description |
|---|---|
| `FourMomentum` | (E/c, pₓ, pᵧ, p_z) |
| `::from_rest_mass_and_beta(m, β)` | Construct from rest mass & velocity |
| `.energy()` | Total energy E |
| `.rest_energy(m)` | E₀ = mc² |
| `.rest_mass()` | From invariant mass: m = √((E/c)² − p²) / c |
| `.kinetic_energy()` | T = E − mc² |
| `.momentum_magnitude()` | \|p\| = √(pₓ² + pᵧ² + p_z²) |
| `.verify_energy_momentum_relation()` | Check E² = (pc)² + (mc²)² |
| `relativistic_kinetic_energy(m, β)` | T = (γ−1)mc² |
| `total_energy(m, β)` | E = γmc² |
| `relativistic_momentum(m, β)` | p = γmβc |
| `beta_from_kinetic_energy(T, m)` | Inverse: velocity from kinetic energy |

### `tensor` — Curved Spacetime

| Function | Description |
|---|---|
| `schwarzschild_radius(M)` | r_s = 2GM/c² |
| `schwarzschild_metric(r, θ, r_s)` | Diagonal metric (g_tt, g_rr, g_θθ, g_φφ) |
| `schwarzschild_christoffel(r, θ, r_s)` | All non-zero Γ^μ_{νλ} as 4×4×4 array |
| `flrw_metric(a, k, r, θ)` | FLRW metric diagonal components |
| `MetricTensor` | General 4×4 metric wrapper with `.inverse_diagonal()` |

### `geodesic` — Geodesic Motion

| Function | Description |
|---|---|
| `effective_potential(r, M, L)` | V_eff = −GM/r + L²/(2r²) − GML²/(c²r³) |
| `geodesic_acceleration(pos, vel, r_s)` | d²x^μ/dτ² = −Γ^μ_{νλ} (dx^ν/dτ)(dx^λ/dτ) |
| `geodesic_step_rk4(pos, vel, r_s, dτ)` | Single RK4 integration step |
| `circular_orbit_radius_newtonian(L, M)` | r = L²/(GM) |
| `isco_radius(r_s)` | ISCO = 3r_s |
| `photon_sphere_radius(r_s)` | Photon sphere = 1.5r_s |
| `SchwarzschildConstants` | Constants of motion (E/m, L/m) |

### `redshift` — Gravitational Redshift

| Function | Description |
|---|---|
| `gravitational_redshift_ratio(r_emit, r_obs, r_s)` | λ_obs/λ_emit |
| `gravitational_redshift_z(r_emit, r_obs, r_s)` | z = Δλ/λ |
| `gravitational_frequency_ratio(r_emit, r_obs, r_s)` | f_obs/f_emit |
| `gravitational_time_dilation(r, r_s)` | dτ/dt = √(1 − r_s/r) |
| `shapiro_delay(r_emit, r_obs, b, M)` | Δt ≈ (4GM/c³)ln(4r₁r₂/b²) |

### `cosmology` — Cosmology

| Function | Description |
|---|---|
| `friedmann_hubble(ρ, k, a, Λ)` | First Friedmann equation: H² = 8πGρ/3 − k/a² + Λ/3 |
| `friedmann_hubble_with_params(H₀, Ω_r, Ω_m, Ω_k, Ω_Λ, a)` | Parameterised version |
| `friedmann_acceleration(ρ, p, Λ)` | Second Friedmann: ä/a = −4πG(ρ+3p/c²)/3 + Λ/3 |
| `critical_density(H)` | ρ_c = 3H²/(8πG) |
| `hubble_distance(H₀)` / `hubble_time(H₀)` | d_H = c/H₀, t_H = 1/H₀ |
| `cosmological_redshift(a)` / `scale_factor_from_redshift(z)` | 1+z = 1/a |
| `luminosity_distance_matter_only(z, H₀)` | d_L for Ω_m=1 flat universe |
| `LCDMParams` | ΛCDM cosmological parameters (H₀, Ω_m, Ω_Λ, Ω_r) |
| `LCDMParams::planck2018()` | Planck 2018 best-fit values |
| `.omega_k()` / `.age_estimate()` / `.h_at(a)` | Derived quantities |

### `agent_spacetime` — Agent Causal Structure

An experimental module applying relativistic analogies to multi-agent systems:

| Type | Description |
|---|---|
| `AgentSpacetime` | Agent with proper time, coordinate time, β, state position |
| `.interval_to(&other)` | Spacetime interval between agents |
| `.can_influence(&other)` | True if timelike or lightlike separation |
| `AgentFrameSystem` | Multi-agent reference frame with Lorentz boosts |
| `.causal_past(idx)` | Agents that can causally influence target |

---

## How It Works

### Minkowski Spacetime

Spacetime 4-vectors use the **(+, −, −, −) metric signature**. The interval s² = (ct)² − x² − y² − z² is the fundamental invariant: it's the same in every inertial frame. The `FourVector` struct stores (ct, x, y, z) and delegates arithmetic to `nalgebra::Vector4`.

### Lorentz Invariance

Every boost matrix Λ satisfies **Λ^T η Λ = η**. The test suite verifies this numerically for multiple β values across all three axes. Velocity addition uses the relativistic formula β = (β₁ + β₂)/(1 + β₁β₂), which guarantees the result never exceeds 1.

### Schwarzschild Geometry

The Schwarzschild metric is encoded in its diagonal form:

- g_tt = −(1 − r_s/r)
- g_rr = (1 − r_s/r)⁻¹
- g_θθ = r²
- g_φφ = r²sin²θ

Christoffel symbols are computed analytically (all non-zero components are hardcoded) and stored in a 4×4×4 array. They satisfy Γ^μ_{νλ} = Γ^μ_{λν} (verified in tests).

### Geodesic Integration

Geodesics are integrated using **4th-order Runge–Kutta (RK4)** on the second-order ODE:

d²x^μ/dτ² = −Γ^μ_{νλ} (dx^ν/dτ)(dx^λ/dτ)

This is decomposed into a first-order system: (position, velocity) → (velocity, acceleration).

### Cosmology

The Friedmann equations are evaluated directly from density parameters. The ΛCDM model uses Planck 2018 best-fit values (H₀ ≈ 67.4 km/s/Mpc, Ω_m = 0.315, Ω_Λ = 0.685). Age is estimated via the integral formula t₀ ≈ 2/(3H₀√Ω_Λ) · asinh(√(Ω_Λ/Ω_m)).

---

## The Math

### Lorentz Factor

γ = 1/√(1 − β²), where β = v/c. At β = 0, γ = 1; as β → 1, γ → ∞.

### Spacetime Interval Classification

| s² | Type | Meaning |
|---|---|---|
| s² > 0 | Timelike | Causal connection possible, proper time τ = √(s²)/c |
| s² = 0 | Lightlike / Null | Connected only by light |
| s² < 0 | Spacelike | No causal connection, proper distance d = √(−s²) |

### Energy-Momentum Relation

**E² = (pc)² + (mc²)²**

At rest (p = 0): E = mc². For photons (m = 0): E = pc.

### Schwarzschild Radius

**r_s = 2GM/c²**

For the Sun: r_s ≈ 2953 m. For Earth: r_s ≈ 8.87 mm.

### Effective Potential

V_eff(r) = −GM/r + L²/(2r²) − GML²/(c²r³)

The last term is the **relativistic correction** that creates the ISCO at 3r_s and the photon sphere at 1.5r_s.

### Friedmann Equations

First: **H² = (8πG/3)ρ − k/a² + Λ/3**

Second: **ä/a = −(4πG/3)(ρ + 3p/c²) + Λ/3**

Matter-only (Λ = 0, ρ > 0): universe decelerates. With Λ > 0 and low ρ: universe accelerates.

### Shapiro Delay

Light passing near a mass M with impact parameter b experiences an extra time delay:

**Δt ≈ (4GM/c³) · ln(4r₁r₂/b²)**

This was one of the classic tests of general relativity (Shapiro 1964).

---

## License

MIT
