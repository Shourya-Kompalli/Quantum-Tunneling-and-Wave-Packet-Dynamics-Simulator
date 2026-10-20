# Numerical method, assumptions, and validation

## Model and units

The simulator solves the one-dimensional, non-relativistic time-dependent
Schrödinger equation

```text
i ℏ ∂ψ/∂t = [-ℏ²/(2m) ∂²/∂x² + V(x)] ψ
```

in dimensionless units with `ℏ = 1`. Length, mass, energy, and time are
dimensionless; the Schrödinger equation supplies their consistent relation.
The packet starts as

```text
ψ(x,0) ∝ exp[-(x - x₀)²/(4σₓ²) + i p₀x],
p₀ = √(2mE₀),  σₚ = 1/(2σₓ).
```

It is normalized by the discrete grid quadrature before propagation. The
barrier is a rectangular potential of the requested height and width, centered
at the configured barrier position. The packet's momentum points toward the
barrier: positive for an initial center to its left and negative for a center
to its right. `E₀` is the central kinetic energy, not the packet's exact mean
energy; the latter includes the finite momentum spread.

## Propagation and boundary treatment

The spatial grid is uniform and periodic for the Fourier representation.
Propagation uses the second-order Strang split-operator method:

```text
ψ(t + Δt) ≈ exp[-iVΔt/(2ℏ)] F⁻¹ exp[-iℏk²Δt/(2m)] F
             exp[-iVΔt/(2ℏ)] ψ(t).
```

The kinetic term is diagonal in Fourier space and the potential term is
diagonal in position space. The global time-discretization error of this
symmetric split is second order in `Δt` for sufficiently smooth evolution.

A quartic complex absorbing potential (CAP) is added only near both grid
edges. It damps outgoing packets before they reach the periodic boundary, which
reduces wraparound reflections. The real barrier is plotted separately from
the CAP. The probability removed by each CAP half-step is accumulated; the
reported balance residual is

```text
grid probability + cumulative CAP absorption - initial probability.
```

The grid probability is partitioned at each output time into left of the
barrier (reflected), right of the barrier (transmitted), and inside the
barrier. Final reflection and transmission are therefore spatially measured
quantities. They are reliable only after the scattered packets have separated
from the barrier. The interface estimates this time by the time for the packet
center to travel from its initial position to four spatial standard deviations
beyond the barrier. It warns if the primary run is shorter than this estimate
or if more than 1% of probability has reached the CAP. Sweep points use this
estimate individually because their group velocities differ. Probability
absorbed at the edges is reported separately, not counted as transmission.

## Stationary benchmark

For a rectangular barrier, the exact stationary transmission is used as a
benchmark. For `E < V₀`:

```text
T(E) = [1 + V₀² sinh²(κa) / (4E(V₀ - E))]⁻¹,
κ = √(2m(V₀ - E))/ℏ.
```

For `E > V₀`, replace `sinh(κa)` with `sin(qa)` and `κ` with
`q = √(2m(E - V₀))/ℏ`; the continuous limit is used at `E = V₀`. Zero height
or zero width gives unit transmission.

A localized packet is not a stationary state. Its displayed benchmark is
therefore the stationary transmission averaged over its Gaussian initial
momentum distribution, weighting each momentum by its probability density and
using `E = p²/(2m)`. The numerical time-dependent result and this benchmark
are distinct outputs: the former is measured after propagation; the latter is
an analytical comparison. Finite duration, energy spread, and any packet not
yet separated at measurement time can cause a difference.

## Validation and sample evidence

Run `python -m pytest` to execute the automated tests. They check:

1. Exact stationary-formula limits (no barrier and the `E = V₀` expression).
2. Free-packet motion against the expected group velocity `p₀/m` from either
   side of the barrier.
3. Probability balance with CAP loss accounted for.
4. Numerical transmission stability under simultaneous refinement from
   `N = 1024, Δt = 0.04` to `N = 2048, Δt = 0.02`, and comparison with the
   packet-averaged rectangular-barrier benchmark for a narrow-energy packet.
5. Rejection of invalid physical/numerical parameters.

The test's coarse and refined values are computed at test runtime; no expected
graph or transmission curve is embedded in the interface. The test thresholds
are tolerances, not fabricated measured outputs. `sample_data/sample_results.csv`
contains independent solver results and its corresponding analytical
benchmarks. Regenerate it with `python generate_sample_data.py`.

For the committed refinement case (`E₀ = 3.5`, `m = 1`, `σₓ = 8`,
`V₀ = 2`, `a = 1.5`, `x₀ = -25`, `x ∈ [-70, 70)`, `t = 23`), the observed
transmissions were:

| Grid | Δt | Numerical T | CAP absorption |
|---:|---:|---:|---:|
| 1024 | 0.04 | 0.944900 | 0.001942 |
| 2048 | 0.02 | 0.927324 | 0.001805 |
| Packet-averaged stationary benchmark | — | 0.946265 | — |

The absolute transmission change under this simultaneous spatial/time
refinement was about `0.0176`; the refined result differed from the
packet-averaged stationary benchmark by about `0.0189`. These are finite-grid,
finite-time computed values, not a claim of exact equality.

The generated sample also shows the computed transmission decreasing with
barrier width at `E₀ = 1.5` (`0.179742`, `0.005466`, `0.000714` for widths
`1`, `2.5`, `4`) and increasing across the selected energy samples at width
`2.5` (`0.000457`, `0.005466`, `0.873702` for energies `0.8`, `1.5`, `3.5`).
Each is paired with a separately computed packet-averaged stationary value in
the CSV. Results depend on the finite packet and setup; the values are not
universal transmission probabilities.

## Assumptions and limitations

- The model is one-dimensional and non-relativistic, with a single particle and
  no environmental decoherence, spin, or particle interactions.
- The Gaussian packet has a finite energy spread. Narrower energy distributions
  require larger `σₓ`, which in turn needs more spatial room and longer
  propagation to separate incident, reflected, and transmitted components.
- The spatial grid is finite and the Fourier kinetic operator is periodic.
  The CAP reduces, but cannot mathematically eliminate, reflection from the
  absorbing region. Check absorbed probability and repeat with a larger domain
  or weaker/earlier CAP if boundary effects are material.
- Grid spacing must resolve the shortest relevant de Broglie wavelength and
  barrier features; time steps must resolve both kinetic and potential phases.
  The automated refinement check covers one representative case, not every
  control setting.
- The stationary benchmark assumes an infinitely extended stationary
  scattering state incident on a single ideal rectangular barrier. It does not
  include finite-time measurement or the CAP.
- Transmission/reflection estimates depend on the chosen final time. Extend
  propagation until the components are spatially separated, while keeping
  absorption small enough that the transmitted/reflected estimates are not
  lost at the edges.
