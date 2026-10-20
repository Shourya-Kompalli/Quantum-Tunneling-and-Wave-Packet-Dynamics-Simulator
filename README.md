# Quantum Tunneling and Wave-Packet Dynamics Simulator

A reproducible, interactive simulator for one-dimensional quantum scattering.
It propagates a normalized Gaussian wave packet through a rectangular
potential barrier, measures reflection and transmission, and compares the
computed results with the stationary-state barrier solution.

This is a numerical experiment—not an animation of a particle passing through
a wall. Probability densities, transmission plots, and sample CSV results are
computed from the time-dependent Schrödinger equation.

**Made With Vibe coding - Github Copilot**

## Features

- Fourier split-operator propagation with a documented second-order Strang
  splitting.
- Adjustable particle mass, incident energy, initial packet width and
  position, barrier height and width, grid spacing, time step, and run time.
- A complex absorbing potential at the spatial boundaries, with absorbed
  probability tracked separately from transmission.
- Playback animation, position-time density map, measured width and energy
  sweeps, and a probability-balance plot.
- CSV download of experiment parameters, measured probabilities, numerical
  resolution, boundary absorption, and stationary benchmarks.
- Automated checks for free-packet motion, conservation, analytic barrier
  limits, and numerical refinement.

## Quick start

Requires Python 3.10 or newer. In PowerShell, from this project directory:

```powershell
py -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

Open the local URL printed by Streamlit (usually `http://localhost:8501`).
Choose the experiment parameters and select **Run experiment**. The initial
packet is directed toward the barrier from whichever side it starts on. Enable
or disable the independent transmission sweeps in the sidebar.

For tests and data regeneration:

```powershell
python -m pip install -r requirements-dev.txt
python -m pytest
python generate_sample_data.py
```

See [HOW_TO_RUN.md](HOW_TO_RUN.md) for detailed setup and interface
instructions, and [REPORT.md](REPORT.md) for the numerical method, assumptions,
limitations, validation cases, and measured refinement results.

## Results included in the repository

[`sample_data/sample_results.csv`](sample_data/sample_results.csv) contains six
computed experiments: three barrier widths and three incident energies. Each
row includes measured reflection and transmission, packet-averaged stationary
benchmark, absorbed probability, conservation residual, and simulation
parameters. Regenerate the file with `python generate_sample_data.py`.

In the documented refinement case, measured transmission was `0.944900` at
1024 grid points and `Δt = 0.04`, and `0.927324` at 2048 grid points and
`Δt = 0.02`. The packet-averaged stationary benchmark was `0.946265`. These
are computed values for that particular finite packet and setup, not
universal probabilities or a claim of exact agreement.

## Project layout

```text
.
├── app.py                       # Interactive Streamlit interface
├── simulator.py                 # Split-operator solver and analytic benchmark
├── generate_sample_data.py      # Reproduce the included sample CSV
├── sample_data/
│   └── sample_results.csv       # Computed sample experiments
├── tests/
│   └── test_simulator.py        # Numerical and scientific validation
├── HOW_TO_RUN.md                # Detailed installation and usage
├── REPORT.md                    # Method, assumptions, limitations, evidence
├── requirements.txt             # Runtime dependencies
└── requirements-dev.txt         # Runtime and test dependencies
```

## Units and scientific interpretation

All quantities are dimensionless, with `ℏ = 1`. The packet's central momentum
is `p₀ = ±√(2mE₀)`, with its sign chosen to point toward the barrier; its
momentum standard deviation is `1/(2σₓ)`. Consequently, the packet has a finite
energy spread and its measured transmission need not equal the stationary
transmission at the central energy.

The stationary comparison averages the exact rectangular-barrier
transmission over the packet's initial momentum distribution. Time-dependent
transmission is instead measured from the propagated probability density.
Results are meaningful only when the outgoing packets have separated from the
barrier and boundary absorption is small. The interface warns when the
propagation time or absorbed probability suggests caution. Consult
[REPORT.md](REPORT.md) before interpreting results.

## License

No license has been specified. All rights remain with the author unless a
license is added.
