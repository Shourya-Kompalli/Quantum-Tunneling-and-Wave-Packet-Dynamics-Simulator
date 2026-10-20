"""One-dimensional split-operator wave-packet scattering in dimensionless units."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

FloatArray = NDArray[np.float64]
ComplexArray = NDArray[np.complex128]

HBAR = 1.0


@dataclass(frozen=True)
class SimulationConfig:
    mass: float = 1.0
    incident_energy: float = 1.5
    sigma: float = 3.5
    barrier_height: float = 2.5
    barrier_width: float = 3.0
    barrier_center: float = 0.0
    x0: float = -25.0
    x_max: float = 60.0
    grid_points: int = 2048
    dt: float = 0.02
    total_time: float = 28.0
    absorber_start_fraction: float = 0.8
    absorber_strength: float = 1.5

    def validate(self) -> None:
        if self.mass <= 0:
            raise ValueError("mass must be positive")
        if self.incident_energy <= 0:
            raise ValueError("incident_energy must be positive")
        if self.sigma <= 0:
            raise ValueError("sigma must be positive")
        if self.barrier_height < 0 or self.barrier_width < 0:
            raise ValueError("barrier height and width cannot be negative")
        if self.x_max <= 0 or not -self.x_max < self.x0 < self.x_max:
            raise ValueError("x0 must lie inside the spatial domain")
        if np.isclose(self.x0, self.barrier_center):
            raise ValueError("x0 must be on one side of the barrier")
        if self.grid_points < 64 or self.grid_points % 2:
            raise ValueError("grid_points must be an even integer of at least 64")
        if self.dt <= 0 or self.total_time <= 0:
            raise ValueError("dt and total_time must be positive")
        if not 0.5 < self.absorber_start_fraction < 1:
            raise ValueError("absorber_start_fraction must be between 0.5 and 1")
        if self.absorber_strength < 0:
            raise ValueError("absorber_strength cannot be negative")
        if abs(self.barrier_center) + self.barrier_width / 2 >= self.x_max:
            raise ValueError("the barrier must fit inside the spatial domain")


@dataclass
class SimulationResult:
    x: FloatArray
    times: FloatArray
    densities: FloatArray
    potential: FloatArray
    reflected: FloatArray
    transmitted: FloatArray
    inside_barrier: FloatArray
    norm: FloatArray
    absorbed: FloatArray
    conservation_error: FloatArray
    final_wavefunction: ComplexArray
    config: SimulationConfig

    @property
    def transmission_probability(self) -> float:
        return float(self.transmitted[-1])

    @property
    def reflection_probability(self) -> float:
        return float(self.reflected[-1])

    @property
    def final_norm(self) -> float:
        return float(self.norm[-1])

    @property
    def final_absorbed_probability(self) -> float:
        return float(self.absorbed[-1])

    @property
    def max_abs_conservation_error(self) -> float:
        return float(np.max(np.abs(self.conservation_error)))


def _normalized_gaussian(
    x: FloatArray, x0: float, sigma: float, momentum: float
) -> ComplexArray:
    psi = np.exp(-((x - x0) ** 2) / (4 * sigma**2) + 1j * momentum * x)
    psi = psi.astype(np.complex128)
    return psi / np.sqrt(np.sum(np.abs(psi) ** 2) * (x[1] - x[0]))


def _probability_parts(
    density: FloatArray, x: FloatArray, config: SimulationConfig, dx: float
) -> tuple[float, float, float, float]:
    left = config.barrier_center - config.barrier_width / 2
    right = config.barrier_center + config.barrier_width / 2
    left_probability = float(np.sum(density[x < left]) * dx)
    right_probability = float(np.sum(density[x > right]) * dx)
    if config.x0 < config.barrier_center:
        reflected, transmitted = left_probability, right_probability
    else:
        reflected, transmitted = right_probability, left_probability
    inside = float(np.sum(density[(x >= left) & (x <= right)]) * dx)
    norm = float(np.sum(density) * dx)
    return reflected, transmitted, inside, norm


def run_simulation(
    config: SimulationConfig, snapshot_stride: int = 10
) -> SimulationResult:
    """Propagate a normalized packet using a Strang-split Fourier operator.

    The complex absorbing potential (CAP) damps outgoing waves near the domain
    edges. Its norm loss is tracked separately so probability balance remains
    measurable rather than treating CAP loss as physical transmission.
    """
    config.validate()
    if snapshot_stride < 1:
        raise ValueError("snapshot_stride must be at least 1")

    x = np.linspace(
        -config.x_max,
        config.x_max,
        config.grid_points,
        endpoint=False,
        dtype=np.float64,
    )
    dx = float(x[1] - x[0])
    momentum = np.sqrt(2 * config.mass * config.incident_energy)
    if config.x0 > config.barrier_center:
        momentum = -momentum
    psi = _normalized_gaussian(x, config.x0, config.sigma, momentum)

    left = config.barrier_center - config.barrier_width / 2
    right = config.barrier_center + config.barrier_width / 2
    barrier = (x >= left) & (x <= right) if config.barrier_width else np.zeros_like(x, dtype=bool)
    real_potential = np.where(barrier, config.barrier_height, 0.0)
    absorber_start = config.absorber_start_fraction * config.x_max
    absorber_coordinate = np.maximum(np.abs(x) - absorber_start, 0.0) / (
        config.x_max - absorber_start
    )
    cap = config.absorber_strength * absorber_coordinate**4
    potential = real_potential.astype(np.complex128) - 1j * cap

    wave_numbers = 2 * np.pi * np.fft.fftfreq(config.grid_points, d=dx)
    kinetic_energy = HBAR**2 * wave_numbers**2 / (2 * config.mass)
    initial_norm = float(np.sum(np.abs(psi) ** 2) * dx)
    absorbed_probability = 0.0

    times: list[float] = []
    densities: list[FloatArray] = []
    reflected: list[float] = []
    transmitted: list[float] = []
    inside: list[float] = []
    norms: list[float] = []
    absorbed: list[float] = []
    conservation_errors: list[float] = []

    def record(time: float) -> None:
        density = np.abs(psi) ** 2
        r, t, b, norm = _probability_parts(density, x, config, dx)
        times.append(time)
        densities.append(density.astype(np.float32))
        reflected.append(r)
        transmitted.append(t)
        inside.append(b)
        norms.append(norm)
        absorbed.append(absorbed_probability)
        conservation_errors.append(norm + absorbed_probability - initial_norm)

    record(0.0)
    elapsed = 0.0
    step = 0
    while elapsed < config.total_time:
        step_dt = min(config.dt, config.total_time - elapsed)
        potential_half_step = np.exp(-1j * potential * step_dt / (2 * HBAR))

        norm_before = float(np.sum(np.abs(psi) ** 2) * dx)
        psi *= potential_half_step
        norm_after = float(np.sum(np.abs(psi) ** 2) * dx)
        absorbed_probability += norm_before - norm_after

        psi = np.fft.ifft(
            np.fft.fft(psi)
            * np.exp(-1j * kinetic_energy * step_dt / HBAR)
        )

        norm_before = float(np.sum(np.abs(psi) ** 2) * dx)
        psi *= potential_half_step
        norm_after = float(np.sum(np.abs(psi) ** 2) * dx)
        absorbed_probability += norm_before - norm_after

        elapsed += step_dt
        step += 1
        if step % snapshot_stride == 0 or elapsed >= config.total_time:
            record(elapsed)

    return SimulationResult(
        x=x,
        times=np.asarray(times),
        densities=np.asarray(densities),
        potential=real_potential.astype(np.float64),
        reflected=np.asarray(reflected),
        transmitted=np.asarray(transmitted),
        inside_barrier=np.asarray(inside),
        norm=np.asarray(norms),
        absorbed=np.asarray(absorbed),
        conservation_error=np.asarray(conservation_errors),
        final_wavefunction=psi,
        config=config,
    )


def rectangular_barrier_transmission(
    energy: float, mass: float, height: float, width: float
) -> float:
    """Stationary transmission for a rectangular barrier with hbar = 1."""
    if energy <= 0 or mass <= 0 or height < 0 or width < 0:
        raise ValueError("energy and mass must be positive; barrier values nonnegative")
    if height == 0 or width == 0:
        return 1.0
    if np.isclose(energy, height, rtol=1e-10, atol=1e-12):
        correction = mass * height * width**2 / (2 * HBAR**2)
        return float(1 / (1 + correction))
    if energy < height:
        kappa = np.sqrt(2 * mass * (height - energy)) / HBAR
        exponent = 2 * kappa * width
        if exponent > 700:
            return 0.0
        numerator = height**2 * np.sinh(kappa * width) ** 2
        denominator = 4 * energy * (height - energy)
    else:
        wave_number = np.sqrt(2 * mass * (energy - height)) / HBAR
        numerator = height**2 * np.sin(wave_number * width) ** 2
        denominator = 4 * energy * (energy - height)
    return float(1 / (1 + numerator / denominator))


def packet_transmission_benchmark(
    energy: float,
    mass: float,
    sigma: float,
    height: float,
    width: float,
    samples: int = 2001,
) -> float:
    """Average the stationary solution over the packet's Gaussian momentum spread."""
    if energy <= 0 or mass <= 0 or sigma <= 0 or samples < 101:
        raise ValueError("invalid packet benchmark parameters")
    central_momentum = np.sqrt(2 * mass * energy)
    momentum_std = HBAR / (2 * sigma)
    momenta = np.linspace(
        central_momentum - 8 * momentum_std,
        central_momentum + 8 * momentum_std,
        samples,
    )
    weights = np.exp(-0.5 * ((momenta - central_momentum) / momentum_std) ** 2)
    stationary = np.array(
        [
            rectangular_barrier_transmission(
                max(float(p * p / (2 * mass)), np.finfo(float).tiny),
                mass,
                height,
                width,
            )
            for p in momenta
        ]
    )
    return float(np.trapezoid(weights * stationary, momenta) / np.trapezoid(weights, momenta))


def recommended_measurement_time(config: SimulationConfig) -> float:
    """Estimate when both outgoing packet centers are four σx from the barrier."""
    config.validate()
    speed = np.sqrt(2 * config.incident_energy / config.mass)
    travel_distance = (
        abs(config.x0 - config.barrier_center)
        + config.barrier_width / 2
        + 4 * config.sigma
    )
    return float(travel_distance / speed)
