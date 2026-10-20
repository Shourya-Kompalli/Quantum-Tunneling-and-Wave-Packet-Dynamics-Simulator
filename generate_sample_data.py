"""Regenerate the small, computed sample data set."""

from __future__ import annotations

import csv
from pathlib import Path

from simulator import (
    SimulationConfig,
    packet_transmission_benchmark,
    recommended_measurement_time,
    run_simulation,
)


def main() -> None:
    base = SimulationConfig(
        mass=1.0,
        incident_energy=1.5,
        sigma=4.0,
        barrier_height=2.5,
        barrier_width=2.5,
        x0=-25.0,
        x_max=60.0,
        grid_points=1024,
        dt=0.03,
        total_time=30.0,
    )
    configurations = [
        ("width_sweep", SimulationConfig(**{**base.__dict__, "barrier_width": width}))
        for width in (1.0, 2.5, 4.0)
    ] + [
        ("energy_sweep", SimulationConfig(**{**base.__dict__, "incident_energy": energy}))
        for energy in (0.8, 1.5, 3.5)
    ]

    rows: list[dict[str, float | str | int]] = []
    for run_kind, config in configurations:
        config = SimulationConfig(
            **{
                **config.__dict__,
                "total_time": recommended_measurement_time(config),
            }
        )
        result = run_simulation(config, snapshot_stride=100)
        rows.append(
            {
                "run_kind": run_kind,
                "mass": config.mass,
                "incident_energy": config.incident_energy,
                "central_momentum": (
                    1 if config.x0 < config.barrier_center else -1
                )
                * (2 * config.mass * config.incident_energy) ** 0.5,
                "momentum_sigma": 1 / (2 * config.sigma),
                "packet_sigma_x": config.sigma,
                "initial_position": config.x0,
                "barrier_height": config.barrier_height,
                "barrier_width": config.barrier_width,
                "barrier_center": config.barrier_center,
                "reflection_probability": result.reflection_probability,
                "transmission_probability": result.transmission_probability,
                "stationary_packet_benchmark": packet_transmission_benchmark(
                    config.incident_energy,
                    config.mass,
                    config.sigma,
                    config.barrier_height,
                    config.barrier_width,
                ),
                "absorbed_probability": result.final_absorbed_probability,
                "final_norm": result.final_norm,
                "max_abs_balance_error": result.max_abs_conservation_error,
                "grid_points": config.grid_points,
                "x_max": config.x_max,
                "dx": 2 * config.x_max / config.grid_points,
                "dt": config.dt,
                "total_time": config.total_time,
                "absorber_start_fraction": config.absorber_start_fraction,
                "absorber_strength": config.absorber_strength,
            }
        )

    output = Path(__file__).parent / "sample_data" / "sample_results.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} computed experiments to {output}")


if __name__ == "__main__":
    main()
