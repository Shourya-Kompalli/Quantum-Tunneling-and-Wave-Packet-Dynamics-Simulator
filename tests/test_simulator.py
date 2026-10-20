import numpy as np

from simulator import (
    SimulationConfig,
    packet_transmission_benchmark,
    rectangular_barrier_transmission,
    run_simulation,
)


def test_rectangular_barrier_limits_and_threshold() -> None:
    assert rectangular_barrier_transmission(1.0, 1.0, 0.0, 2.0) == 1.0
    assert rectangular_barrier_transmission(1.0, 1.0, 3.0, 0.0) == 1.0
    assert 0 < rectangular_barrier_transmission(1.0, 1.0, 2.0, 1.0) < 1
    threshold_transmission = rectangular_barrier_transmission(2.0, 1.0, 2.0, 1.0)
    assert np.isclose(threshold_transmission, 0.5)


def test_free_packet_moves_at_expected_group_velocity() -> None:
    config = SimulationConfig(
        incident_energy=1.125,
        sigma=3.0,
        barrier_height=0.0,
        barrier_width=0.0,
        x0=-10.0,
        x_max=40.0,
        grid_points=512,
        dt=0.02,
        total_time=2.0,
    )
    result = run_simulation(config, snapshot_stride=10)
    dx = result.x[1] - result.x[0]
    initial_center = np.sum(result.x * result.densities[0]) * dx
    final_center = np.sum(result.x * result.densities[-1]) * dx
    expected_center = -10.0 + np.sqrt(2 * config.incident_energy / config.mass) * 2.0
    assert abs(initial_center + 10.0) < 1e-3
    assert abs(final_center - expected_center) < 0.08
    assert result.max_abs_conservation_error < 1e-10


def test_packet_started_right_of_barrier_moves_toward_it() -> None:
    config = SimulationConfig(
        incident_energy=1.125,
        sigma=3.0,
        barrier_height=0.0,
        barrier_width=0.0,
        x0=10.0,
        x_max=40.0,
        grid_points=512,
        dt=0.02,
        total_time=2.0,
    )
    result = run_simulation(config, snapshot_stride=100)
    dx = result.x[1] - result.x[0]
    final_center = np.sum(result.x * result.densities[-1]) * dx
    expected_center = 10.0 - np.sqrt(2 * config.incident_energy / config.mass) * 2.0
    assert abs(final_center - expected_center) < 0.08
    assert result.reflection_probability > result.transmission_probability


def test_scattering_matches_packet_averaged_stationary_solution() -> None:
    common = {
        "mass": 1.0,
        "incident_energy": 3.5,
        "sigma": 8.0,
        "barrier_height": 2.0,
        "barrier_width": 1.5,
        "x0": -25.0,
        "x_max": 70.0,
        "total_time": 23.0,
    }
    coarse = run_simulation(
        SimulationConfig(**common, grid_points=1024, dt=0.04),
        snapshot_stride=100,
    )
    refined = run_simulation(
        SimulationConfig(**common, grid_points=2048, dt=0.02),
        snapshot_stride=100,
    )
    benchmark = packet_transmission_benchmark(3.5, 1.0, 8.0, 2.0, 1.5)

    assert abs(coarse.transmission_probability - refined.transmission_probability) < 0.03
    assert abs(refined.transmission_probability - benchmark) < 0.04
    assert refined.final_absorbed_probability < 0.01
    assert refined.max_abs_conservation_error < 1e-10
    assert abs(
        refined.reflection_probability
        + refined.transmission_probability
        + refined.inside_barrier[-1]
        + refined.final_absorbed_probability
        - 1
    ) < 1e-10


def test_invalid_configuration_is_rejected() -> None:
    try:
        run_simulation(SimulationConfig(mass=0))
    except ValueError as error:
        assert "mass" in str(error)
    else:
        raise AssertionError("invalid mass should be rejected")
