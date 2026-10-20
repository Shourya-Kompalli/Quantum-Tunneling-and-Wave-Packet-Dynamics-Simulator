"""Interactive Streamlit front end for the quantum scattering simulator."""

from __future__ import annotations

import io

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from plotly.subplots import make_subplots

from simulator import (
    SimulationConfig,
    SimulationResult,
    packet_transmission_benchmark,
    recommended_measurement_time,
    run_simulation,
)


st.set_page_config(page_title="Quantum Tunneling Simulator", layout="wide")
st.title("Quantum Tunneling & Wave-Packet Dynamics")
st.caption(
    "Numerical results from the time-dependent Schrödinger equation. "
    "Dimensionless units are used with ℏ = 1."
)


def measurement_row(result: SimulationResult, run_kind: str) -> dict[str, float | str | int]:
    config = result.config
    dx = 2 * config.x_max / config.grid_points
    return {
        "run_kind": run_kind,
        "mass": config.mass,
        "incident_energy": config.incident_energy,
        "central_momentum": (
            1 if config.x0 < config.barrier_center else -1
        )
        * np.sqrt(2 * config.mass * config.incident_energy),
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
        "dx": dx,
        "dt": config.dt,
        "total_time": config.total_time,
        "absorber_start_fraction": config.absorber_start_fraction,
        "absorber_strength": config.absorber_strength,
    }


def animated_wavefunction(result: SimulationResult) -> go.Figure:
    fig = go.Figure(
        data=[
            go.Scatter(x=result.x, y=result.densities[0], name="|ψ(x,t)|²"),
            go.Scatter(
                x=result.x,
                y=result.potential,
                name="V(x)",
                yaxis="y2",
                line={"color": "#e45756", "dash": "dash"},
            ),
        ],
        layout={
            "title": "Computed probability density and barrier",
            "xaxis_title": "Position x (length units)",
            "yaxis": {"title": "Probability density"},
            "yaxis2": {
                "title": "Potential energy",
                "overlaying": "y",
                "side": "right",
                "showgrid": False,
            },
            "legend": {"orientation": "h"},
            "updatemenus": [
                {
                    "type": "buttons",
                    "showactive": False,
                    "buttons": [
                        {
                            "label": "Play",
                            "method": "animate",
                            "args": [
                                None,
                                {
                                    "frame": {"duration": 60, "redraw": True},
                                    "fromcurrent": True,
                                },
                            ],
                        },
                        {
                            "label": "Pause",
                            "method": "animate",
                            "args": [
                                [None],
                                {"frame": {"duration": 0, "redraw": False}, "mode": "immediate"},
                            ],
                        },
                    ],
                }
            ],
            "sliders": [
                {
                    "currentvalue": {"prefix": "Time (units): "},
                    "steps": [
                        {
                            "label": f"{time:.2f}",
                            "method": "animate",
                            "args": [
                                [str(index)],
                                {
                                    "mode": "immediate",
                                    "frame": {"duration": 0, "redraw": True},
                                    "transition": {"duration": 0},
                                },
                            ],
                        }
                        for index, time in enumerate(result.times)
                    ],
                }
            ],
        },
        frames=[
            go.Frame(
                data=[go.Scatter(y=density)],
                traces=[0],
                name=str(index),
            )
            for index, density in enumerate(result.densities)
        ],
    )
    return fig


def density_time_figure(result: SimulationResult) -> go.Figure:
    fig = go.Figure(
        go.Heatmap(
            x=result.x,
            y=result.times,
            z=result.densities,
            colorscale="Viridis",
            colorbar={"title": "|ψ|²"},
        )
    )
    fig.update_layout(
        title="Probability density versus position and time",
        xaxis_title="Position x (length units)",
        yaxis_title="Time (units)",
    )
    return fig


def sweep_result(
    base: SimulationConfig, *, energy: float | None = None, width: float | None = None
) -> SimulationResult:
    values = {
        "incident_energy": base.incident_energy if energy is None else energy,
        "barrier_width": base.barrier_width if width is None else width,
    }
    config = SimulationConfig(**{**base.__dict__, **values})
    config = SimulationConfig(
        **{
            **config.__dict__,
            "total_time": recommended_measurement_time(config),
        }
    )
    return run_simulation(config, snapshot_stride=max(1, int(np.ceil(config.total_time / config.dt))))


with st.sidebar:
    st.header("Experiment controls")
    mass = st.number_input("Particle mass", min_value=0.1, max_value=10.0, value=1.0, step=0.1)
    energy = st.number_input(
        "Incident energy", min_value=0.05, max_value=20.0, value=1.5, step=0.1
    )
    sigma = st.number_input(
        "Initial packet width σₓ", min_value=0.5, max_value=12.0, value=3.5, step=0.5
    )
    x0 = st.number_input(
        "Initial packet center x₀", min_value=-40.0, max_value=40.0, value=-25.0, step=1.0
    )
    height = st.number_input(
        "Barrier height", min_value=0.0, max_value=20.0, value=2.5, step=0.1
    )
    width = st.number_input(
        "Barrier width", min_value=0.0, max_value=20.0, value=3.0, step=0.25
    )
    grid_points = st.selectbox("Spatial grid points", [1024, 2048, 4096], index=1)
    dt = st.number_input(
        "Time step Δt", min_value=0.002, max_value=0.1, value=0.02, step=0.002, format="%.3f"
    )
    total_time = st.number_input(
        "Total propagation time", min_value=5.0, max_value=100.0, value=28.0, step=1.0
    )
    run_sweeps = st.checkbox("Run numerical width and energy sweeps", value=True)
    run_button = st.button("Run experiment", type="primary")

if run_button and np.isclose(x0, 0.0):
    st.error("The initial packet center must be on one side of the barrier, not at x = 0.")

if run_button and not np.isclose(x0, 0.0):
    config = SimulationConfig(
        mass=float(mass),
        incident_energy=float(energy),
        sigma=float(sigma),
        barrier_height=float(height),
        barrier_width=float(width),
        x0=float(x0),
        grid_points=int(grid_points),
        dt=float(dt),
        total_time=float(total_time),
    )
    config.validate()
    with st.spinner("Propagating the wave packet and measuring probabilities..."):
        main_result = run_simulation(config)
        results = [measurement_row(main_result, "primary")]

        width_results: list[dict[str, float | str | int]] = []
        energy_results: list[dict[str, float | str | int]] = []
        if run_sweeps:
            width_values = np.linspace(max(0.25, width * 0.25), max(1.0, width * 1.75), 6)
            energy_values = np.linspace(max(0.3, energy * 0.5), max(energy * 1.5, height * 1.4), 6)
            for test_width in width_values:
                result = sweep_result(config, width=float(test_width))
                row = measurement_row(result, "width_sweep")
                width_results.append(row)
                results.append(row)
            for test_energy in energy_values:
                result = sweep_result(config, energy=float(test_energy))
                row = measurement_row(result, "energy_sweep")
                energy_results.append(row)
                results.append(row)

    st.session_state["main_result"] = main_result
    st.session_state["measurement_rows"] = results
    st.session_state["width_results"] = width_results
    st.session_state["energy_results"] = energy_results

if "main_result" not in st.session_state:
    st.info("Choose physical parameters and click **Run experiment** to calculate results.")
else:
    result = st.session_state["main_result"]
    final_row = measurement_row(result, "primary")
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Measured transmission", f"{result.transmission_probability:.5f}")
    col2.metric("Measured reflection", f"{result.reflection_probability:.5f}")
    col3.metric("Probability absorbed at edges", f"{result.final_absorbed_probability:.3e}")
    col4.metric("Max balance error", f"{result.max_abs_conservation_error:.3e}")

    st.caption(
        f"Packet-averaged stationary benchmark: "
        f"{final_row['stationary_packet_benchmark']:.5f}. "
        "The benchmark averages the exact rectangular-barrier stationary solution "
        "over the packet's initial momentum distribution; it is not simulation data."
    )
    recommended_time = recommended_measurement_time(result.config)
    if abs(result.config.x0 - result.config.barrier_center) < 3 * result.config.sigma:
        st.warning(
            "The initial packet is within three spatial widths of the barrier; "
            "the starting state may already overlap the barrier."
        )
    if result.config.total_time < recommended_time:
        st.warning(
            f"The selected propagation time ({result.config.total_time:.2f}) is shorter "
            f"than the estimated packet-separation time ({recommended_time:.2f}). "
            "Extend the run before interpreting final reflection/transmission."
        )
    if result.final_absorbed_probability > 0.01:
        st.warning(
            "More than 1% of the probability reached the absorbing boundaries before "
            "measurement. Increase the domain or shorten the run and verify convergence."
        )
    tab_animation, tab_map, tab_sweeps, tab_conservation, tab_data = st.tabs(
        ["Animation", "Density vs time", "Transmission sweeps", "Probability balance", "CSV data"]
    )
    with tab_animation:
        st.plotly_chart(animated_wavefunction(result), width="stretch")
    with tab_map:
        st.plotly_chart(density_time_figure(result), width="stretch")
    with tab_sweeps:
        width_results = st.session_state["width_results"]
        energy_results = st.session_state["energy_results"]
        if width_results or energy_results:
            fig = make_subplots(
                rows=1,
                cols=2,
                subplot_titles=("Transmission vs barrier width", "Transmission vs incident energy"),
            )
            if width_results:
                frame = pd.DataFrame(width_results).sort_values("barrier_width")
                fig.add_trace(
                    go.Scatter(
                        x=frame["barrier_width"],
                        y=frame["transmission_probability"],
                        mode="lines+markers",
                        name="Numerical wave-packet measurements",
                    ),
                    row=1,
                    col=1,
                )
                fig.add_trace(
                    go.Scatter(
                        x=frame["barrier_width"],
                        y=frame["stationary_packet_benchmark"],
                        mode="lines",
                        name="Stationary packet benchmark",
                        line={"dash": "dash"},
                    ),
                    row=1,
                    col=1,
                )
            if energy_results:
                frame = pd.DataFrame(energy_results).sort_values("incident_energy")
                fig.add_trace(
                    go.Scatter(
                        x=frame["incident_energy"],
                        y=frame["transmission_probability"],
                        mode="lines+markers",
                        name="Numerical wave-packet measurements",
                        showlegend=False,
                    ),
                    row=1,
                    col=2,
                )
                fig.add_trace(
                    go.Scatter(
                        x=frame["incident_energy"],
                        y=frame["stationary_packet_benchmark"],
                        mode="lines",
                        name="Stationary packet benchmark",
                        line={"dash": "dash"},
                        showlegend=False,
                    ),
                    row=1,
                    col=2,
                )
            fig.update_xaxes(title_text="Barrier width (length units)", row=1, col=1)
            fig.update_xaxes(title_text="Incident energy (energy units)", row=1, col=2)
            fig.update_yaxes(title_text="Transmission probability", range=[0, 1], row=1, col=1)
            fig.update_yaxes(title_text="Transmission probability", range=[0, 1], row=1, col=2)
            st.plotly_chart(fig, width="stretch")
            st.caption(
                "Solid markers are transmission measured from independent time-dependent "
                "simulations. Dashed lines are the stationary rectangular-barrier benchmark."
            )
        else:
            st.info("Enable numerical sweeps and rerun the experiment to populate these graphs.")
    with tab_conservation:
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=result.times, y=result.norm, name="Probability remaining in grid"))
        fig.add_trace(go.Scatter(x=result.times, y=result.absorbed, name="Cumulative CAP absorption"))
        fig.add_trace(
            go.Scatter(
                x=result.times,
                y=result.conservation_error,
                name="Balance residual (remaining + absorbed − 1)",
                yaxis="y2",
            )
        )
        fig.update_layout(
            title="Probability accounting during propagation",
            xaxis_title="Time (units)",
            yaxis={"title": "Probability"},
            yaxis2={"title": "Balance residual", "overlaying": "y", "side": "right"},
        )
        st.plotly_chart(fig, width="stretch")
        st.caption(
            "The CAP removes outgoing probability before it reaches the periodic FFT boundary. "
            "The residual compares remaining grid norm plus accumulated CAP loss with the "
            "initial normalized probability."
        )
    with tab_data:
        frame = pd.DataFrame(st.session_state["measurement_rows"])
        st.dataframe(frame, width="stretch")
        csv_buffer = io.StringIO()
        frame.to_csv(csv_buffer, index=False)
        st.download_button(
            "Download experimental results (CSV)",
            csv_buffer.getvalue(),
            file_name="quantum_tunneling_results.csv",
            mime="text/csv",
        )
