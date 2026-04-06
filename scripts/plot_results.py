from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
FIGURES = ROOT / "report" / "figures"


def setup_style() -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 10,
            "legend.fontsize": 9,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "figure.dpi": 200,
            "savefig.dpi": 300,
        }
    )


def savefig(fig: plt.Figure, name: str) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(FIGURES / name, bbox_inches="tight")
    plt.close(fig)


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"missing required file: {path}")
    return pd.read_csv(path)


def plot_free_fall() -> None:
    df = read_csv(RESULTS / "verification" / "free_fall" / "diagnostics.csv")
    g = 9.81
    z0 = 2.0
    t = df["time"].to_numpy()
    z_num = df["center_of_mass_z"].to_numpy()
    z_exact = z0 - 0.5 * g * t**2
    err = np.abs(z_num - z_exact)
    vz_mag_num = df["max_speed"].to_numpy()
    vz_mag_exact = g * t

    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    ax.plot(t, z_num, label="Numerical", linewidth=2.2, color="#0b5fa5")
    ax.plot(t, z_exact, "--", label="Analytical", linewidth=2.0, color="#d95f02")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Particle height z (m)")
    ax.set_title("Free-fall verification")
    ax.legend(frameon=True)
    savefig(fig, "free_fall_numerical_vs_analytical.png")

    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    ax.plot(t, err, color="#b22222", linewidth=1.8)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Absolute error (m)")
    ax.set_title("Free-fall trajectory error")
    savefig(fig, "free_fall_error_vs_time.png")

    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    ax.plot(t, vz_mag_num, label="Numerical", linewidth=2.2, color="#0b5fa5")
    ax.plot(t, vz_mag_exact, "--", label="Analytical", linewidth=2.0, color="#d95f02")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel(r"$|v_z|$ (m/s)")
    ax.set_title("Free-fall velocity verification")
    ax.legend(frameon=True)
    savefig(fig, "free_fall_velocity_numerical_vs_analytical.png")


def plot_timestep_sensitivity() -> None:
    df = read_csv(RESULTS / "verification" / "timestep_sensitivity" / "timestep_scan.csv")
    dt_ms = 1.0e3 * df["dt"].to_numpy()
    err_mm = 1.0e3 * df["final_z_error"].to_numpy()

    fig, ax = plt.subplots(figsize=(5.8, 3.8))
    ax.plot(dt_ms, err_mm, marker="o", linewidth=2.0, color="#b22222")
    ax.set_xticks(dt_ms)
    ax.set_xlabel(r"Timestep $\Delta t$ (ms)")
    ax.set_ylabel("Final height error (mm)")
    ax.set_title("Free-fall final height error versus timestep")
    savefig(fig, "timestep_sensitivity.png")


def plot_constant_velocity() -> None:
    df = read_csv(RESULTS / "verification" / "constant_velocity" / "diagnostics.csv")
    t = df["time"].to_numpy()
    x_num = df["center_of_mass_x"].to_numpy()
    y_num = df["center_of_mass_y"].to_numpy()
    z_num = df["center_of_mass_z"].to_numpy()

    x0, y0, z0 = 0.5, 0.6, 0.7
    vx, vy, vz = 0.8, -0.35, 0.5
    x_exact = x0 + vx * t
    y_exact = y0 + vy * t
    z_exact = z0 + vz * t
    err = np.sqrt((x_num - x_exact) ** 2 + (y_num - y_exact) ** 2 + (z_num - z_exact) ** 2)

    fig, axes = plt.subplots(3, 1, figsize=(6.2, 6.8), sharex=True)
    components = [
        ("x", x_num, x_exact, "#0b5fa5"),
        ("y", y_num, y_exact, "#1b9e77"),
        ("z", z_num, z_exact, "#d95f02"),
    ]
    for ax, (label, num, exact, color) in zip(axes, components):
        ax.plot(t, num, linewidth=2.0, color=color, label=f"{label} numerical")
        ax.plot(t, exact, "--", linewidth=1.8, color="black", label=f"{label} analytical")
        ax.set_ylabel(f"{label}(t) (m)")
        ax.legend(frameon=True, loc="best")
    axes[-1].set_xlabel("Time (s)")
    fig.suptitle("Constant-velocity verification: position comparison", y=1.02)
    savefig(fig, "constant_velocity_numerical_vs_analytical.png")

    snapshot_dir = RESULTS / "verification" / "constant_velocity"
    snapshots = sorted(snapshot_dir.glob("particles_step_*.csv"))
    if snapshots:
        rows = [read_csv(path).iloc[0] for path in snapshots]
        times = []
        vx_num = []
        vy_num = []
        vz_num = []
        for path, row in zip(snapshots, rows):
            step = int(path.stem.split("_")[-1])
            times.append(step * 2.0e-4)
            vx_num.append(row["vx"])
            vy_num.append(row["vy"])
            vz_num.append(row["vz"])

        fig, axes = plt.subplots(3, 1, figsize=(6.2, 6.8), sharex=True)
        v_components = [
            ("v_x", np.array(vx_num), vx, "#0b5fa5"),
            ("v_y", np.array(vy_num), vy, "#1b9e77"),
            ("v_z", np.array(vz_num), vz, "#d95f02"),
        ]
        for ax, (label, num, exact, color) in zip(axes, v_components):
            ax.plot(times, num, marker="o", linewidth=1.8, color=color, label=f"{label} numerical")
            ax.axhline(exact, linestyle="--", linewidth=1.8, color="black", label=f"{label} analytical")
            ax.set_ylabel(f"{label} (m/s)")
            ax.legend(frameon=True, loc="best")
        axes[-1].set_xlabel("Time (s)")
        fig.suptitle("Constant-velocity verification: velocity comparison", y=1.02)
        savefig(fig, "constant_velocity_velocity_check.png")


def _find_bounce_peaks(z: np.ndarray) -> np.ndarray:
    peaks: list[float] = []
    for i in range(1, len(z) - 1):
        if z[i] > z[i - 1] and z[i] >= z[i + 1]:
            peaks.append(z[i])
    return np.array(peaks)


def plot_bounce_and_damping() -> None:
    light = read_csv(RESULTS / "verification" / "bounce_light_damping" / "diagnostics.csv")
    heavy = read_csv(RESULTS / "verification" / "bounce_heavy_damping" / "diagnostics.csv")

    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    ax.plot(light["time"], light["center_of_mass_z"], linewidth=1.8, color="#1f4e79", label="Light damping")
    ax.plot(heavy["time"], heavy["center_of_mass_z"], linewidth=1.8, color="#b22222", label="Heavy damping")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Particle height z (m)")
    ax.set_title("Single-particle bounce comparison")
    ax.legend(frameon=True)
    savefig(fig, "bounce_height_vs_time.png")

    fig, axes = plt.subplots(2, 1, figsize=(6.2, 5.4), sharex=True)
    max_ke = max(light["kinetic_energy"].max(), heavy["kinetic_energy"].max())
    ke_series = [
        (axes[0], light, "Light damping", "#1f4e79"),
        (axes[1], heavy, "Heavy damping", "#b22222"),
    ]
    for ax, data, label, color in ke_series:
        time = data["time"].to_numpy()
        ke = data["kinetic_energy"].to_numpy()
        stride = max(1, len(time) // 140)
        ax.plot(time, ke, linewidth=1.8, color=color)
        ax.plot(
            time[::stride],
            ke[::stride],
            linestyle="none",
            marker="o",
            markersize=2.8,
            color=color,
            alpha=0.75,
        )
        ax.set_ylabel("KE (J)")
        ax.set_ylim(0.0, 1.05 * max_ke)
        ax.set_title(label)
        ax.grid(True, alpha=0.35)
    axes[-1].set_xlabel("Time (s)")
    fig.suptitle("Bounce verification: kinetic energy history", y=1.02)
    savefig(fig, "bounce_kinetic_energy_vs_time.png")

    radius = 0.04
    light_overlap = np.maximum(0.0, radius - light["center_of_mass_z"].to_numpy())
    heavy_overlap = np.maximum(0.0, radius - heavy["center_of_mass_z"].to_numpy())
    light_contact = light_overlap > 0.0
    heavy_contact = heavy_overlap > 0.0

    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    ax.scatter(
        light.loc[light_contact, "time"],
        light_overlap[light_contact],
        s=14,
        color="#0b5fa5",
        label="Light damping",
    )
    ax.scatter(
        heavy.loc[heavy_contact, "time"],
        heavy_overlap[heavy_contact],
        s=14,
        color="#d95f02",
        marker="s",
        label="Heavy damping",
    )
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Wall overlap (m)")
    ax.set_title("Bounce verification: wall-contact overlap")
    ax.legend(frameon=True)
    savefig(fig, "bounce_wall_overlap.png")

    scan = read_csv(RESULTS / "science_bonus" / "damping_scan.csv")
    cloud = read_csv(RESULTS / "science_bonus" / "cloud_settling_by_damping.csv")
    gammas = [5, 20, 50, 100]
    gamma_colors = {
        5: "#0b5fa5",
        20: "#1b9e77",
        50: "#d95f02",
        100: "#b22222",
    }

    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    for gamma in gammas:
        data = read_csv(RESULTS / "science_bonus" / f"bounce_gamma_{gamma}" / "diagnostics.csv")
        ax.plot(
            data["time"],
            data["center_of_mass_z"],
            linewidth=1.8,
            color=gamma_colors[gamma],
            label=fr"$\gamma_n={gamma}$",
        )
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Particle height z (m)")
    ax.set_title("Single-particle height history for varying damping")
    ax.legend(frameon=True, ncol=2)
    savefig(fig, "damping_height_vs_time.png")

    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    for gamma in gammas:
        data = read_csv(RESULTS / "science_bonus" / f"bounce_gamma_{gamma}" / "diagnostics.csv")
        ax.plot(
            data["time"],
            data["kinetic_energy"],
            linewidth=1.8,
            color=gamma_colors[gamma],
            label=fr"$\gamma_n={gamma}$",
        )
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Kinetic energy (J)")
    ax.set_title("Single-particle kinetic energy history for varying damping")
    ax.legend(frameon=True, ncol=2)
    savefig(fig, "damping_ke_vs_time.png")

    fig, ax = plt.subplots(figsize=(5.8, 3.8))
    ax.plot(scan["gamma_n"], scan["final_ke"], marker="o", linewidth=1.8, label="Single particle")
    ax.plot(cloud["gamma_n"], cloud["final_ke"], marker="s", linewidth=1.8, label="Particle cloud")
    ax.set_xlabel("Damping coefficient $\\gamma_n$")
    ax.set_ylabel("Final kinetic energy")
    ax.set_title("Effect of damping on energy dissipation")
    ax.legend(frameon=True)
    savefig(fig, "damping_vs_energy.png")

    fig, ax = plt.subplots(figsize=(5.8, 3.8))
    rebound = []
    labels = []
    for gamma in gammas:
        data = read_csv(RESULTS / "science_bonus" / f"bounce_gamma_{gamma}" / "diagnostics.csv")
        peaks = _find_bounce_peaks(data["center_of_mass_z"].to_numpy())
        rebound.append(peaks[:4] if len(peaks) >= 4 else np.pad(peaks, (0, max(0, 4 - len(peaks))), constant_values=np.nan))
        labels.append(gamma)
    rebound_arr = np.array(rebound)
    for i in range(rebound_arr.shape[0]):
        ax.plot(np.arange(1, rebound_arr.shape[1] + 1), rebound_arr[i], marker="o", linewidth=1.6, label=f"$\\gamma_n={labels[i]}$")
    ax.set_xlabel("Bounce number")
    ax.set_ylabel("Rebound peak height (m)")
    ax.set_title("Rebound height decay with damping")
    ax.legend(frameon=True, ncol=2)
    savefig(fig, "rebound_height_vs_bounce_number.png")


def plot_neighbor_bonus() -> None:
    df = read_csv(RESULTS / "neighbor_bonus" / "neighbor_search_comparison.csv")

    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    for mode, group in df.groupby("mode"):
        ax.plot(group["count"], group["runtime"], marker="o", linewidth=1.8, label=mode.replace("_", " "))
    ax.set_xlabel("Particle count")
    ax.set_ylabel("Wall-clock runtime (s)")
    ax.set_title("Neighbor-search runtime comparison")
    ax.legend(frameon=True)
    savefig(fig, "neighbor_runtime_comparison.png")

    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    for mode, group in df.groupby("mode"):
        ax.plot(group["count"], group["contact_stage_time"], marker="o", linewidth=1.8, label=mode.replace("_", " "))
    ax.set_xlabel("Particle count")
    ax.set_ylabel("Contact-stage time (s)")
    ax.set_title("Neighbor-search contact-stage runtime")
    ax.legend(frameon=True)
    savefig(fig, "neighbor_contact_stage_runtime.png")

    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    for mode, group in df.groupby("mode"):
        ax.plot(group["count"], group["candidate_pairs"], marker="o", linewidth=1.8, label=mode.replace("_", " "))
    ax.set_xlabel("Particle count")
    ax.set_ylabel("Candidate pairs examined")
    ax.set_title("Neighbor-search candidate-pair reduction")
    ax.legend(frameon=True)
    savefig(fig, "neighbor_candidate_pairs.png")

    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    for mode, group in df.groupby("mode"):
        linestyle = "--" if mode == "all_pairs" else "-"
        marker = "s" if mode == "all_pairs" else "o"
        ax.plot(
            group["count"],
            group["active_contacts"],
            marker=marker,
            linestyle=linestyle,
            linewidth=1.8,
            label=mode.replace("_", " "),
        )
    ax.set_xlabel("Particle count")
    ax.set_ylabel("Active contacts detected")
    ax.set_title("Detected-contact correctness check")
    ax.legend(frameon=True)
    savefig(fig, "neighbor_active_contacts.png")


def plot_scaling() -> None:
    strong_1000 = read_csv(RESULTS / "openmp_scaling_n1000" / "speedup.csv")
    strong_5000 = read_csv(RESULTS / "openmp_scaling_n5000" / "speedup.csv")
    weak = read_csv(RESULTS / "weak_scaling" / "weak_scaling.csv")

    fig, ax = plt.subplots(figsize=(5.8, 3.8))
    ax.plot(strong_1000["threads"], strong_1000["speedup"], marker="o", markersize=5, linewidth=2.1, color="#0b5fa5", label="Measured speedup (N=1000)")
    ax.plot(strong_5000["threads"], strong_5000["speedup"], marker="s", markersize=5, linewidth=2.1, color="#2d6a4f", label="Measured speedup (N=5000)")
    ax.plot(
        strong_5000["threads"],
        strong_5000["threads"],
        linestyle="--",
        marker="^",
        markersize=4.5,
        linewidth=1.8,
        color="#b22222",
        label="Ideal linear speedup",
    )
    ax.set_xlabel("Threads")
    ax.set_ylabel("Speedup")
    ax.set_title("OpenMP speedup")
    ax.legend(frameon=True)
    savefig(fig, "openmp_speedup.png")

    fig, ax = plt.subplots(figsize=(5.8, 3.8))
    ax.plot(strong_1000["threads"], strong_1000["efficiency"], marker="o", linewidth=1.8, color="#0b5fa5", label="N=1000")
    ax.plot(strong_5000["threads"], strong_5000["efficiency"], marker="s", linewidth=1.8, color="#2d6a4f", label="N=5000")
    ax.set_xlabel("Threads")
    ax.set_ylabel("Efficiency")
    ax.set_title("OpenMP parallel efficiency")
    ax.legend(frameon=True)
    savefig(fig, "openmp_efficiency.png")

    fig, ax = plt.subplots(figsize=(5.8, 3.8))
    baseline = weak["runtime"].iloc[0]
    ax.plot(
        weak["threads"],
        weak["runtime"],
        marker="o",
        markersize=5,
        linewidth=2.0,
        color="#0b5fa5",
        label="Measured weak-scaling runtime",
    )
    ax.axhline(baseline, linestyle="--", linewidth=1.8, color="#b22222", label="Ideal constant runtime")
    ax.set_xlabel("Threads")
    ax.set_ylabel("Runtime (s)")
    ax.set_title("OpenMP weak scaling")
    ax.legend(frameon=True)
    savefig(fig, "openmp_weak_scaling.png")

    fig, ax = plt.subplots(figsize=(5.8, 3.8))
    ax.plot(
        weak["threads"],
        weak["weak_efficiency"],
        marker="o",
        markersize=5,
        linewidth=2.0,
        color="#2d6a4f",
    )
    ax.set_xlabel("Threads")
    ax.set_ylabel("Weak-scaling efficiency")
    ax.set_title("Weak-scaling efficiency")
    savefig(fig, "openmp_weak_efficiency.png")

    fig, ax = plt.subplots(figsize=(5.8, 3.8))
    ax.plot(
        weak["threads"],
        weak["overhead"],
        marker="o",
        markersize=5,
        linewidth=2.0,
        color="#b22222",
    )
    ax.set_xlabel("Threads")
    ax.set_ylabel("Overhead (s)")
    ax.set_title("Weak-scaling overhead")
    savefig(fig, "openmp_weak_overhead.png")

    counts = [1000, 5000]
    fig, axes = plt.subplots(2, 3, figsize=(8.8, 5.6), sharex="col")
    colors = {"v_x": "#0b5fa5", "v_y": "#1b9e77", "v_z": "#d95f02"}

    for row_idx, count in enumerate(counts):
        serial_dir = RESULTS / "parallel_correctness" / f"serial_reference_n{count}"
        parallel_dir = RESULTS / "parallel_correctness" / f"parallel_reference_n{count}"
        serial_snaps = sorted(serial_dir.glob("particles_step_*.csv"))
        parallel_snaps = sorted(parallel_dir.glob("particles_step_*.csv"))

        times = []
        serial_data = {"v_x": [], "v_y": [], "v_z": []}
        parallel_data = {"v_x": [], "v_y": [], "v_z": []}

        for s_path, p_path in zip(serial_snaps, parallel_snaps):
            step = int(s_path.stem.split("_")[-1])
            times.append(step * 1.0e-4)
            s_row = read_csv(s_path).iloc[0]
            p_row = read_csv(p_path).iloc[0]
            serial_data["v_x"].append(s_row["vx"])
            serial_data["v_y"].append(s_row["vy"])
            serial_data["v_z"].append(s_row["vz"])
            parallel_data["v_x"].append(p_row["vx"])
            parallel_data["v_y"].append(p_row["vy"])
            parallel_data["v_z"].append(p_row["vz"])

        for col_idx, component in enumerate(["v_x", "v_y", "v_z"]):
            ax = axes[row_idx, col_idx]
            ax.plot(times, serial_data[component], linewidth=1.9, color=colors[component], label="Serial")
            ax.plot(
                times,
                parallel_data[component],
                linestyle="--",
                linewidth=1.7,
                color="black",
                label="Parallel",
            )
            if row_idx == 0:
                ax.set_title(component)
            if col_idx == 0:
                ax.set_ylabel(f"N={count}\nVelocity (m/s)")
            if row_idx == len(counts) - 1:
                ax.set_xlabel("Time (s)")
            if row_idx == 0 and col_idx == 0:
                ax.legend(frameon=True, loc="best")

    fig.suptitle("Serial versus OpenMP velocity-component histories", y=1.02)
    savefig(fig, "parallel_velocity_comparison.png")


def plot_particle_count_study() -> None:
    counts = [200, 1000, 5000]
    runtimes = []
    contact_times = []
    candidate_pairs = []
    final_ke = []
    max_speed = []

    for count in counts:
        summary_path = RESULTS / "profiling" / f"n{count}_serial_pairs" / "summary.txt"
        diag_path = RESULTS / "profiling" / f"n{count}_serial_pairs" / "diagnostics.csv"

        values: dict[str, float] = {}
        with summary_path.open() as f:
            for line in f:
                if ":" not in line:
                    continue
                key, value = line.split(":", 1)
                values[key.strip()] = float(value.strip())

        diag = read_csv(diag_path)
        last = diag.iloc[-1]

        runtimes.append(values["runtime"])
        contact_times.append(values["contact_time"])
        candidate_pairs.append(values["candidate_pairs"])
        final_ke.append(last["kinetic_energy"])
        max_speed.append(last["max_speed"])

    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    ax.plot(
        counts,
        runtimes,
        marker="o",
        markersize=5,
        linewidth=2.2,
        color="#0b5fa5",
        label="Total runtime",
    )
    ax.plot(
        counts,
        contact_times,
        linestyle="--",
        marker="s",
        markersize=5,
        linewidth=2.2,
        color="#b22222",
        label="Contact stage",
    )
    ax.set_xlabel("Particle count N")
    ax.set_ylabel("Runtime (s)")
    ax.set_title("Effect of increasing particle count on runtime")
    ax.legend(frameon=True)
    savefig(fig, "particle_count_runtime.png")

    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    ax.plot(counts, candidate_pairs, marker="o", linewidth=2.0, color="#1b9e77")
    ax.set_xlabel("Particle count N")
    ax.set_ylabel("Candidate pairs")
    ax.set_title("Growth of pair checks with particle count")
    savefig(fig, "particle_count_candidate_pairs.png")

    fig, axes = plt.subplots(2, 1, figsize=(6.0, 6.0), sharex=True)
    axes[0].plot(counts, final_ke, marker="o", linewidth=2.0, color="#0b5fa5")
    axes[0].set_ylabel("Final kinetic energy")
    axes[0].set_title("System response for increasing particle count")
    axes[1].plot(counts, max_speed, marker="s", linewidth=2.0, color="#d95f02")
    axes[1].set_xlabel("Particle count N")
    axes[1].set_ylabel("Final max speed (m/s)")
    savefig(fig, "particle_count_system_response.png")


def plot_compiler_optimization() -> None:
    counts = [200, 1000, 5000]
    runtime_o0 = []
    runtime_o3 = []
    speedup = []

    for count in counts:
        o0_path = RESULTS / "compiler_opt_O0" / "experiments" / f"n{count}_serial_pairs" / "summary.txt"
        o3_path = RESULTS / "compiler_opt_O3" / "experiments" / f"n{count}_serial_pairs" / "summary.txt"

        values_o0: dict[str, float] = {}
        values_o3: dict[str, float] = {}
        with o0_path.open() as f:
            for line in f:
                if ":" not in line:
                    continue
                key, value = line.split(":", 1)
                values_o0[key.strip()] = float(value.strip())
        with o3_path.open() as f:
            for line in f:
                if ":" not in line:
                    continue
                key, value = line.split(":", 1)
                values_o3[key.strip()] = float(value.strip())

        runtime_o0.append(values_o0["runtime"])
        runtime_o3.append(values_o3["runtime"])
        speedup.append(values_o0["runtime"] / values_o3["runtime"])

    fig, axes = plt.subplots(2, 1, figsize=(6.0, 6.0), sharex=True)
    axes[0].plot(counts, runtime_o0, marker="o", linewidth=2.0, color="#b22222", label="-O0")
    axes[0].plot(counts, runtime_o3, marker="s", linewidth=2.0, color="#0b5fa5", label="-O3")
    axes[0].set_ylabel("Runtime (s)")
    axes[0].set_title("Effect of compiler optimization")
    axes[0].legend(frameon=True)

    axes[1].plot(counts, speedup, marker="o", linewidth=2.0, color="#2d6a4f")
    axes[1].set_xlabel("Particle count N")
    axes[1].set_ylabel(r"Runtime ratio $T_{O0} / T_{O3}$")
    axes[1].set_title("Measured optimization speedup")
    savefig(fig, "compiler_optimization_effect.png")


def plot_cloud_settling_snapshot() -> None:
    snapshot_dir = RESULTS / "science_bonus" / "cloud_settling"
    snapshots = sorted(snapshot_dir.glob("particles_step_*.csv"))
    if not snapshots:
        return
    diag = read_csv(snapshot_dir / "diagnostics.csv")

    fig, axes = plt.subplots(2, 1, figsize=(6.0, 5.8), sharex=True)
    axes[0].plot(diag["time"], diag["center_of_mass_z"], linewidth=2.0, color="#0b5fa5")
    axes[0].set_ylabel("COM height z (m)")
    axes[0].set_title("Particle-cloud settling: center-of-mass evolution")
    axes[0].grid(True, alpha=0.35)

    axes[1].plot(diag["time"], diag["kinetic_energy"], linewidth=2.0, color="#b22222")
    axes[1].set_xlabel("Time (s)")
    axes[1].set_ylabel("Total kinetic energy (J)")
    axes[1].set_title("Particle-cloud settling: kinetic-energy evolution")
    axes[1].grid(True, alpha=0.35)
    savefig(fig, "cloud_settling_diagnostics.png")

    fig, axes = plt.subplots(1, min(4, len(snapshots)), figsize=(10.0, 2.8), sharex=True, sharey=True)
    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])
    sample_indices = np.linspace(0, len(snapshots) - 1, num=len(axes), dtype=int)
    for ax, snap_idx in zip(axes, sample_indices):
        snap = read_csv(snapshots[snap_idx])
        step = int(snapshots[snap_idx].stem.split("_")[-1])
        ax.scatter(snap["x"], snap["z"], s=8, c=snap["y"], cmap="viridis", vmin=0.0, vmax=1.0, alpha=0.85)
        time_row = diag.iloc[min(step, len(diag) - 1)]
        ax.set_title(f"t = {time_row['time']:.3f} s")
        ax.set_xlabel("x (m)")
    axes[0].set_ylabel("z (m)")
    fig.suptitle("Particle-cloud configurations at selected times", y=1.05)
    savefig(fig, "cloud_settling_configurations.png")

    cloud = read_csv(RESULTS / "science_bonus" / "cloud_settling_by_damping.csv")
    fig, axes = plt.subplots(2, 1, figsize=(5.8, 5.8), sharex=True)
    axes[0].plot(cloud["gamma_n"], cloud["final_ke"], marker="o", linewidth=1.8, color="#b22222")
    axes[0].set_ylabel("Final KE")
    axes[0].set_title("Cloud-settling response versus damping")
    axes[1].plot(cloud["gamma_n"], cloud["com_z"], marker="s", linewidth=1.8, color="#0b5fa5")
    axes[1].set_xlabel("Damping coefficient $\\gamma_n$")
    axes[1].set_ylabel("Final COM height z (m)")
    savefig(fig, "cloud_damping_response.png")

    df = read_csv(snapshots[-1])

    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    sc = ax.scatter(df["x"], df["z"], c=df["y"], s=10, cmap="viridis", alpha=0.85)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("z (m)")
    ax.set_title("Particle cloud settling snapshot")
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("y (m)")
    savefig(fig, "cloud_settling_snapshot.png")


def main() -> None:
    setup_style()
    plot_free_fall()
    plot_timestep_sensitivity()
    plot_constant_velocity()
    plot_bounce_and_damping()
    plot_particle_count_study()
    plot_compiler_optimization()
    plot_neighbor_bonus()
    plot_scaling()
    plot_cloud_settling_snapshot()
    print(f"saved figures to {FIGURES}")


if __name__ == "__main__":
    main()
