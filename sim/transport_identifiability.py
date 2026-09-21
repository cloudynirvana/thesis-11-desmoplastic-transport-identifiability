#!/usr/bin/env python3
"""Three-shell transport versus a lumped burden ODE.

Seeded structural and Fisher ranks for a fibrotic delivery barrier.
Research sketch only. Not a device, a dose, or a patient measurement.

Seed 20260921.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

SEED = 20260921
RNG = np.random.default_rng(SEED)

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

# Declared constants (known, not estimated).
TAU = 4.0  # h, plasma decay
K_OUT = 0.25  # 1/h, interstitial clearance
P_V = 20.0  # mmHg, microvascular pressure
S_F = 0.015  # 1/h per mmHg, filtration scale
S_G = 1.20  # 1/h, radial coupling scale
T_END = 24.0
TIMES = np.linspace(1.0, T_END, 24)
NAMES = ("kappa", "phi", "alpha", "p_if")

# Truth: a stroma-rich setting with elevated interstitial pressure.
THETA = np.array([1.0, 0.40, 0.45, 12.0])  # kappa, phi, alpha, p_if

# Observation noise used only to scale the Fisher matrix.
# Rank itself is read from singular values of the sensitivity, before scaling.
SIGMA_C = 0.008  # concentration units
SIGMA_P = 1.0  # mmHg
SIGMA_A = 0.05  # anisotropy index


def plasma(t):
    return np.exp(-np.asarray(t, dtype=float) / TAU)


def conductances(theta):
    """Map barrier coordinates to the two conductances that enter the ODE."""
    kappa, phi, alpha, p_if = [float(v) for v in theta]
    dp = P_V - p_if
    if dp <= 0.0 or not (0.0 < phi < 1.0) or not (0.0 <= alpha < 1.0) or kappa <= 0.0:
        return np.nan, np.nan
    filt = S_F * kappa * dp * (1.0 - phi)
    diff = S_G * (1.0 - phi) ** 2 * (1.0 - alpha)
    return filt, diff


def conductance_jacobian(theta):
    """Analytic 2 x 4 derivative of (F, G) with respect to theta."""
    kappa, phi, alpha, p_if = [float(v) for v in theta]
    dp = P_V - p_if
    one_m_phi = 1.0 - phi
    one_m_alpha = 1.0 - alpha
    filt = S_F * kappa * dp * one_m_phi
    diff = S_G * one_m_phi**2 * one_m_alpha
    jac = np.zeros((2, 4))
    jac[0, 0] = filt / kappa
    jac[0, 1] = -filt / one_m_phi
    jac[0, 2] = 0.0
    jac[0, 3] = -filt / dp
    jac[1, 0] = 0.0
    jac[1, 1] = -2.0 * diff / one_m_phi
    jac[1, 2] = -diff / one_m_alpha
    jac[1, 3] = 0.0
    return jac, filt, diff


def rhs_spatial(_t, y, filt, diff):
    cr, cm, cc = y
    cp = float(plasma(_t))
    dcr = filt * (cp - cr) + diff * (cm - cr) - K_OUT * cr
    dcm = diff * (cr - 2.0 * cm + cc) - K_OUT * cm
    dcc = diff * (cm - cc) - K_OUT * cc
    return (dcr, dcm, dcc)


def rhs_lumped(_t, y, k_del):
    c = float(y[0])
    cp = float(plasma(_t))
    return (k_del * (cp - c) - K_OUT * c,)


def simulate_spatial(theta, times=TIMES):
    filt, diff = conductances(theta)
    if not np.isfinite(filt):
        return None
    sol = solve_ivp(
        rhs_spatial,
        (0.0, float(times[-1])),
        (0.0, 0.0, 0.0),
        t_eval=np.asarray(times, dtype=float),
        args=(filt, diff),
        method="BDF",
        rtol=1e-8,
        atol=1e-10,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y  # shape (3, n_t)


def simulate_lumped(k_del, times=TIMES):
    sol = solve_ivp(
        rhs_lumped,
        (0.0, float(times[-1])),
        (0.0,),
        t_eval=np.asarray(times, dtype=float),
        args=(float(k_del),),
        method="BDF",
        rtol=1e-8,
        atol=1e-10,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y[0]


def mean_conc(shells):
    return shells.mean(axis=0)


def fit_lumped(target_mean, times=TIMES):
    def residual(x):
        return simulate_lumped(x[0], times) - target_mean

    fit = least_squares(residual, np.array([0.05]), bounds=(1e-6, 5.0), ftol=1e-14, xtol=1e-14)
    pred = simulate_lumped(fit.x[0], times)
    rmse = float(np.sqrt(np.mean((pred - target_mean) ** 2)))
    return float(fit.x[0]), pred, rmse


def exact_pressure_alias(theta, p_new):
    """Keep F and G fixed by raising kappa as interstitial pressure rises."""
    kappa, phi, alpha, p_if = [float(v) for v in theta]
    dp0 = P_V - p_if
    dp1 = P_V - p_new
    kappa_new = kappa * dp0 / dp1
    return np.array([kappa_new, phi, alpha, p_new])


def exact_collagen_alias(theta, phi_new):
    """Keep F, G, and p_if fixed. phi trades against kappa and alpha."""
    kappa, phi, alpha, p_if = [float(v) for v in theta]
    g_factor = (1.0 - phi) ** 2 * (1.0 - alpha)
    one = 1.0 - phi_new
    kappa_new = kappa * (1.0 - phi) / one
    alpha_new = 1.0 - g_factor / one**2
    return np.array([kappa_new, phi_new, alpha_new, p_if])


def stacked_output(theta, schedule, times=TIMES):
    """Return the observation vector for a named schedule."""
    shells = simulate_spatial(theta, times)
    pieces = []
    if schedule in ("L", "P", "PI", "PIA", "AUC"):
        if schedule == "L":
            pieces.append(mean_conc(shells))
        elif schedule == "AUC":
            pieces.append(np.array([np.trapezoid(mean_conc(shells), times)]))
        else:
            pieces.append(shells.reshape(-1))
    if schedule in ("PI", "PIA"):
        pieces.append(np.array([float(theta[3])]))
    if schedule == "PIA":
        pieces.append(np.array([float(theta[2])]))
    if schedule == "I":
        pieces.append(np.array([float(theta[3])]))
    return np.concatenate(pieces)


def noise_vector(schedule, times=TIMES):
    n_t = len(times)
    if schedule == "L":
        return np.full(n_t, SIGMA_C)
    if schedule == "AUC":
        # Independent noise on the mean, propagated through trapezoidal weights.
        dt = np.diff(times)
        weights = np.zeros(n_t)
        weights[0] = dt[0] / 2.0
        weights[-1] = dt[-1] / 2.0
        weights[1:-1] = 0.5 * (dt[:-1] + dt[1:])
        return np.array([SIGMA_C * float(np.sqrt(np.sum(weights**2)))])
    if schedule == "P":
        return np.full(3 * n_t, SIGMA_C)
    if schedule == "I":
        return np.array([SIGMA_P])
    if schedule == "PI":
        return np.concatenate([np.full(3 * n_t, SIGMA_C), [SIGMA_P]])
    if schedule == "PIA":
        return np.concatenate([np.full(3 * n_t, SIGMA_C), [SIGMA_P, SIGMA_A]])
    raise KeyError(schedule)


def sensitivity(theta, schedule, step=1e-5):
    """Central differences of the observation map. Shape (n_obs, 4)."""
    base = stacked_output(theta, schedule)
    cols = []
    for j in range(4):
        h = step * max(abs(theta[j]), 1.0)
        plus = theta.copy()
        minus = theta.copy()
        plus[j] += h
        minus[j] -= h
        cols.append((stacked_output(plus, schedule) - stacked_output(minus, schedule)) / (2.0 * h))
    return np.column_stack(cols), base


def fisher_report(theta, schedule):
    sens, _ = sensitivity(theta, schedule)
    sigma = noise_vector(schedule)
    weighted = sens / sigma[:, None]
    # Relative columns: parameter-scaled, then noise-weighted.
    rel = weighted * theta[None, :]
    _u, s_raw, vt = np.linalg.svd(rel, full_matrices=True)
    singular = np.zeros(4)
    singular[: len(s_raw)] = np.clip(s_raw, 0.0, None)
    smax = float(singular[0]) if singular[0] > 0 else 1.0
    structural_rank = int(np.sum(singular > 1e-8 * smax))
    # Practical rank: principal-axis relative SE below one half.
    # SE = 1/singular value on the noise-weighted, parameter-scaled map.
    # Fisher eigenvalues on the relative scale (same nonzero spectrum as s^2).
    eig = singular**2
    # Analytic conductance rank is a property of the map, independent of schedule.
    jac, filt, diff = conductance_jacobian(theta)
    jac_rank = int(np.linalg.matrix_rank(jac, tol=1e-10))
    # Null-space basis of the relative sensitivity, for the write-up.
    null_tol = 1e-8 * smax
    null_basis = vt[singular <= null_tol].tolist()
    rel_se = []
    for val in eig:
        if val > 1e-8 * (smax**2):
            rel_se.append(float(1.0 / np.sqrt(val)))
        else:
            rel_se.append(None)
    practical_rank = int(sum(se is not None and se < 0.5 for se in rel_se))
    return {
        "schedule": schedule,
        "n_obs": int(sens.shape[0]),
        "structural_rank": structural_rank,
        "practical_rank": practical_rank,
        "n_param": 4,
        "singular_values_relative": [float(v) for v in singular],
        "singular_value_ratios": [float(v / smax) for v in singular],
        "fisher_eigenvalues_relative": [float(v) for v in eig],
        "approx_relative_se": rel_se,
        "null_basis_rows_of_V": null_basis,
        "conductance_jacobian_rank": jac_rank,
        "F": float(filt),
        "G": float(diff),
    }


def starved_core_case(theta):
    """Raise collagen fraction and anisotropy, then retune kappa so the mean matches."""
    phi = 0.58
    alpha = 0.82
    p_if = float(theta[3])
    target = mean_conc(simulate_spatial(theta))

    def residual(x):
        trial = np.array([float(x[0]), phi, alpha, p_if])
        shells = simulate_spatial(trial)
        if shells is None:
            return np.ones_like(target) * 10.0
        return mean_conc(shells) - target

    fit = least_squares(residual, np.array([1.5]), bounds=(0.05, 8.0))
    tuned = np.array([float(fit.x[0]), phi, alpha, p_if])
    shells = simulate_spatial(tuned)
    truth = simulate_spatial(theta)
    k_del, lumped, mean_rmse = fit_lumped(mean_conc(truth))
    core_rmse = float(np.sqrt(np.mean((lumped - truth[2]) ** 2)))
    rim_rmse = float(np.sqrt(np.mean((lumped - truth[0]) ** 2)))
    return {
        "theta": tuned.tolist(),
        "F_G": [float(v) for v in conductances(tuned)],
        "mean_rmse_vs_truth": float(np.sqrt(np.mean((mean_conc(shells) - target) ** 2))),
        "core_rmse_vs_truth": float(np.sqrt(np.mean((shells[2] - truth[2]) ** 2))),
        "rim_at_8h_truth": float(np.interp(8.0, TIMES, truth[0])),
        "core_at_8h_truth": float(np.interp(8.0, TIMES, truth[2])),
        "rim_at_8h_starved": float(np.interp(8.0, TIMES, shells[0])),
        "core_at_8h_starved": float(np.interp(8.0, TIMES, shells[2])),
        "mean_at_8h_truth": float(np.interp(8.0, TIMES, target)),
        "mean_at_8h_starved": float(np.interp(8.0, TIMES, mean_conc(shells))),
        "lumped_k_del_on_truth_mean": k_del,
        "lumped_rmse_vs_mean": mean_rmse,
        "lumped_rmse_vs_core": core_rmse,
        "lumped_rmse_vs_rim": rim_rmse,
        "core_over_rim_at_8h_truth": float(
            np.interp(8.0, TIMES, truth[2]) / np.interp(8.0, TIMES, truth[0])
        ),
        "lumped_value_at_8h": float(np.interp(8.0, TIMES, lumped)),
    }


def multistart(schedule, n_starts=24):
    """Least squares from scattered starts. Dispersion among near-best fits."""
    y = stacked_output(THETA, schedule)
    sigma = noise_vector(schedule)
    # One noisy draw, so the target is not exactly the generating mean.
    y_obs = y + RNG.normal(0.0, sigma)

    lo = np.array([0.15, 0.08, 0.05, 4.0])
    hi = np.array([4.0, 0.72, 0.90, 18.0])

    def residual(x):
        pred = stacked_output(x, schedule)
        if pred is None or not np.all(np.isfinite(pred)):
            return np.ones_like(y_obs) * 1e3
        return (pred - y_obs) / sigma

    starts = [THETA.copy()]
    for _ in range(n_starts - 1):
        starts.append(lo + RNG.random(4) * (hi - lo))

    rows = []
    for x0 in starts:
        fit = least_squares(residual, x0, bounds=(lo, hi), ftol=1e-12, xtol=1e-12, gtol=1e-12, max_nfev=80)
        cost = float(np.sum(fit.fun**2))
        rows.append({"theta": fit.x.tolist(), "cost": cost, "success": bool(fit.success)})

    costs = np.array([r["cost"] for r in rows])
    best = float(costs.min())
    # Keep fits whose excess cost is within a chi-square-ish band of the best.
    # 3.84 is the 95% point of chi-square with 1 d.f.; used only as a window.
    keep = [r for r in rows if r["cost"] <= best + 3.84]
    mat = np.array([r["theta"] for r in keep])
    return {
        "n_starts": n_starts,
        "n_near_best": int(len(keep)),
        "best_cost": best,
        "cost_span_near_best": float(mat[:, 0].shape[0] and (costs[costs <= best + 3.84].max() - best)),
        "theta_min": mat.min(axis=0).tolist(),
        "theta_max": mat.max(axis=0).tolist(),
        "theta_std": mat.std(axis=0).tolist(),
        "names": list(NAMES),
    }


def step_check(theta, schedule="P"):
    """Halve the finite-difference step and compare singular values."""
    s1, _ = sensitivity(theta, schedule, step=1e-5)
    s2, _ = sensitivity(theta, schedule, step=5e-6)
    sig = noise_vector(schedule)
    r1 = np.linalg.svd((s1 / sig[:, None]) * theta, compute_uv=False)
    r2 = np.linalg.svd((s2 / sig[:, None]) * theta, compute_uv=False)
    keep = r1 > 1e-6
    rel = float(np.max(np.abs(r1[keep] - r2[keep]) / r1[keep]))
    return {
        "max_relative_gap_retained_singular_values": rel,
        "singular_1em5": r1.tolist(),
        "singular_5em6": r2.tolist(),
    }


def trajectory_dump(theta):
    shells = simulate_spatial(theta)
    return {
        "t": TIMES.tolist(),
        "rim": shells[0].tolist(),
        "mid": shells[1].tolist(),
        "core": shells[2].tolist(),
        "mean": mean_conc(shells).tolist(),
    }


def make_figures(truth_traj, lumped_pred, alias_p, alias_c, starved, spectra, path_p, path_c):
    t = TIMES

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.plot(t, truth_traj["rim"], label="rim", color="#1b4f72")
    ax.plot(t, truth_traj["mid"], label="mid-stroma", color="#148f77")
    ax.plot(t, truth_traj["core"], label="core", color="#b03a2e")
    ax.plot(t, truth_traj["mean"], label="spatial mean", color="#1a1a1a", lw=1.6)
    ax.plot(t, lumped_pred, label="lumped ODE fit to the mean", color="#7d3c98", ls="--")
    ax.set_xlabel("time (h)")
    ax.set_ylabel("concentration (arbitrary)")
    ax.set_title("Three shells and the lumped fit to their mean")
    ax.legend(frameon=False, fontsize=8)
    ax.set_xlim(1, 24)
    fig.tight_layout()
    fig.savefig(FIG / "fig_trajectories.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    truth_shells = simulate_spatial(THETA)
    alias_shells = simulate_spatial(alias_p)
    starved_shells = simulate_spatial(np.array(starved["theta"]))
    ax.plot(t, truth_shells[2], label="core, truth", color="#b03a2e")
    ax.plot(t, alias_shells[2], label="core, pressure-conductivity alias", color="#1a1a1a", ls="--")
    ax.plot(t, starved_shells[2], label="core, fibrotic retune", color="#b9770e")
    ax.plot(t, mean_conc(truth_shells), label="mean, truth", color="#5d6d7e", ls=":")
    ax.plot(t, mean_conc(starved_shells), label="mean, fibrotic retune", color="#b9770e", ls=":")
    ax.set_xlabel("time (h)")
    ax.set_ylabel("concentration (arbitrary)")
    ax.set_title("Core exposure under an alias and under a matched-mean barrier")
    ax.legend(frameon=False, fontsize=8)
    ax.set_xlim(1, 24)
    fig.tight_layout()
    fig.savefig(FIG / "fig_core_alias.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    labels = ["L mean", "AUC", "P perfusion", "PI perfusion+IFP", "PIA +anisotropy"]
    keys = ["L", "AUC", "P", "PI", "PIA"]
    x = np.arange(4)
    width = 0.15
    for i, key in enumerate(keys):
        ratios = np.array(spectra[key]["singular_value_ratios"])
        ax.bar(x + (i - 2) * width, ratios, width=width, label=labels[i])
    ax.set_xticks(x)
    ax.set_xticklabels(["s1", "s2", "s3", "s4"])
    ax.set_ylabel("singular value / largest")
    ax.set_yscale("log")
    ax.set_ylim(1e-16, 2)
    ax.set_title("Relative sensitivity spectra")
    ax.legend(frameon=False, fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(FIG / "fig_spectra.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.8), sharey=True)
    axes[0].plot(path_p["p_if"], path_p["max_abs_mean_gap"], label="mean gap", color="#1a1a1a")
    axes[0].plot(path_p["p_if"], path_p["max_abs_core_gap"], label="core gap", color="#b03a2e")
    axes[0].plot(path_p["p_if"], path_p["ifp_gap"], label="|p − p truth| / 10", color="#1b4f72", ls="--")
    axes[0].set_xlabel("alias interstitial pressure (mmHg)")
    axes[0].set_ylabel("absolute gap")
    axes[0].set_title("Pressure-conductivity level set")
    axes[0].legend(frameon=False, fontsize=7)
    axes[1].plot(path_c["alpha"], path_c["max_abs_shell_gap"], label="shell gap", color="#1a1a1a")
    axes[1].plot(path_c["alpha"], path_c["alpha_gap"], label="|α − α truth|", color="#117a65", ls="--")
    axes[1].set_xlabel("alias anisotropy α")
    axes[1].set_title("Collagen level set at fixed IFP")
    axes[1].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "fig_level_sets.png", dpi=160)
    plt.close(fig)


def level_path_pressure():
    p_grid = np.linspace(8.0, 17.5, 20)
    rows = {k: [] for k in ("p_if", "kappa", "max_abs_mean_gap", "max_abs_core_gap", "ifp_gap")}
    truth = simulate_spatial(THETA)
    for p in p_grid:
        alias = exact_pressure_alias(THETA, float(p))
        shells = simulate_spatial(alias)
        rows["p_if"].append(float(p))
        rows["kappa"].append(float(alias[0]))
        rows["max_abs_mean_gap"].append(float(np.max(np.abs(mean_conc(shells) - mean_conc(truth)))))
        rows["max_abs_core_gap"].append(float(np.max(np.abs(shells[2] - truth[2]))))
        rows["ifp_gap"].append(abs(float(p) - float(THETA[3])) / 10.0)
    return rows


def level_path_collagen():
    phi_grid = np.linspace(0.18, 0.55, 16)
    rows = {k: [] for k in ("phi", "alpha", "kappa", "max_abs_shell_gap", "alpha_gap")}
    truth = simulate_spatial(THETA)
    for phi in phi_grid:
        alias = exact_collagen_alias(THETA, float(phi))
        if not (0.0 < alias[2] < 0.98):
            continue
        shells = simulate_spatial(alias)
        rows["phi"].append(float(phi))
        rows["alpha"].append(float(alias[2]))
        rows["kappa"].append(float(alias[0]))
        rows["max_abs_shell_gap"].append(float(np.max(np.abs(shells - truth))))
        rows["alpha_gap"].append(abs(float(alias[2]) - float(THETA[2])))
    return rows


def noiseless_refit(schedule, x0):
    """Fit from a distant start to noiseless output. Recovery is a property of the map."""
    y = stacked_output(THETA, schedule)
    sigma = noise_vector(schedule)
    lo = np.array([0.15, 0.08, 0.05, 4.0])
    hi = np.array([4.0, 0.72, 0.90, 18.0])

    def residual(x):
        return (stacked_output(x, schedule) - y) / sigma

    fit = least_squares(residual, np.asarray(x0, dtype=float), bounds=(lo, hi), ftol=1e-14, xtol=1e-14, gtol=1e-14)
    return {
        "start": [float(v) for v in x0],
        "estimate": fit.x.tolist(),
        "cost": float(np.sum(fit.fun**2)),
        "abs_error": np.abs(fit.x - THETA).tolist(),
    }


def pia_coordinate_se(theta):
    """Diagonal Cramér-Rao sketch on the relative scale, full-rank schedule only."""
    sens, _ = sensitivity(theta, "PIA")
    sigma = noise_vector("PIA")
    rel = (sens / sigma[:, None]) * theta[None, :]
    fim = rel.T @ rel
    cov = np.linalg.inv(fim)
    se = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    return {"relative_coordinate_se": se.tolist(), "names": list(NAMES), "condition_number": float(np.linalg.cond(fim))}


def main():
    jac, filt, diff = conductance_jacobian(THETA)
    jac_rank = int(np.linalg.matrix_rank(jac, tol=1e-10))
    # Kernel check: H v = 0 for the two analytic local null directions.
    kappa, phi, alpha, p_if = THETA
    dp = P_V - p_if
    v_pressure = np.array([kappa, 0.0, 0.0, dp])
    v_collagen = np.array([
        kappa / (1.0 - phi),
        1.0,
        -2.0 * (1.0 - alpha) / (1.0 - phi),
        0.0,
    ])
    kernel_residuals = {
        "pressure_conductivity": float(np.linalg.norm(jac @ v_pressure)),
        "collagen_trade": float(np.linalg.norm(jac @ v_collagen)),
    }

    schedules = {}
    for key in ("L", "AUC", "P", "I", "PI", "PIA"):
        schedules[key] = fisher_report(THETA, key)

    shells = simulate_spatial(THETA)
    k_del, lumped_pred, mean_rmse = fit_lumped(mean_conc(shells))
    starved = starved_core_case(THETA)

    alias_p = exact_pressure_alias(THETA, 16.0)
    alias_c = exact_collagen_alias(THETA, 0.22)
    alias_checks = {
        "pressure_alias_theta": alias_p.tolist(),
        "pressure_alias_FG": [float(v) for v in conductances(alias_p)],
        "pressure_alias_max_abs_shell": float(
            np.max(np.abs(simulate_spatial(alias_p) - shells))
        ),
        "collagen_alias_theta": alias_c.tolist(),
        "collagen_alias_FG": [float(v) for v in conductances(alias_c)],
        "collagen_alias_max_abs_shell": float(
            np.max(np.abs(simulate_spatial(alias_c) - shells))
        ),
        "truth_FG": [float(filt), float(diff)],
    }

    path_p = level_path_pressure()
    path_c = level_path_collagen()
    fd = step_check(THETA, "P")
    distant = np.array([2.4, 0.25, 0.72, 16.5])
    refit_L = noiseless_refit("L", distant)
    refit_PIA = noiseless_refit("PIA", distant)
    coord_se = pia_coordinate_se(THETA)

    # Multistart is the practical companion. L should wander; PIA should cluster.
    print("multistart L")
    multi_L = multistart("L", n_starts=20)
    print("multistart PIA")
    multi_PIA = multistart("PIA", n_starts=16)

    truth_traj = trajectory_dump(THETA)
    # Downsample stored trajectories for the JSON (figures use the arrays directly).
    make_figures(truth_traj, lumped_pred, alias_p, alias_c, starved, schedules, path_p, path_c)

    peak = {
        "rim": float(shells[0].max()),
        "mid": float(shells[1].max()),
        "core": float(shells[2].max()),
        "mean": float(mean_conc(shells).max()),
        "rim_at_8h": float(np.interp(8.0, TIMES, shells[0])),
        "mid_at_8h": float(np.interp(8.0, TIMES, shells[1])),
        "core_at_8h": float(np.interp(8.0, TIMES, shells[2])),
        "mean_at_8h": float(np.interp(8.0, TIMES, mean_conc(shells))),
    }

    results = {
        "seed": SEED,
        "theta_truth": THETA.tolist(),
        "names": list(NAMES),
        "constants": {
            "tau_h": TAU,
            "k_out_per_h": K_OUT,
            "P_v_mmHg": P_V,
            "s_F": S_F,
            "s_G": S_G,
            "sigma_c": SIGMA_C,
            "sigma_p_mmHg": SIGMA_P,
            "sigma_alpha": SIGMA_A,
            "times_h": TIMES.tolist(),
        },
        "conductances": {"F": float(filt), "G": float(diff)},
        "conductance_jacobian": jac.tolist(),
        "conductance_jacobian_rank": jac_rank,
        "kernel_residuals": kernel_residuals,
        "kernel_vectors": {"pressure_conductivity": v_pressure.tolist(), "collagen_trade": v_collagen.tolist()},
        "finite_difference_check": fd,
        "schedules": schedules,
        "lumped": {
            "k_del": k_del,
            "rmse_vs_mean": mean_rmse,
            "rmse_vs_core": float(np.sqrt(np.mean((lumped_pred - shells[2]) ** 2))),
            "rmse_vs_rim": float(np.sqrt(np.mean((lumped_pred - shells[0]) ** 2))),
            "value_at_8h": float(np.interp(8.0, TIMES, lumped_pred)),
        },
        "peaks": peak,
        "aliases": alias_checks,
        "starved": starved,
        "level_path_pressure_max_mean_gap": float(np.max(path_p["max_abs_mean_gap"])),
        "level_path_pressure_max_core_gap": float(np.max(path_p["max_abs_core_gap"])),
        "level_path_collagen_max_shell_gap": float(np.max(path_c["max_abs_shell_gap"])),
        "multistart_L": multi_L,
        "multistart_PIA": multi_PIA,
        "noiseless_refit_L": refit_L,
        "noiseless_refit_PIA": refit_PIA,
        "pia_coordinate_se": coord_se,
        "trajectory_truth": {
            "t": truth_traj["t"],
            "rim": truth_traj["rim"],
            "mid": truth_traj["mid"],
            "core": truth_traj["core"],
            "mean": truth_traj["mean"],
            "lumped": lumped_pred.tolist(),
        },
    }
    out = ROOT / "results.json"
    out.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps({
        "F": filt,
        "G": diff,
        "jac_rank": jac_rank,
        "kernel": kernel_residuals,
        "ranks": {k: {"structural": v["structural_rank"], "practical": v["practical_rank"], "ratios": v["singular_value_ratios"]} for k, v in schedules.items()},
        "lumped": results["lumped"],
        "peaks": peak,
        "aliases": alias_checks,
        "starved_subset": {k: starved[k] for k in ("theta", "mean_rmse_vs_truth", "core_rmse_vs_truth", "core_at_8h_truth", "core_at_8h_starved", "mean_at_8h_truth", "mean_at_8h_starved", "core_over_rim_at_8h_truth")},
        "fd": fd["max_relative_gap_retained_singular_values"],
        "refit_L": refit_L,
        "refit_PIA": refit_PIA,
        "coord_se": coord_se,
        "multi_L_std": multi_L["theta_std"],
        "multi_L_n": multi_L["n_near_best"],
        "multi_PIA_std": multi_PIA["theta_std"],
        "multi_PIA_n": multi_PIA["n_near_best"],
        "multi_PIA_range": [multi_PIA["theta_min"], multi_PIA["theta_max"]],
        "multi_L_range": [multi_L["theta_min"], multi_L["theta_max"]],
    }, indent=2))


if __name__ == "__main__":
    main()
