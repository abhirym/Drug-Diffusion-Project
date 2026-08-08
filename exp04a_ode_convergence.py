import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# --------------------------------------------------
# 1. Parameters / initial conditions
# --------------------------------------------------

# Initial state vector
y0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

params = {
    "EGF" : 1.0,
    "G" : 0.0,
    "IC_50" : 0.3,
    "n" : 1.0,
    "k_EGF" : 1.0,
    "k_EGFR_off" : 0.5,
    "k_G" : 0.5,
    "k_RAS" : 1.0,
    "k_RAS_off" : 0.5,
    "k_RAF" : 1.0,
    "k_RAF_off" : 0.5,
    "k_MEK" : 1.0,
    "k_MEK_off" : 0.5,
    "k_ERK" : 1.0,
    "k_ERK_off" : 0.5,
    "k_T" : 1.0,
    "k_deg" : 0.5,
    "k_P" : 1.0,
    "k_loss" : 0.5,
}

species = [
    "EGFR*",
    "RAS*",
    "RAF*",
    "MEK*",
    "ERK*",
    "Cyclin D",
    "P"
]

T_values = [40, 80, 160, 320, 640]

RATE_TOL = 1e-6
STATE_TOL = 1e-6

def pathway(t, y, params):
    # External factors
    EGF = params["EGF"]
    G = params["G"]
    IC_50 = params["IC_50"]
    n = params["n"]

    # Extract activated protein and Cyclin D/Proliferation rates from array y
    EGFR_star, RAS_star, RAF_star, MEK_star, ERK_star, CycD, P = y

    # Total protein concentrations
    EGFR_tot = 1.0
    RAS_tot = 1.0
    RAF_tot = 1.0
    MEK_tot = 1.0
    ERK_tot = 1.0

    # Rate constants (normalized currently)
    k_EGF = params["k_EGF"]
    k_EGFR_off = params["k_EGFR_off"]
    k_G = params["k_G"]

    k_RAS = params["k_RAS"]
    k_RAS_off = params["k_RAS_off"]

    k_RAF = params["k_RAF"]
    k_RAF_off = params["k_RAF_off"]

    k_MEK = params["k_MEK"]
    k_MEK_off = params["k_MEK_off"]

    k_ERK = params["k_ERK"]
    k_ERK_off = params["k_ERK_off"]

    k_T = params["k_T"]
    k_deg = params["k_deg"]

    k_P = params["k_P"]
    k_loss = params["k_loss"]

    # Calculating derivatives
    EGFR_gefit = 1 / (1 + (G / IC_50) ** n)
    EGFR_act = k_EGF * EGF * (EGFR_tot - EGFR_star) * EGFR_gefit
    EGFR_dect = k_EGFR_off * EGFR_star

    dEGFR = EGFR_act - EGFR_dect

    RAS_act = k_RAS * EGFR_star * (RAS_tot - RAS_star)
    RAS_dect = k_RAS_off * RAS_star

    dRAS = RAS_act - RAS_dect

    RAF_act = k_RAF * RAS_star * (RAF_tot - RAF_star)
    RAF_dect = k_RAF_off * RAF_star

    dRAF = RAF_act - RAF_dect

    MEK_act = k_MEK * RAF_star * (MEK_tot - MEK_star)
    MEK_dect = k_MEK_off * MEK_star

    dMEK = MEK_act - MEK_dect

    ERK_act = k_ERK * MEK_star * (ERK_tot - ERK_star)
    ERK_dect = k_ERK_off * ERK_star

    dERK = ERK_act - ERK_dect

    CycD_act = k_T * ERK_star
    CycD_dect = k_deg * CycD

    dCycD = CycD_act - CycD_dect

    P_act = k_P * CycD
    P_dect = k_loss * P

    dP = P_act - P_dect

    return [dEGFR, dRAS, dRAF, 
                dMEK, dERK, dCycD, dP]

# --------------------------------------------------
# 2. Convergence test
# --------------------------------------------------

previous_final_state = None

for T in T_values:

    t_eval = np.linspace(0, T, 2000)

    solution = solve_ivp(
        pathway,
        [0, T],
        y0,
        args=(params,),
        t_eval=t_eval,
        rtol=1e-9,
        atol=1e-12
    )

    if not solution.success:
        raise RuntimeError(solution.message)

    # Final state
    final_state = solution.y[:, -1]

    # ----------------------------------------------
    # Criterion 1: relative rate of change
    # ----------------------------------------------

    final_derivatives = pathway(
        solution.t[-1],
        final_state,
        params
    )

    relative_rates = (
        np.abs(final_derivatives)
        / np.maximum(np.abs(final_state), 1e-12)
    )

    max_relative_rate = np.max(relative_rates)

    # ----------------------------------------------
    # Criterion 2: relative state change
    # ----------------------------------------------

    max_state_change = None

    if previous_final_state is not None:

        relative_state_change = (
            np.abs(final_state - previous_final_state)
            / np.maximum(np.abs(final_state), 1e-12)
        )

        max_state_change = np.max(relative_state_change)

    # ----------------------------------------------
    # Print results
    # ----------------------------------------------

    print(f"\nT = {T}")
    print(f"max |dy/dt| = "
          f"{np.max(np.abs(final_derivatives)):.3e}")
    print(f"max relative rate = "
          f"{max_relative_rate:.3e}")

    if max_state_change is not None:
        print(f"max relative state change = "
              f"{max_state_change:.3e}")

    # ----------------------------------------------
    # Check convergence
    # ----------------------------------------------

    rate_converged = max_relative_rate < RATE_TOL

    state_converged = (
        max_state_change is not None
        and max_state_change < STATE_TOL
    )

    if state_converged and rate_converged:
        print(">>> CONVERGED")

    previous_final_state = final_state

results = []

previous_final_state = None

for T in T_values:

    t_eval = np.linspace(0, T, 2000)

    solution = solve_ivp(
        pathway,
        [0, T],
        y0,
        args=(params,),
        t_eval=t_eval,
        rtol=1e-9,
        atol=1e-12
    )

    final_state = solution.y[:, -1]

    final_derivatives = pathway(
        solution.t[-1],
        final_state,
        params
    )

    relative_rates = (
        np.abs(final_derivatives)
        / np.maximum(np.abs(final_state), 1e-12)
    )

    max_relative_rate = np.max(relative_rates)

    if previous_final_state is None:
        max_state_change = np.nan
    else:
        relative_state_change = (
            np.abs(final_state - previous_final_state)
            / np.maximum(np.abs(final_state), 1e-12)
        )

        max_state_change = np.max(relative_state_change)

    results.append({
        "T": T,
        "max_relative_rate": max_relative_rate,
        "max_state_change": max_state_change
    })

    previous_final_state = final_state