# run_optimization.py
import math
import csv

from cycle_solver import solve_cycle
from properties import to_SI_T

# =========================
# USER INPUTS
# =========================
Q_available = 10e6          # 10 MW from helium IHX cycle
P_low_MPa = 8.0             # Chosen above CO2 critical pressure to avoid two-phase behavior
P_low = P_low_MPa * 1e6

eta_c = 0.90                # Chosen isentropic efficiency for compressor
eta_t = 0.90                # Chosen isentropic efficiency for turbine
T_cold_in_C = 35.0          # Chosen above CO2 critical temperature (31C) to avoid two-phase behavior (typically 27-47C)

# Design space
P_high_list_MPa = [21, 24, 27, 30, 33, 35]      # try different compressor outlet pressures (P2)
                                                # Chosen based on mechanical design limits (typical range 20-35MPa)

T_turbine_list_C = [650, 700, 750, 800, 850, 900]    # try different turbine inlet temperatures (T3)
                                                # Chosen based on mechanical design limits
 
mdot_min = 0.1
mdot_max = 300.0        # large range to ensure actual mdot lies within the range during bisection search

# =========================
# SOLVER
# =========================
def solve_for_mdot(P_high, T_turbine):
    def residual(mdot):
        res = solve_cycle(
            "CO2",
            P_high,
            P_low,
            to_SI_T(T_turbine),
            to_SI_T(T_cold_in_C),
            mdot,
            eta_c,
            eta_t
        )
        return mdot * res["q_in_Jkg"] - Q_available, res

    r_lo, _ = residual(mdot_min)
    r_hi, _ = residual(mdot_max)

    if r_lo * r_hi > 0:
        return None

    lo, hi = mdot_min, mdot_max
    res = None

    for _ in range(50):
        mid = 0.5 * (lo + hi)
        r_mid, res_mid = residual(mid)
        if abs(r_mid) < 1e-3:
            res = res_mid
            mdot = mid
            break
        if r_mid * r_lo < 0:
            hi = mid
        else:
            lo = mid
        res = res_mid
        mdot = mid

    net_power = mdot * res["net_w_Jkg"]
    eta = net_power / Q_available

    return {
        "P_high_MPa": P_high / 1e6,
        "P_low_MPa": P_low_MPa,
        "T_turbine_in_C": T_turbine,
        "mdot_kg_s": mdot,
        "net_power_MW": net_power / 1e6,
        "eta": eta,
        "T1_C": res["T1_K"] - 273.15,
        "T2_C": res["T2_K"] - 273.15,
        "T3_C": res["T3_K"] - 273.15,
        "T4_C": res["T4_K"] - 273.15
    }

# =========================
# RUN GRID
# =========================
results = []

for P_high_MPa in P_high_list_MPa:
    for Tt in T_turbine_list_C:
        out = solve_for_mdot(P_high_MPa * 1e6, Tt)
        if out:
            results.append(out)

# =========================
# SAVE CSV
# =========================
with open("optimization_results.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=results[0].keys())
    writer.writeheader()
    writer.writerows(results)

best = max(results, key=lambda r: r["eta"])
print("Best design point:")
for k, v in best.items():
    print(f"{k:20s} : {v}")
