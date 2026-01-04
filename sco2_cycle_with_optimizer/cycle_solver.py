
from properties import to_SI_T, h_pT, s_pT, cp_pT
from components import compressor_outlet_state, turbine_outlet_state, recuperator_counterflow

def solve_cycle(fluid, P_high_pa, P_low_pa, T_turbine_in_K, T_cold_in_K, mdot,
                eta_c=0.80, eta_t=0.85, recup_N=50):
    p1 = P_low_pa
    T1 = T_cold_in_K
    h1 = h_pT(fluid, p1, T1)
    T2, h2, w_comp = compressor_outlet_state(fluid, p1, T1, P_high_pa, eta_c)
    T3 = T_turbine_in_K
    h3 = h_pT(fluid, P_high_pa, T3)
    T4, h4, w_turb = turbine_outlet_state(fluid, P_high_pa, T3, P_low_pa, eta_t)
    Th_out, Tc_out, pinch, pinch_viol = recuperator_counterflow(fluid, fluid, P_high_pa, T4, P_low_pa,
                                                                 P_low_pa, T2, P_high_pa, mdot, mdot,
                                                                 N_nodes=recup_N)
    net_work_specific = w_turb - w_comp
    h_before_heater = h_pT(fluid, P_high_pa, Tc_out)
    q_in_specific = h3 - h_before_heater
    eta_thermal = net_work_specific / q_in_specific if q_in_specific>0 else float('nan')
    results = {
        "T1_K": T1, "p1_Pa": p1, "h1_Jkg": h1,
        "T2_K": T2, "h2_Jkg": h2,
        "T3_K": T3, "h3_Jkg": h3,
        "T4_K": T4, "h4_Jkg": h4,
        "w_comp_Jkg": w_comp, "w_turb_Jkg": w_turb,
        "net_w_Jkg": net_work_specific,
        "q_in_Jkg": q_in_specific,
        "eta_thermal": eta_thermal,
        "recup_pinch_K": pinch, "recup_pinch_violated": pinch_viol
    }
    return results
