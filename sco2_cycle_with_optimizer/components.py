
import math
from properties import h_pT, s_pT, T_ps, cp_pT

def compressor_outlet_state(fluid, p_in, T_in, p_out, eta_isentropic):
    h_in = h_pT(fluid, p_in, T_in)
    s_in = s_pT(fluid, p_in, T_in)
    T_out_ideal = T_ps(fluid, p_out, s_in)
    h_out_ideal = h_pT(fluid, p_out, T_out_ideal)
    h_out = h_in + (h_out_ideal - h_in) / eta_isentropic
    cp = cp_pT(fluid, p_out, T_out_ideal)
    T_out = T_out_ideal + (h_out - h_out_ideal) / cp
    work_specific = h_out - h_in
    return T_out, h_out, work_specific

def turbine_outlet_state(fluid, p_in, T_in, p_out, eta_isentropic):
    h_in = h_pT(fluid, p_in, T_in)
    s_in = s_pT(fluid, p_in, T_in)
    T_out_ideal = T_ps(fluid, p_out, s_in)
    h_out_ideal = h_pT(fluid, p_out, T_out_ideal)
    h_out = h_in - eta_isentropic*(h_in - h_out_ideal)
    cp = cp_pT(fluid, p_out, T_out_ideal)
    T_out = T_out_ideal + (h_out - h_out_ideal) / cp
    work_specific = h_in - h_out
    return T_out, h_out, work_specific

def recuperator_counterflow(fluid_hot, fluid_cold, p_hot_in, T_hot_in, p_hot_out_target,
                            p_cold_in, T_cold_in, p_cold_out_target, mdot_hot, mdot_cold,
                            N_nodes=50, min_dT_allowed=3.0):
    cp_hot = cp_pT(fluid_hot, p_hot_in, T_hot_in)
    cp_cold = cp_pT(fluid_cold, p_cold_in, T_cold_in)
    C_hot = mdot_hot * cp_hot
    C_cold = mdot_cold * cp_cold
    Cmin = min(C_hot, C_cold)
    NTU = 5.0
    UA = NTU * Cmin
    Th = T_hot_in
    Tc = T_cold_in
    pinch_min = float('inf')
    for i in range(N_nodes):
        cp_h = cp_pT(fluid_hot, p_hot_in, Th)
        cp_c = cp_pT(fluid_cold, p_cold_in, Tc)
        Ch = mdot_hot * cp_h
        Cc = mdot_cold * cp_c
        Cmin_loc = min(Ch, Cc)
        dQ = UA / N_nodes * (Th - Tc) * (Cmin_loc / Cmin) * 1.0
        max_possible = Ch * (Th - 1.0)
        dQ = min(dQ, max_possible)
        Th_new = Th - dQ / Ch
        Tc_new = Tc + dQ / Cc
        pinch_min = min(pinch_min, Th_new - Tc_new)
        Th, Tc = Th_new, Tc_new
    pinch_violated = pinch_min < min_dT_allowed
    return Th, Tc, pinch_min, pinch_violated
