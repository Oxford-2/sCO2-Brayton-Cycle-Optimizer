
import math
try:
    from CoolProp.CoolProp import PropsSI
    COOLPROP_AVAILABLE = True
except Exception as e:
    COOLPROP_AVAILABLE = False

def to_SI_T(t_c):
    return t_c + 273.15

def h_pT(fluid, p_pa, T_k):
    if COOLPROP_AVAILABLE:
        return PropsSI('H','P',p_pa,'T',T_k,fluid)
    else:
        cp = cp_pT(fluid, p_pa, T_k)
        return cp * T_k

def s_pT(fluid, p_pa, T_k):
    if COOLPROP_AVAILABLE:
        return PropsSI('S','P',p_pa,'T',T_k,fluid)
    else:
        R = gas_constant(fluid)
        Tref = 273.15
        pref = 101325.0
        cp = cp_pT(fluid, p_pa, T_k)
        return cp * math.log(T_k / Tref) - R * math.log(p_pa / pref)

def T_ps(fluid, p_pa, s_jpkg):
    if COOLPROP_AVAILABLE:
        return PropsSI('T','P',p_pa,'S',s_jpkg,fluid)
    else:
        R = gas_constant(fluid)
        cp = cp_pT(fluid, p_pa, 300.0)
        Tref = 273.15
        pref = 101325.0
        T = Tref * math.exp((s_jpkg + R * math.log(p_pa/pref)) / cp)
        return T

def cp_pT(fluid, p_pa, T_k):
    if COOLPROP_AVAILABLE:
        return PropsSI('C','P',p_pa,'T',T_k,fluid)
    else:
        if fluid.upper() == "CO2":
            T = T_k
            cp = 800.0 + 0.1*(T-300.0)
            return cp
        else:
            return 1000.0

def rho_pT(fluid, p_pa, T_k):
    if COOLPROP_AVAILABLE:
        return PropsSI('D','P',p_pa,'T',T_k,fluid)
    else:
        R = gas_constant(fluid)
        return p_pa / (R * T_k)

def gas_constant(fluid):
    R_universal = 8.31446261815324
    mm = 0.04401 if fluid.upper()=="CO2" else 0.02897
    return R_universal / mm
