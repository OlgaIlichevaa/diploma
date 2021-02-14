from constants import *
from t_model import TModel


model = TModel(d, T0, Lya, d_sloy, p_nasypnoe, por, ads_begin_vlaga, h, rho, t1, Cp, G, vozd_begin_vlaga, p, C_ads,
               delta_x, delta_t, G_reg, C_reg, T_gas_reg)
model.go()
