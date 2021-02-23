from constants import *
from t_model import TModel


model = TModel(d, T0, Lya, d_sloy, p_nasypnoe, por, ads_begin_vlaga, h, rho_ads, t1, G, vozd_begin_vlaga, p,
               C_ads, delta_x, delta_t, G_reg, C_reg, A_reg, Tvoz_reg, p_reg, R_vozd)
model.go()
