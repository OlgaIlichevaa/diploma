from math import pi, ceil


# todo: Add interface and initialize constants
path = '/home/darklord/PycharmProjects/diploma/graph/'

# adsorber
h = 2
d_sloy = 1
T0 = 308
ads_begin_vlaga = 0.32
delta_x = 0.01
delta_t = 100

# adsorbent
d = 0.004
p_nasypnoe = 800
por = 0.682
C_ads = 920
Lya = 0.0244
rho = 2200

# gas
t1 = 308
Cp = 1000
G = 1.2
vozd_begin_vlaga = 3.5
p = 1000000

# regeneration
G_reg = 1.2
C_reg = 3.5
T_gas_reg = 423
p_reg = 1000000

# calculating
s = (1 / 4) * pi * d_sloy * d_sloy
s_zer = (1 / 4) * pi * d * d
s_pov_gran = pi * d * d
n_zer = (s * 0.785) / ((d / delta_x) * s_zer)
n = ceil(h / delta_x)
v_sl = s * delta_x
m_sl = v_sl * p_nasypnoe

G_obyem = G * 1000000 / p
v = G_obyem / s

v_reg = G_reg / s

# true constants
total_time = 4 * 3600

#
c_wod = 4200
Rgaspost = 287
k_ads = 2600
A = 0.5
D_dif = 0.0001
eps = 0.9
D_dif_gran = 0.0001
dmax = 0.01

