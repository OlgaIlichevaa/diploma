from math import pi, ceil


# todo: Add interface and initialize constants
path = "C:\\Users\\olgai\\PycharmProjects\\diploma\\graph_v2\\"
# /home/darklord/PycharmProjects/diploma/graph/
# adsorber
h = 2
d_sloy = 1
T0 = 308
ads_begin_vlaga = 0.32
delta_x = 0.002
delta_t = 1  #100

# adsorbent
d = 0.002
p_nasypnoe = 600 #700
por = 0.59
C_ads = 920
Lya = 0.0244
rho_ads = 1200

# gas
t1 = 308
G = 10
vozd_begin_vlaga = 3.5
p = 1000000


# regeneration
G_reg = 2.4
C_reg = 3.5
A_reg = 3.5
Tvoz_reg = 423
p_reg = 1000000
R_vozd = 287

# calculating
s = (1 / 4) * pi * d_sloy * d_sloy
s_zer = (1 / 4) * pi * d * d
s_pov_gran = pi * d * d
n_zer = (s * 0.785) / ((d / delta_x) * s_zer)
n = ceil(h / delta_x)
v_sl = s * delta_x
m_sl = v_sl * p_nasypnoe

G_obyem = G * 100000 / p
v = G_obyem / s

v_reg = G_reg / s

# true constants
total_time = 4 * 3600

#
c_wod = 4200

k_ads = 2600
A = 0.5
D_dif = 0.0001
eps = 0.682
D_dif_gran = 0.0001
dmax = 0.01

