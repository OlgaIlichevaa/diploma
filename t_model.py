from table import Table
from math import pi, ceil
from constants import c_wod, Rgaspost, k_ads, A, D_dif, dmax, total_time, delta_t
from visualization import plot_data


class TModel:
    def __init__(self, ad, aT0, aLya, ad_sloy, ap_nasypnoe, apor, a_ads_begin_vlaga, ah, arho, at1, aCp, aG,
                 avozd_begin_vlaga, ap, aC_ads, adelta_x, adelta_t, aG_reg, aC_reg, aTvoz_reg):
        self.d = ad
        self.T0 = aT0
        self.Lya = aLya
        self.d_sloy = ad_sloy
        self.p_nasypnoe = ap_nasypnoe
        self.por = apor
        self.ads_begin_vlaga = a_ads_begin_vlaga
        self.h = ah
        self.rho = arho
        self.t1 = at1
        self.Cp = aCp
        self.G = aG
        self.vozd_begin_vlaga = avozd_begin_vlaga
        self.p = ap
        self.C_ads = aC_ads
        self.delta_x = adelta_x
        self.delta_t = adelta_t
        self.G_reg = aG_reg
        self.C_reg = aC_reg
        self.Tvoz_reg = aTvoz_reg
        self.s = (1 / 4) * pi * self.d_sloy * self.d_sloy
        self.s_zer = (1 / 4) * pi * self.d * self.d
        self.s_pov_gran = pi * self.d * self.d
        self.n_zer = (self.s * 0.785) / ((self.d / self.delta_x) * self.s_zer)
        self.n = ceil(self.h / self.delta_x)
        self.v_sl = self.s * self.delta_x
        self.m_sl = self.v_sl * self.p_nasypnoe
        self.G_obyem = self.G * 1000000 / self.p
        self.v = self.G_obyem / self.s
        self.v_reg = self.G_reg / self.s

    def go(self):
        arho = self.p / (287.4 * self.t1)
        avyaz = 0.00001717 * (self.t1 / 273) ** 0.82
        aL = 0.0244 * (self.t1 / 273) ** 0.82
        aRe = arho * self.v * self.d / avyaz
        output = [[Table() for layer in range(self.n)] for dt_step in range(total_time // delta_t)]
        plot_data(output, layer=5, field='VlagosoderganieOfVozd', title='layer=5_VlagosoderganieOfVozd.png')
        # todo: output[0][0].TemperatureOfVozd and output[0][0].VlagosoderganieOfVozd initialized two times

        # start conditions
        for dt_step in range(total_time // delta_t):
            output[dt_step][0].TemperatureOfVozd = self.T0
            output[dt_step][0].VlagosoderganieOfVozd = self.vozd_begin_vlaga
        for layer in range(self.n):
            output[0][layer].TemperatureOfAds = self.T0
            output[0][layer].TemperatureOfVozd = self.T0
            output[0][layer].VlagosoderganieOfAds = self.ads_begin_vlaga
            output[0][layer].VlagosoderganieOfVozd = self.vozd_begin_vlaga
        plot_data(output, layer=5, field='VlagosoderganieOfVozd', title='layer=5_VlagosoderganieOfVozd.png')

        for dt_step in range(total_time // delta_t):
            for layer in range(1, self.n):
                if aRe < 1000:
                    aC = 0.56
                    an = 0.5
                else:
                    aC = 0.28
                    an = 0.6
                aPr = avyaz * self.Cp / aL
                aNu = aC * aRe ** an * aPr ** 0.33
                aAlf = aNu * aL / self.d
        pass
