from table import Table
from math import pi, ceil, exp
from constants import c_wod, k_ads, A, D_dif, dmax, total_time, delta_t
from visualization import plot_data, plot_all_layers


class TModel:
    def __init__(self, ad, aT0, aLya, ad_sloy, ap_nasypnoe, apor, a_ads_begin_vlaga, ah, arho_ads, at1, aG,
                 avozd_begin_vlaga, ap, aC_ads, adelta_x, adelta_t, aG_reg, aC_reg, aA_reg, aTvoz_reg,
                 ap_reg, aR_vozd):
        self.d = ad
        self.T0 = aT0
        self.Lya = aLya
        self.d_sloy = ad_sloy
        self.p_nasypnoe = ap_nasypnoe
        self.por = apor
        self.ads_begin_vlaga = a_ads_begin_vlaga
        self.h = ah
        self.rho_ads = arho_ads
        self.t1 = at1
        self.G = aG
        self.vozd_begin_vlaga = avozd_begin_vlaga
        self.p = ap
        self.C_ads = aC_ads
        self.delta_x = adelta_x
        self.delta_t = adelta_t
        self.G_reg = aG_reg
        self.C_reg = aC_reg
        self.A_reg = aA_reg
        self.Tvoz_reg = aTvoz_reg
        self.p_reg = ap_reg
        self.R_vozd = aR_vozd
        self.s = (1 / 4) * pi * self.d_sloy * self.d_sloy
        self.s_zer = (1 / 4) * pi * self.d * self.d
        self.s_pov_gran = pi * self.d * self.d
        self.n_zer = (self.s * 0.785) / ((self.d / self.delta_x) * self.s_zer)
        self.n = ceil(self.h / self.delta_x)
        self.v_sl = self.s * self.delta_x
        self.m_sl = self.v_sl * self.p_nasypnoe
        self.G_obyem = self.G * 100000 / self.p
        self.v = self.G_obyem / self.s
        self.v_reg = self.G_reg / self.s


    # adsorption

    def go(self):
        # arho_voz = self.p / (287.4 * self.t1)
        # avyaz_voz = 0.00001717 * (self.t1 / 273) ** 0.683
        # aLya_voz = 0.0244 * (self.t1 / 273) ** 0.82
        # aCp_voz = 1000
        # avyaz_kin = avyaz_voz / arho_voz
        output = [[Table() for layer in range(self.n)] for dt_step in range(int(total_time / delta_t))]
        # plot_data(output, layer=5, field='VlagosoderganieOfVozd', title='layer=5_VlagosoderganieOfVozd.png')
        # todo: output[0][0].TemperatureOfVozd and output[0][0].VlagosoderganieOfVozd initialized two times

        # start conditions adsorption
        for layer in range(self.n):
            output[0][layer].TemperatureOfAds = self.T0
            output[0][layer].TemperatureOfVozd = self.T0
            output[0][layer].VlagosoderganieOfAds = self.ads_begin_vlaga
            output[0][layer].VlagosoderganieOfVozd = self.vozd_begin_vlaga
        # plot_data(output, layer=5, field='VlagosoderganieOfVozd', title='layer=5_VlagosoderganieOfVozd.png')

        # boundary conditions adsorption
        for dt_step in range(1, int(total_time / delta_t)):
            output[dt_step][0].TemperatureOfVozd = self.Tvoz_reg
            output[dt_step][0].VlagosoderganieOfVozd = self.vozd_begin_vlaga


        for dt_step in range(int(total_time / delta_t) - 1):
            # for layer in range(self.n - 1, 0, -1):
            # for layer in range(1, self.n):
            for layer in range(1, self.n):
                arho_voz = self.p / (287.4 * output[dt_step][layer].TemperatureOfVozd)
                aCp_voz = 10 ** 3 * (1.0005 + 1.1904 * 10 ** (-4) * (output[dt_step][layer].TemperatureOfVozd - 273))
                avyaz_voz = 0.00001717 * (output[dt_step][layer].TemperatureOfVozd / 273) ** 0.683
                avyaz_kin = avyaz_voz / arho_voz
                aLya_voz = 0.0244 * (output[dt_step][layer].TemperatureOfVozd / 273) ** 0.82

                aRe = self.v * self.d / avyaz_kin
                if isinstance(aRe, complex):
                    print(aRe, output[dt_step][layer].TemperatureOfVozd)

                if aRe < 1000:
                    aC = 0.56
                    an = 0.5
                else:
                    aC = 0.28
                    an = 0.6
                aPr = avyaz_voz * aCp_voz / aLya_voz
                aNu = aC * aRe ** an * aPr ** 0.33
                aAlf = aNu * aLya_voz / self.d
                alpha = aAlf * (self.s_pov_gran * self.n_zer * self.m_sl * self.delta_t / (self.C_ads * 700 * aCp_voz))
                c_a = self.C_ads * self.m_sl
                c_v = aCp_voz * (self.G * arho_voz * self.delta_x / self.v)

                # aPr = avyaz * aCp / aL
                # aNu = aC * aRe ** an * aPr ** 0.33
                # aAlf = aNu * aL / self.d
                #
                # alpha = aAlf * self.s_pov_gran * self.n_zer
                # c_a = self.C_ads * self.m_sl
                # # c_v = aCp * self.G * arho * self.delta_x / self.v
                # c_v = aCp * self.G

                #print(self.v * delta_t / self.delta_x * (output[dt_step][layer].TemperatureOfVozd
                #                                                   - output[dt_step][layer - 1].TemperatureOfVozd))
                # output[dt_step + 1][layer].TemperatureOfVozd = (output[dt_step][layer].TemperatureOfVozd
                #                                                 - alpha * delta_t
                #                                                 * ((output[dt_step][layer].TemperatureOfVozd
                #                                                 - output[dt_step][layer].TemperatureOfAds) / c_v)
                #                                                 - self.v * delta_t / self.delta_x #/ 13636.8
                #                                                 * (output[dt_step][layer].TemperatureOfVozd
                #                                                    - output[dt_step][layer - 1].TemperatureOfVozd))

                output[dt_step + 1][layer].TemperatureOfVozd = (output[dt_step][layer].TemperatureOfVozd
                                                                - alpha * self.C_ads
                                                                * (output[dt_step][layer].TemperatureOfVozd
                                                                    - output[dt_step][layer].TemperatureOfAds)
                                                                - self.v * delta_t / self.delta_x  # / 13636.8
                                                                * (output[dt_step][layer].TemperatureOfVozd
                                                                   - output[dt_step][layer - 1].TemperatureOfVozd))

                output[dt_step + 1][layer].TemperatureOfAds = (output[dt_step][layer].TemperatureOfAds
                                                               + alpha
                                                               * (output[dt_step][layer].TemperatureOfVozd
                                                               - output[dt_step][layer].TemperatureOfAds))




        plot_data(output, layer=20, field='TemperatureOfVozd', title='layer=20_TemperatureOfVozd_ads_v2.png')
        plot_data(output, layer=20, field='TemperatureOfAds', title='layer=20_TemperatureOfAds_ads_v2.png')

        # Plot layers
        plot_all_layers(output, n_layers=self.n, field='TemperatureOfVozd', dt_list=[30, 200, 1000],
                        title='dt=[30, 200, 1000]_TemperatureOfVozd_ads_v2.png')
        plot_all_layers(output, n_layers=self.n, field='TemperatureOfAds', dt_list=[30, 200, 1000],
                        title='dt=[30, 200, 1000]_TemperatureOfAds_ads_v2.png')
        pass
