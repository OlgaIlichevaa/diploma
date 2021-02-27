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

    def go(self):
        output = [[Table() for layer in range(self.n)] for dt_step in range(int(total_time / delta_t))]

        # Adsorption

        # start conditions
        for layer in range(self.n):
            output[0][layer].TemperatureOfAds = self.T0
            output[0][layer].VlagosoderganieOfAds = self.ads_begin_vlaga
            output[0][layer].TemperatureOfVozd = self.T0
            output[0][layer].VlagosoderganieOfVozd = self.ads_begin_vlaga
        # boundary conditions
        for dt_step in range(1, int(total_time / delta_t)):
            output[dt_step][0].TemperatureOfVozd = self.T0
            output[dt_step][0].VlagosoderganieOfVozd = self.vozd_begin_vlaga


        for dt_step in range(int(total_time / delta_t) - 1):
            for layer in range(self.n - 1):
                arho_voz = self.p / (287.4 * self.T0)
                avyaz_voz = 0.00001717 * (self.T0 / 273) ** 0.683
                aLya_voz = 0.0244 * ((self.T0 / 273) ** 0.82)
                aCp_voz = 1000
                aRe = arho_voz * self.v * self.d / avyaz_voz
                if aRe < 1000:
                    aC = 0.56
                    an = 0.5
                else:
                    aC = 0.28
                    an = 0.6
                aPr = avyaz_voz * aCp_voz / aLya_voz
                aNu = aC * aRe ** an * aPr ** 0.33
                aAlf = aNu * aLya_voz / self.d

                alpha = aAlf * self.s_pov_gran * self.n_zer


                output[dt_step][layer + 1].TemperatureOfVozd = output[dt_step][layer].TemperatureOfVozd - alpha \
                                                           * (output[dt_step][layer].TemperatureOfVozd
                                                              - output[dt_step][layer].TemperatureOfAds)\
                                                               / (self.G * aCp_voz)

                aQ = self.G * aCp_voz * (output[dt_step][layer].TemperatureOfVozd
                                         - output[dt_step][layer + 1].TemperatureOfVozd) * delta_t

                output[dt_step + 1][layer].TemperatureOfAds = aQ / (self.C_ads * self.m_sl) \
                                                              + output[dt_step][layer].TemperatureOfAds

                b = 0.5 * exp((k_ads) / (self.R_vozd * output[dt_step + 1][layer].TemperatureOfAds))
                X = (A * b * (output[dt_step][layer].VlagosoderganieOfVozd)) \
                    / (1 + b * (output[dt_step][layer].VlagosoderganieOfVozd))
                aB_mass_pov = aAlf / aCp_voz
                X_max_vlag_sloy = self.m_sl * X * 1000
                k = X_max_vlag_sloy / output[dt_step][layer].VlagosoderganieOfVozd
                d_ads = output[dt_step][layer].VlagosoderganieOfAds / k
                output[dt_step][layer + 1].VlagosoderganieOfVozd = output[dt_step][layer].VlagosoderganieOfVozd - (
                        (aB_mass_pov * (output[dt_step][layer].VlagosoderganieOfVozd - d_ads) * self.s_pov_gran * self.n_zer) / self.G)
                output[dt_step + 1][layer].VlagosoderganieOfAds = output[dt_step][layer].VlagosoderganieOfAds + (
                        output[dt_step][layer].VlagosoderganieOfVozd - output[dt_step][layer + 1].VlagosoderganieOfVozd)
                if (output[dt_step + 1][layer].VlagosoderganieOfAds > X_max_vlag_sloy):
                    output[dt_step][layer + 1].VlagosoderganieOfVozd = output[dt_step][layer + 1].VlagosoderganieOfVozd + \
                                                                   (output[dt_step + 1][layer].VlagosoderganieOfAds-X_max_vlag_sloy)
                    output[dt_step + 1][layer].VlagosoderganieOfAds = X_max_vlag_sloy

                output[dt_step + 1][layer].TemperatureOfAds = (aQ + k_ads * (output[dt_step + 1][layer].VlagosoderganieOfAds - output[dt_step][layer].VlagosoderganieOfAds) + (
                                                              (c_wod / 1000) * (output[dt_step][layer].TemperatureOfVozd - output[dt_step][layer].TemperatureOfAds)
                                                              * (output[dt_step + 1][layer].VlagosoderganieOfAds
                                                                 - output[dt_step][layer].VlagosoderganieOfAds))) / ((self.C_ads * self.m_sl)
                                                            + (output[dt_step + 1][layer].VlagosoderganieOfAds * (c_wod / 1000))) + output[dt_step][layer].TemperatureOfAds


        # # Desorption
        #
        # # start conditions
        # for layer in range(self.n):
        #     output[0][layer].TempOfAdsReg = self.t1
        #     output[0][layer].TempOfVozdReg = self.Tvoz_reg
        #     output[0][layer].VlagOfAdsReg = output[int(total_time / delta_t - 1)][self.n - layer].VlagosoderganieOfAds
        #     output[0][layer].VlagOfVozdReg = self.C_reg
        #
        # # boundary conditions
        # for dt_step in range(int(total_time / delta_t)):
        #     output[dt_step][0].TempOfVozdReg = self.Tvoz_reg
        #     output[dt_step][0].VlagOfVozdReg = output[int(total_time / delta_t - 1)][self.n - 1].VlagosoderganieOfVozd
        #
        #
        # for dt_step in range(1, int(total_time / delta_t)):
        #     for layer in range(1, self.n):
        #         arho_voz = self.p / (287.4 * output[dt_step - 1][layer - 1].TempOfVozdReg)
        #         avyaz_voz = 0.00001717 * (output[dt_step - 1][layer - 1].TempOfVozdReg / 273) ** 0.683
        #         aLya_voz = 0.0244 * (output[dt_step - 1][layer - 1].TempOfVozdReg / 273) ** 0.82
        #         aRe = arho_voz * self.v_reg * self.d / avyaz_voz
        #         aCp_voz = 1000
        #         if aRe < 1000:
        #             aC = 0.56
        #             an = 0.5
        #         else:
        #             aC = 0.28
        #             an = 0.6
        #         aPr = avyaz_voz * aCp_voz / aLya_voz
        #         aNu = aC * aRe ** an * aPr ** 0.33
        #         aAlf = aNu * aLya_voz / self.d
        #
        #
        #         output[dt_step - 1][layer].TempOfVozdReg = output[dt_step - 1][layer - 1].TempOfVozdReg - ((aAlf
        #                     * (output[dt_step - 1][layer - 1].TempOfVozdReg - output[dt_step - 1][layer - 1].TempOfAdsReg)
        #                     * self.s_pov_gran * self.n_zer) / (self.G_reg * aCp_voz))
        #         aQ = self.G_reg * aCp_voz * (output[dt_step - 1][layer - 1].TempOfVozdReg - output[dt_step - 1][layer].TempOfVozdReg) \
        #                     * delta_t
        #         output[dt_step][layer - 1].TempOfAdsReg = aQ / (self.C_ads * self.m_sl) + output[dt_step - 1][layer - 1].TempOfAdsReg
        #         aB_mass_pov = aAlf / aCp_voz
        #         Max_vlag_vozd = 0.00021 * (output[dt_step - 1][layer - 1].TempOfVozdReg ** 3) - 0.016 \
        #                     * (output[dt_step - 1][layer - 1].TempOfVozdReg ** 2) + 0.176 * (output[dt_step - 1][layer - 1].TempOfVozdReg) + 0.0054
        #
        #         output[dt_step - 1][layer].VlagOfVozdReg = output[dt_step - 1][layer - 1].VlagOfVozdReg - ((aB_mass_pov * (
        #                 output[dt_step - 1][layer - 1].VlagOfVozdReg - output[dt_step - 1][layer - 1].VlagOfAdsReg) * self.s_pov_gran * self.n_zer) / self.G_reg)
        #         output[dt_step][layer - 1].VlagOfVozdReg = output[dt_step - 1][layer - 1].VlagOfAdsReg + (
        #                 output[dt_step - 1][layer - 1].VlagOfVozdReg - output[dt_step - 1][layer].VlagOfVozdReg)
        #
        #         if (output[dt_step - 1][layer].VlagOfVozdReg > Max_vlag_vozd):
        #             output[dt_step][layer - 1].VlagOfAdsReg = output[dt_step][layer - 1].VlagOfAdsReg + (output[dt_step - 1][layer].VlagOfVozdReg
        #                                                                                   - Max_vlag_vozd)
        #             output[dt_step - 1][layer].VlagOfVozdReg = Max_vlag_vozd
        #
        #
        #         output[dt_step][layer - 1].TempOfAdsReg = (aQ + ((k_ads / 1000) * (output[dt_step][layer - 1].VlagOfAdsReg
        #                 - output[dt_step - 1][layer - 1].VlagOfAdsReg)) + ((c_wod / 1000) * (output[dt_step - 1][layer - 1].TempOfVozdReg
        #                 - output[dt_step - 1][layer - 1].TempOfAdsReg) * (output[dt_step][layer - 1].VlagOfAdsReg - output[dt_step - 1][layer - 1].VlagOfAdsReg))) \
        #                 / ((self.C_ads * self.m_sl) + (output[dt_step][layer - 1].VlagOfAdsReg * (c_wod / 1000))) + output[dt_step - 1][layer - 1].TempOfAdsReg



        plot_data(output, layer=999, field='TemperatureOfVozd', title='layer=999_TemperatureOfVozd_ads_v2.png')
        plot_data(output, layer=999, field='TemperatureOfAds', title='layer=999_TemperatureOfAds_ads_v2.png')

        plot_data(output, layer=0, field='TemperatureOfVozd', title='layer=0_TemperatureOfVozd_ads_v2.png')
        plot_data(output, layer=0, field='TemperatureOfAds', title='layer=0_TemperatureOfAds_ads_v2.png')

        plot_data(output, layer=1, field='TemperatureOfVozd', title='layer=1_TemperatureOfVozd_ads_v2.png')
        plot_data(output, layer=1, field='TemperatureOfAds', title='layer=1_TemperatureOfAds_ads_v2.png')

        plot_data(output, layer=100, field='TemperatureOfVozd', title='layer=100_TemperatureOfVozd_ads_v2.png')
        plot_data(output, layer=100, field='TemperatureOfAds', title='layer=100_TemperatureOfAds_ads_v2.png')




        # Plot layers
        plot_all_layers(output, n_layers=self.n, field='TemperatureOfVozd', dt_list=[30, 120, 300, 3600, 7200, 10800, 14399],
                        title='KOSTYA_dt=[30с, 2 мин, 5 мин, 1 час, 2 часа, 3 часаб 4 ч]_TemperatureOfVozd_ads_v2.png')
        plot_all_layers(output, n_layers=self.n, field='TemperatureOfAds', dt_list=[30, 120, 300, 3600, 7200, 10800, 14399],
                        title='KOSTYA_dt=[30с, 2 мин, 5 мин, 1 час, 2 часа, 3 часа, 4 ч]_TemperatureOfAds_ads_v2.png')
        #
        #
        #
        #
        # # plot_data(output, layer=20, field='TempOfVozdReg', title='layer=20_TempOfVozd_Reg.png')
        # # plot_data(output, layer=20, field='TempOfAdsReg', title='layer=20_TempOfAds_Reg.png')
        #
        # # Plot layers
        # plot_all_layers(output, n_layers=self.n, field='TempOfVozdReg', dt_list=[30, 1200, 2400, 3600, 5400],
        #                 title='dt=[30с, 20 мин, 40 мин, 1 ч, 1,5 ч]_TempOfVozd_Reg.png')
        # plot_all_layers(output, n_layers=self.n, field='TempOfAdsReg', dt_list=[30, 60, 180, 300],
        #                 title='dt=[30с, 2 мин, 3 мин, 5 мин]_TempOfAds_Reg.png')
        #
        #
        #
        #
        #
        #
        #
        #
        plot_data(output, layer=999, field='VlagosoderganieOfVozd', title='layer=999_VlagosoderganieOfVozd_ads.png')
        plot_data(output, layer=999, field='VlagosoderganieOfAds', title='layer=999_VlagosoderganieOfAds_ads.png')


        plot_data(output, layer=0, field='VlagosoderganieOfVozd', title='layer=0_VlagosoderganieOfVozd_ads.png')
        plot_data(output, layer=0, field='VlagosoderganieOfAds', title='layer=0_VlagosoderganieOfAds_ads.png')

        # Plot layers
        plot_all_layers(output, n_layers=self.n, field='VlagosoderganieOfVozd', dt_list=[30, 120, 300, 3600, 7200, 10800, 14399],
                        title='dt=[30с, 2 мин, 5 мин, 1 час, 2 часа, 3 часаб 4 ч]_VlagosoderganieOfVozd_ads.png')
        plot_all_layers(output, n_layers=self.n, field='VlagosoderganieOfAds', dt_list=[30, 120, 300, 3600, 7200, 10800, 14399],
                        title='dt=[30с, 2 мин, 5 мин, 1 час, 2 часа, 3 часа, 4 ч]_VlagosoderganieOfAds_ads.png')
        #
        #
        # # plot_data(output, layer=20, field='VlagOfVozdReg', title='layer=20_VlagOfVozdReg_Reg.png')
        # # plot_data(output, layer=20, field='VlagOfAdsReg', title='layer=20_VlagOfAdsReg_Reg.png')
        #
        # # Plot layers
        # plot_all_layers(output, n_layers=self.n, field='VlagOfVozdReg', dt_list=[30, 3600, 5400, 10800, 12600],
        #                 title='dt=[30с, 1 ч, 1.5 ч, 3 ч, 3.5 ч]_VlagOfVozdReg_Reg.png')
        # plot_all_layers(output, n_layers=self.n, field='VlagOfAdsReg', dt_list=[30, 900, 3600, 5400, 9000, 10800],
        #                 title='dt=[30с, 15 мин, 1 ч, 1.5 ч, 2.5 ч, 3 ч]_VlagOfAdsReg_Reg.png')
        #
        #
        #
        # pass
