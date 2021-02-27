from table import Table
from math import pi, ceil, exp
from constants import c_wod, k_ads, A, D_dif, dmax, total_time, delta_t
from visualization import plot_data, plot_all_layers, plot_all_layers_no_zero_layer


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
        output = [[Table() for layer in range(self.n + 2)] for dt_step in range(int(total_time / delta_t))]

        # Adsorption
        # boundary conditions
        for dt_step in range(int(total_time / delta_t)):
            output[dt_step][0].TemperatureOfVozd = self.t1
            output[dt_step][0].VlagosoderganieOfVozd = self.vozd_begin_vlaga

        # start conditions
        for layer in range(self.n):
            output[0][layer].TemperatureOfAds = self.T0
            output[0][layer].VlagosoderganieOfAds = self.ads_begin_vlaga

        for layer in range(1, self.n + 1):
            output[0][layer].TemperatureOfVozd = self.T0
            output[0][layer].VlagosoderganieOfVozd = self.ads_begin_vlaga

        for dt_step in range(int(total_time / delta_t) - 1):
            for layer in range(self.n):
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

                #todo: Fix this hardcode
                # if dt_step == int(total_time / delta_t) - 2:
                #     output[dt_step + 1][layer + 1].TemperatureOfVozd = output[dt_step + 1][layer].TemperatureOfVozd - alpha \
                #                                                    * (output[dt_step + 1][layer].TemperatureOfVozd
                #                                                       - output[dt_step + 1][layer].TemperatureOfAds) \
                #                                                    / (self.G * aCp_voz)

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


        # Desorption

        # boundary conditions
        for dt_step in range(int(total_time / delta_t)):
            output[dt_step][self.n - 1].TempOfVozdReg = self.Tvoz_reg
            output[dt_step][self.n - 1].VlagOfVozdReg = self.C_reg


        # start conditions
        for layer in range(self.n - 1, 0, -1):
            output[0][layer].TempOfAdsReg = self.t1
            output[0][layer].VlagOfAdsReg = output[int(total_time / delta_t - 2)][layer].VlagosoderganieOfAds


        for layer in range(self.n - 2, -1, -1):
            output[0][layer].TempOfVozdReg = self.t1
            output[0][layer].VlagOfVozdReg = output[int(total_time / delta_t - 2)][layer + 1].VlagosoderganieOfVozd

        for dt_step in range(int(total_time / delta_t - 1)):
            for layer in range(self.n - 1, 0, -1):
                arho_voz = self.p / (287.4 * output[dt_step][layer].TempOfVozdReg)
                avyaz_voz = 0.00001717 * (output[dt_step][layer].TempOfVozdReg / 273) ** 0.683
                aLya_voz = 0.0244 * (output[dt_step][layer].TempOfVozdReg / 273) ** 0.82
                aRe = arho_voz * self.v_reg * self.d / avyaz_voz
                aCp_voz = 1000
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

                output[dt_step][layer - 1].TempOfVozdReg = output[dt_step][layer].TempOfVozdReg - (alpha
                            * (output[dt_step][layer].TempOfVozdReg - output[dt_step][layer].TempOfAdsReg) / (self.G_reg * aCp_voz))
                aQ = self.G_reg * aCp_voz * (output[dt_step][layer].TempOfVozdReg - output[dt_step][layer - 1].TempOfVozdReg) \
                            * delta_t
                output[dt_step + 1][layer].TempOfAdsReg = aQ / (self.C_ads * self.m_sl) + output[dt_step][layer].TempOfAdsReg
                aB_mass_pov = aAlf / aCp_voz
                Max_vlag_vozd = 0.00021 * (output[dt_step][layer].TempOfVozdReg ** 3) - 0.016 \
                            * (output[dt_step][layer].TempOfVozdReg ** 2) + 0.176 * (output[dt_step][layer].TempOfVozdReg) + 0.0054

                output[dt_step][layer - 1].VlagOfVozdReg = output[dt_step][layer].VlagOfVozdReg - ((aB_mass_pov * (
                        output[dt_step][layer].VlagOfVozdReg - output[dt_step][layer].VlagOfAdsReg) * self.s_pov_gran * self.n_zer) / self.G_reg)
                output[dt_step + 1][layer].VlagOfAdsReg = output[dt_step][layer].VlagOfAdsReg + (
                        output[dt_step][layer].VlagOfVozdReg - output[dt_step][layer - 1].VlagOfVozdReg)

                if output[dt_step][layer - 1].VlagOfVozdReg > Max_vlag_vozd:
                    output[dt_step + 1][layer].VlagOfAdsReg = output[dt_step + 1][layer].VlagOfAdsReg + \
                                                              (output[dt_step][layer - 1].VlagOfVozdReg - Max_vlag_vozd)
                    output[dt_step][layer - 1].VlagOfVozdReg = Max_vlag_vozd


                output[dt_step + 1][layer].TempOfAdsReg = (aQ + ((k_ads / 1000) * (output[dt_step + 1][layer].VlagOfAdsReg
                        - output[dt_step][layer].VlagOfAdsReg)) + ((c_wod / 1000) * (output[dt_step][layer].TempOfVozdReg
                        - output[dt_step][layer].TempOfAdsReg) * (output[dt_step + 1][layer].VlagOfAdsReg - output[dt_step][layer].VlagOfAdsReg))) \
                        / ((self.C_ads * self.m_sl) + (output[dt_step + 1][layer].VlagOfAdsReg * (c_wod / 1000))) + output[dt_step][layer].TempOfAdsReg



        # plot_data(output, layer=999, field='TemperatureOfVozd', title='layer=999_TemperatureOfVozd_ads_v2.png')
        # plot_data(output, layer=999, field='TemperatureOfAds', title='layer=999_TemperatureOfAds_ads_v2.png')
        #
        # plot_data(output, layer=0, field='TemperatureOfVozd', title='layer=0_TemperatureOfVozd_ads_v2.png')
        # plot_data(output, layer=0, field='TemperatureOfAds', title='layer=0_TemperatureOfAds_ads_v2.png')
        #
        # plot_data(output, layer=1, field='TemperatureOfVozd', title='layer=1_TemperatureOfVozd_ads_v2.png')
        # plot_data(output, layer=1, field='TemperatureOfAds', title='layer=1_TemperatureOfAds_ads_v2.png')
        #
        # plot_data(output, layer=100, field='TemperatureOfVozd', title='layer=100_TemperatureOfVozd_ads_v2.png')
        # plot_data(output, layer=100, field='TemperatureOfAds', title='layer=100_TemperatureOfAds_ads_v2.png')


        ################################################################################################################
        ################################################################################################################
        ################################################################################################################
        # ADSORBTION
        # Plot layers
        plot_all_layers(output, n_layers=self.n, field='TemperatureOfVozd', dt_list=[30, 120, 300, 3600, 7200, 10800, 14398],
                        title='KOSTYA_dt=[30с, 2 мин, 5 мин, 1 час, 2 часа, 3 часаб 4 ч]_TemperatureOfVozd_ads_v2.png')
        plot_all_layers(output, n_layers=self.n, field='TemperatureOfAds', dt_list=[30, 120, 300, 3600, 7200, 10800, 14398],
                        title='KOSTYA_dt=[30с, 2 мин, 5 мин, 1 час, 2 часа, 3 часа, 4 ч]_TemperatureOfAds_ads_v2.png')

        # Plot layers
        plot_all_layers(output, n_layers=self.n, field='VlagosoderganieOfVozd', dt_list=[300, 3600, 7200, 10800, 12600, 14390, 14398],
                         title='KOSTYA_dt=[5 мин, 1 час, 2 часа, 3 часа, 3.5 ч, 4 ч]_VlagosoderganieOfVozd_ads.png')
        plot_all_layers(output, n_layers=self.n, field='VlagosoderganieOfAds', dt_list=[300, 3600, 7200, 10800, 12600, 14390, 14398],
                        title='KOSTYA_dt=[5 мин, 1 час, 2 часа, 3 часа, 3.5 ч, 4 ч]_VlagosoderganieOfAds_ads.png')

        ################################################################################################################
        ################################################################################################################
        ################################################################################################################
        # DESORBTION

        # # plot_data(output, layer=20, field='TempOfVozdReg', title='layer=20_TempOfVozd_Reg.png')
        # # plot_data(output, layer=20, field='TempOfAdsReg', title='layer=20_TempOfAds_Reg.png')
        #
        # Plot layers
        plot_all_layers_no_zero_layer(output, n_layers=self.n, field='TempOfVozdReg', dt_list=[300, 600, 900, 1800, 3600, 7200, 10800, 12600, 14390, 14398],
                        title='KOSTYA_dt=[5 мин, 10 мин, 15 мин, 1 час, 1.5 часа, 2 часа, 3 часа, 3.5 ч, 4 ч]_TempOfVozd_Reg.png')
        plot_all_layers_no_zero_layer(output, n_layers=self.n, field='TempOfAdsReg', dt_list=[300, 600, 900, 1800, 3600, 7200, 10800, 12600, 14390, 14398],
                        title='KOSTYA_dt=[5 мин, 10 мин, 15 мин, 1 час, 1.5 часа, 2 часа, 3 часа, 3.5 ч, 4 ч]_TempOfAds_Reg.png')

        # Plot layers
        plot_all_layers_no_zero_layer(output, n_layers=self.n, field='VlagOfVozdReg', dt_list=[300, 600, 900, 1800, 3600, 7200, 10800, 12600, 14390, 14398],
                        title='KOSTYA_dt=[5 мин, 10 мин, 15 мин, 1 час, 1.5 часа, 2 часа, 3 часа, 3.5 ч, 4 ч]_VlagOfVozdReg_Reg.png')
        plot_all_layers_no_zero_layer(output, n_layers=self.n, field='VlagOfAdsReg', dt_list=[300, 600, 900, 1800, 3600, 7200, 10800, 12600, 14390, 14398],
                        title='KOSTYA_dt=[5 мин, 10 мин, 15 мин, 1 час, 1.5 часа, 2 часа, 3 часа, 3.5 ч, 4 ч]_VlagOfAdsReg_Reg.png')

        pass
