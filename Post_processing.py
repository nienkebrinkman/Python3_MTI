# ---------------------------------------------------------------------------------------------------------------------#
#                                                 Post Processing                                                      #
# ---------------------------------------------------------------------------------------------------------------------#
# from Plots import Plots
import os
import itertools
import glob
import pandas as pd
import matplotlib.pyplot as plt
import yaml
from obspy.imaging.beachball import aux_plane
from pandas.plotting import autocorrelation_plot
from obspy.imaging.beachball import beachball
from pandas.plotting import scatter_matrix
import pylab
import obspy

# from mpl_toolkits.basemap import Basemap
import geographiclib.geodesic as geo
# import mplstereonet
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from obspy.imaging.beachball import beach
import matplotlib.image as mpimg
import io

from Get_Parameters import Get_Parameters

from pyrocko import moment_tensor as mtm


def main():


    # NPZ_path = '/home/nienke/Documents/Master/Data/Database/TAYAK.npz'
    # file = np.load(NPZ_path)

    ## Post - Processing [processing the results from inversion]
    result = Post_processing_sdr()

    strike, dip, rake = aux_plane(238, 80, 143) # check quickly any moment tensor

    directory = '/home/nienke/Documents/Master/Data/Mars/S0235b/waveforms/Output/'
    # directory = '/home/nienke/Documents/Master/Data/Mars/S0173a/waveforms/Output/'
    # path_to_file = directory + 'TAYAK_Shift_1.txt'
    # path_to_file = directory + 'DWTHot_Amp_Shift.txt'
    path_to_file = directory + 'EH45Tcold_Amp_Shift.txt'

    savename = 'Trials'
    show = False  # Choose True for direct show, choose False for saving
    skiprows = 56#48
    # column_names = ["Epi", "Depth", "Strike", "Dip", "Rake", "M0", "Total_misfit", "p_z", "p_r", "s_z", "s_r", "s_t",
    #                 'Xi_amp', 'Shift_S', 'Shift_P', 'accept']

    column_names = ["Epi", "Depth", "Strike", "Dip", "Rake", "M0", "Total_misfit", "p_z", "p_r", "s_z", "s_r",
                         "s_t", 'Npz', 'Npr', 'Nsz', 'Nsr', 'Nst', 'PZ_Amplitude', 'Shift_S', 'Shift_P', 'accept']
    burnin = 0

    File = np.loadtxt(path_to_file, delimiter=',', skiprows=skiprows)
    strike, dip, rake = aux_plane(np.mean(File[:,2]), np.mean(File[:,3]), np.mean(File[:,4]))

    # strike, dip, rake = aux_plane(183, 55, -146)
    # real_v = np.array([None,None, strike, dip, rake, None])
    # real_v = np.array([88.4756, 38438, s2, d2, r2, exp])  # MSS event 5.0
    real_v = None
    #
    # result.Scatter_3d(filepath=path_to_file, savename=savename, directory=directory, skiprows=skiprows,
    #                column_names=column_names,  burnin=burnin,real_v=real_v)
    result.trace(filepath=path_to_file, savename=savename, directory=directory, skiprows=skiprows,
                 column_names=column_names,burnin=burnin ,lowest = 'Total_misfit',real_v=real_v)
    # result.get_BBB(filepath=path_to_file, savename=savename, directory=directory, skiprows=skiprows,
    #                column_names=column_names,  burnin=burnin,lowest = 'Total_misfit',real_v=real_v)
    # result.Cluster(filepath=path_to_file, savename=savename, directory=directory, skiprows=skiprows,
    #                column_names=column_names,  burnin=burnin,real_v=real_v)

    # result.marginal_grid( savename=savename, directory=directory,samples=File[:,:], dimensions_list = [6,7,8,9,10,11], show = False)
    # result.get_convergence(filepath=path_to_file, savename = savename, directory = directory, skiprows = skiprows, burnin=burnin,column_names = column_names, show=False)
    # result.Get_cross_plots(filepath=path_to_file, savename = savename, directory = directory, skiprows = skiprows, burnin=burnin,column_names = column_names, show=False)

    # result.event_plot(savename = savename,directory = directory,la_receiver = 4.502384, lo_receiver = 135.623447 , la_source = 3.45, lo_source = 164.68)


class Post_processing_sdr:
    def Scatter_3d(self, filepath, savename, directory, skiprows, column_names, burnin, real_v):
        dir = directory + '/%s' % (savename.strip('.yaml'))
        if not os.path.exists(dir):
            os.makedirs(dir)

        data = np.loadtxt(filepath, delimiter=',', skiprows=skiprows)

        df = pd.DataFrame(data, columns=column_names)
        n_lowest = 10000
        lowest_indices = df['Total_misfit'].values.argsort()[0:n_lowest]

        ## TAYAK
        # # PZ = (df['p_z'].values)
        # # PR = (df['p_r'].values) * 0.017
        # # SZ = (df['s_z'].values) * 0.019
        # # SR = (df['s_r'].values) * 0.0556
        # # ST = (df['s_t'].values) * 2.33
        # # AMP = (( (np.abs(np.log10(df['Npz'].values ) - np.log10(df['Nst'].values ))) / np.log(2)) ) ** (1/4) / 2
        # # Misfit = PZ + ST + PR + SZ + SR + AMP
        #
        # # indices = np.where(Misfit < 0.13)
        # # indices = np.where(Misfit < 0.27)
        # # indices = np.where(Misfit < 1.25)
        # # indices = np.where(Misfit < 1.55)
        # ## DWTHot_Shift_1
        # # indices = np.where(Misfit < 0.184) # PZ
        # # indices = np.where(Misfit < 0.00000002) # AMP --> Low treshold
        # # indices = np.where(Misfit < 0.000002) # AMP --> High treshold
        # # indices = np.where(Misfit < 0.42) # PZ + ST
        #
        # misfits = Misfit[indices]
        #
        # strike = df['Strike'].values[indices]
        # dip = df['Dip'].values[indices]
        # rake = df['Rake'].values[indices]

        ## EVENT_235_2
        ####PZ
        # G1_S = np.hstack((strike[(strike < 70) ] , strike[(strike > 349) ] ))
        # G1_D = np.hstack(( dip[(strike < 70) ] ,  dip[(strike > 349) ] ))
        # G1_R = np.hstack(( rake[(strike < 70) ] , rake[(strike > 349) ] ))
        # M1 = np.hstack(( misfits[(strike < 70)] , misfits[(strike > 349) ] ))
        #
        # G21_S = strike[( ( (70 < strike) & (strike < 250)) & ((-180 < rake ) & (rake <-100)) )]
        # G21_D = dip[( ( (70 < strike) & (strike < 250)) & ((-180 < rake ) & (rake <-100)) )]
        # G21_R = rake[( ( (70 < strike) & (strike < 250)) & ((-180 < rake ) & (rake <-100)) )]
        # M21 = misfits[( ( (70 < strike) & (strike < 250)) & ((-180 < rake ) & (rake <-100)) )]
        #
        # G22_S = strike[( ( (70 < strike) & (strike < 250)) & ((-100 < rake ) & (rake <0)) )]
        # G22_D = dip[( ( (70 < strike) & (strike < 250)) & ((-100 < rake ) & (rake <0)) )]
        # G22_R = rake[( ( (70 < strike) & (strike < 250)) & ((-100 < rake ) & (rake <0)) )]
        # M22 = misfits[( ( (70 < strike) & (strike < 250)) & ((-100 < rake ) & (rake <0)) )]
        #
        # G3_S = strike[(250< strike) & (strike < 349) ]
        # G3_D = dip[(250< strike) & (strike < 349) ]
        # G3_R = rake[(250< strike) & (strike < 349) ]
        # M3 = misfits[(250< strike) & (strike < 349) ]

        ####PZ+ST
        # G1_S = strike[(strike < 250) ]
        # G1_D = dip[(strike < 250) ]
        # G1_R = rake[(strike < 250) ]
        # M1 = misfits[(strike < 250) ]
        #
        # G2_S = strike[(strike > 250) ]
        # G2_D = dip[(strike > 250) ]
        # G2_R = rake[(strike > 250) ]
        # M2 = misfits[(strike > 250) ]

        ####TOTAL
        # G1_S = np.hstack(( strike[(strike < 100) ], strike[(strike>255) ] ))
        # G1_D = np.hstack((dip[(strike < 100)  ] ,dip[(strike>255)]))
        # G1_R = np.hstack(( rake[(strike < 100) ] , rake[(strike>255)] ))
        # M1 = np.hstack(( misfits[(strike < 100) ], misfits[(strike>255)] ))
        #
        # G2_S = strike[((strike > 100)  & (strike<255) )]
        # G2_D = dip[((strike > 100)  & (strike<255) )]
        # G2_R = rake[((strike > 100)  & (strike<255) )]
        # M2 = misfits[((strike > 100)  & (strike<255) )]

        ## DWTHot
        # PZ = (df['p_z'].values)
        # PR = (df['p_r'].values) * 0.1
        # SZ = (df['s_z'].values) * 0.055
        # SR = (df['s_r'].values) * 0.21
        # ST = (df['s_t'].values) * 0.2
        # AMP =(( (np.abs(np.log10(df['Npz'].values ) - np.log10(df['Nst'].values ))) / np.log(2)) ) ** (1/4)
        #
        # Misfit = PZ + ST + PR + SZ + SR #+ AMP
        #
        # indices = np.where(Misfit < 0.13)
        # indices = np.where(Misfit < 5.95)
        #
        #
        # misfits = Misfit[indices]
        #
        # strike = df['Strike'].values[indices]
        # dip = df['Dip'].values[indices]
        # rake = df['Rake'].values[indices]

        ###PZ
        # G1_S = strike[(rake < -100) ]
        # G1_D =  dip[(rake < -100) ]
        # G1_R = rake[(rake < -100) ]
        # M1 =  misfits[(rake < -100) ]
        #
        # G2_S =  strike[((rake > -100) & (rake < 0))]
        # G2_D =  dip[((rake > -100) & (rake < 0))]
        # G2_R =  rake[((rake > -100) & (rake < 0))]
        # M2   =  misfits[((rake > -100) & (rake < 0))]
        #
        # G3_S =  strike[((rake > 0) & (rake < 100))]
        # G3_D =  dip[((rake > 0) & (rake < 100))]
        # G3_R =  rake[((rake > 0) & (rake < 100))]
        # M3   =  misfits[((rake > 0) & (rake < 100))]
        #
        # G4_S = strike[(rake > 100) ]
        # G4_D = dip[(rake > 100) ]
        # G4_R = rake[(rake > 100) ]
        # M4   = misfits[(rake > 100) ]


        ### PZ + ST

        # G1_S = strike[(rake < -100) ]
        # G1_D =  dip[(rake < -100) ]
        # G1_R = rake[(rake < -100) ]
        # M1 =  misfits[(rake < -100) ]
        #
        # G2_S =  strike[((rake > -100) & (rake < 0))]
        # G2_D =  dip[((rake > -100) & (rake < 0))]
        # G2_R =  rake[((rake > -100) & (rake < 0))]
        # M2   =  misfits[((rake > -100) & (rake < 0))]
        #
        # G3_S =  strike[((rake > 0) & (rake < 100))]
        # G3_D =  dip[((rake > 0) & (rake < 100))]
        # G3_R =  rake[((rake > 0) & (rake < 100))]
        # M3   =  misfits[((rake > 0) & (rake < 100))]
        #
        # G4_S = strike[(rake > 100) ]
        # G4_D = dip[(rake > 100) ]
        # G4_R = rake[(rake > 100) ]
        # M4   = misfits[(rake > 100) ]


        #### TOTAL (INCLUDING AMPLITUDE)

        # G1_S = strike[(strike > 225) ]
        # G1_D =  dip[(strike > 225) ]
        # G1_R = rake[(strike > 225) ]
        # M1 =  misfits[(strike > 225) ]
        #
        # G2_S =  strike[(strike < 225) ]
        # G2_D =  dip[(strike < 225) ]
        # G2_R =  rake[(strike < 225) ]
        # M2   =  misfits[(strike < 225) ]


        ## EH34Tcold
        PZ = (df['p_z'].values)
        PR = (df['p_r'].values) * 0.5
        SZ = (df['s_z'].values) * 0.16
        SR = (df['s_r'].values) * 0.15
        ST = (df['s_t'].values) * 0.94
        AMP =(( (np.abs(np.log10(df['Npz'].values ) - np.log10(df['Nst'].values ))) / np.log(2)) ) ** (1/4) * 1.3

        Misfit = PZ + ST + PR + SZ + SR #+ AMP

        indices = np.where(Misfit < 0.5)
        indices = np.where(Misfit < 0.7)
        indices = np.where(Misfit < 8.75)


        misfits = Misfit[indices]

        strike = df['Strike'].values[indices]
        dip = df['Dip'].values[indices]
        rake = df['Rake'].values[indices]

        ## PZ
        # G1_S =  strike[(strike < 160) ]
        # G1_D =  dip[(strike < 160) ]
        # G1_R =  rake[(strike < 160) ]
        # M1   =  misfits[(strike < 160) ]
        #
        # G2_S =  strike[(strike > 160) ]
        # G2_D =  dip[(strike > 160) ]
        # G2_R =  rake[(strike > 160) ]
        # M2   =  misfits[(strike > 160) ]

        ## PZ + ST

        # G1_S =  strike[(strike < 200) ]
        # G1_D =  dip[(strike < 200) ]
        # G1_R =  rake[(strike < 200) ]
        # M1   =  misfits[(strike < 200) ]
        #
        # G2_S =  strike[(strike > 200) ]
        # G2_D =  dip[(strike > 200) ]
        # G2_R =  rake[(strike > 200) ]
        # M2   =  misfits[(strike > 200) ]

        ## TOTAL

        G1_S =  strike[(strike < 125) ]
        G1_D =  dip[(strike < 125) ]
        G1_R =  rake[(strike < 125) ]
        M1   =  misfits[(strike < 125) ]

        G2_S =  strike[(strike > 125) ]
        G2_D =  dip[(strike > 125) ]
        G2_R =  rake[(strike > 125) ]
        M2   =  misfits[(strike > 125) ]

        # fig = plt.figure(figsize=(5,2.5))
        # ax = plt.subplot(111)
        # plt.scatter(G4_R,M4, c = 'red', label = 'Group 4')
        # plt.ylabel('Misfit')
        # plt.xlabel('Rake')
        # plt.xlim(-180,180)
        # plt.legend()
        # plt.tight_layout()
        # plt.savefig(dir + '/Xi.pdf')
        #
        # plt.close()


        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure(figsize=(25,10))
        ax = plt.subplot(111, projection = '3d')
        # Scat = ax.scatter(strike, dip, rake, c=misfits,cmap = 'rainbow')
        Scat = ax.scatter(G1_S, G1_D, G1_R,label = 'Group 1')
        Scat = ax.scatter(G2_S, G2_D, G2_R,label = 'Group 2')
        # Scat = ax.scatter(G3_S, G3_D, G3_R,label = 'Group 3')
        # Scat = ax.scatter(G4_S, G4_D, G4_R, label='Group 4')
        ax.set_xlabel('Strike')
        ax.set_xlim(0,360)
        ax.set_ylabel('Dip')
        ax.set_ylim(0,90)
        ax.set_zlabel('Rake')
        ax.set_zlim(-180,180)
        plt.legend(loc='upper left')
        # cbar = plt.colorbar()
        # cbar.set_label('Misfit', rotation=270, labelpad=20)
        fig.colorbar(Scat)
        plt.show()

    def get_convergence(self, filepath, savename, directory, skiprows, burnin,column_names, show=True):
        dir = directory + '/%s' % (savename)
        if not os.path.exists(dir):
            os.makedirs(dir)

        if filepath.endswith('.yaml') == True:
            with open(filepath, 'r') as stream:
                data_file = yaml.load(stream)
                stream.close()
                data = data_file['data']
                # parameters = data_file['parameters']
        else:
            # parameters = open(filepath, "r").readlines()[:33]
            data = np.loadtxt(filepath, delimiter=',', skiprows=skiprows)

        df = pd.DataFrame(data,
                          columns=column_names)
        # plt.figure(1)
        # ax = plt.subplot(111)
        # ax.plot(np.arange(0, len(df['Total_misfit'])), df['Total_misfit'])
        # # ax.plot(np.arange(0, len(df['Total_misfit'])), df['p_z'],c='r')
        # # ax.plot(np.arange(0, len(df['Total_misfit'])), df['s_t'],c='g')
        # plt.yscale('log')
        # plt.xlabel('Iteration')
        # plt.ylabel('-Log(likelihood)')
        # # ax.invert_yaxis()
        # ax.xaxis.tick_top()
        #
        # plt.tight_layout()
        # if show == True:
        #     plt.show()
        # else:
        #     plt.savefig(dir + '/Convergence.pdf')
        #     plt.close()

        #### PLOT:
        # amount_of_samples = len(df['Total_misfit'])
        # loop_length = int(amount_of_samples/10 -1)
        # epi_mean    = np.zeros(loop_length)
        # depth_mean  = np.zeros(loop_length)
        # strike_mean = np.zeros(loop_length)
        # dip_mean    = np.zeros(loop_length)
        # rake_mean   = np.zeros(loop_length)
        # M0_mean     = np.zeros(loop_length)
        #
        #
        # for i in range(loop_length):
        #     epi_mean[i] = np.mean(df['Epi'].values[burnin:int(i * 10 + 1)])
        #     depth_mean[i] = np.mean(df['Depth'].values[burnin:int(i * 10 + 1)])
        #     strike_mean[i] = np.mean(df['Strike'].values[burnin:int(i * 10 + 1)])
        #     dip_mean[i] = np.mean(df['Dip'].values[burnin:int(i * 10 + 1)])
        #     rake_mean[i] = np.mean(df['Rake'].values[burnin:int(i * 10 + 1)])
        #     M0_mean[i] = np.mean(df['M0'].values[burnin:int(i * 10 + 1)])
        #
        # plt.close()
        # plt.figure(figsize=(20,10))
        # ax1 = plt.subplot(231)
        # ax1.plot(np.arange(loop_length) * 10, epi_mean , label = 'Epicentral distance')
        # ax1.set_ylabel("Mean", fontsize=25)
        # ax1.set_xlabel("Sample", fontsize=25)
        # ax1.tick_params(axis='x', labelsize=20)
        # ax1.tick_params(axis='y', labelsize=20)
        # ax1.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        # plt.legend(fontsize=20)
        # plt.tight_layout()
        # # plt.yscale('log')
        #
        #
        # ax2 = plt.subplot(232)
        # ax2.plot(np.arange(loop_length) * 10, depth_mean , label = 'Depth')
        # ax2.set_ylabel("Mean", fontsize=25)
        # ax2.set_xlabel("Sample", fontsize=25)
        # ax2.tick_params(axis='x', labelsize=20)
        # ax2.tick_params(axis='y', labelsize=20)
        # ax2.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        # plt.legend(fontsize=20)
        # plt.tight_layout()
        # # plt.yscale('log')
        #
        #
        # ax3 = plt.subplot(233)
        # ax3.plot(np.arange(loop_length) * 10, M0_mean , label = 'M0')
        # ax3.set_ylabel("Mean", fontsize=25)
        # ax3.set_xlabel("Sample", fontsize=25)
        # ax3.tick_params(axis='x', labelsize=20)
        # ax3.tick_params(axis='y', labelsize=20)
        # ax3.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        # plt.legend(fontsize=20)
        # plt.tight_layout()
        # # plt.yscale('log')
        #
        # ax4 = plt.subplot(234)
        # ax4.plot(np.arange(loop_length) * 10, strike_mean , label = 'Strike')
        # ax4.set_ylabel("Mean", fontsize=25)
        # ax4.set_xlabel("Sample", fontsize=25)
        # ax4.tick_params(axis='x', labelsize=20)
        # ax4.tick_params(axis='y', labelsize=20)
        # ax4.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        # plt.legend(fontsize=20)
        # plt.tight_layout()
        # # plt.yscale('log')
        #
        # ax5 = plt.subplot(235)
        # ax5.plot(np.arange(loop_length) * 10, dip_mean , label = 'Dip')
        # ax5.set_ylabel("Mean", fontsize=25)
        # ax5.set_xlabel("Sample", fontsize=25)
        # ax5.tick_params(axis='x', labelsize=20)
        # ax5.tick_params(axis='y', labelsize=20)
        # ax5.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        # plt.legend(fontsize=20)
        # plt.tight_layout()
        # # plt.yscale('log')
        #
        # ax6 = plt.subplot(236)
        # ax6.plot(np.arange(loop_length) * 10, rake_mean , label = 'Rake')
        # ax6.set_ylabel("Mean", fontsize=25)
        # ax6.set_xlabel("Sample", fontsize=25)
        # ax6.tick_params(axis='x', labelsize=20)
        # ax6.tick_params(axis='y', labelsize=20)
        # ax6.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        # plt.legend(fontsize=20)
        # plt.tight_layout()
        # # plt.yscale('log')
        #
        #
        # if show == True:
        #
        #     plt.show()
        # else:
        #     plt.savefig(dir + '/Parameter_convergence.pdf')
        #     plt.close()

        #### PLOT DWTHot 1:
        # PZ = (df['p_z'].values) * 1.85
        # PR = (df['p_r'].values)   # * 0.005#/ 0.1 * 0.017
        # SZ = (df['s_z'].values) * 0.073  #* 0.010#/ 0.14 * 0.019
        # SR = (df['s_r'].values) * 0.073   #* 0.03# * 0.0556)#/ 0.35) * 0.0556
        # ST = (df['s_t'].values)     # * 2.33
        # AMP = df['Xi_amp'].values * 1000 *2

        #### PLOT DWTHot 2:
        # PZ = (df['p_z'].values)
        # PR = (df['p_r'].values) * 0.1
        # SZ = (df['s_z'].values) * 0.055
        # SR = (df['s_r'].values) * 0.21
        # ST = (df['s_t'].values) * 0.2
        # AMP =(( (np.abs(np.log10(df['Npz'].values ) - np.log10(df['Nst'].values ))) / np.log(2)) ) ** (1/4)

        #TAYAK:
        # PZ = (df['p_z'].values)
        # PR = (df['p_r'].values) * 0.017
        # SZ = (df['s_z'].values) * 0.019
        # SR = (df['s_r'].values) * 0.0556
        # ST = (df['s_t'].values) * 2.33
        # AMP = (( (np.abs(np.log10(df['Npz'].values ) - np.log10(df['Nst'].values ))) / np.log(2)) ) ** (1/4) / 2

        #### EVENT 173
        # PZ = (df['p_z'].values)
        # PR = (df['p_r'].values)
        # SZ = (df['s_z'].values)
        # SR = (df['s_r'].values)
        # ST = (df['s_t'].values)
        # AMP = (np.abs(np.log10(df['Npz'].values) - np.log10(df['Nst'].values)) / np.log(2))**2


        ## EH45Tcold
        PZ = (df['p_z'].values)
        PR = (df['p_r'].values)# * 0.5
        SZ = (df['s_z'].values)# * 0.16
        SR = (df['s_r'].values)# * 0.15
        ST = (df['s_t'].values)# * 0.94
        AMP =(( (np.abs(np.log10(df['Npz'].values ) - np.log10(df['Nst'].values ))) / np.log(2)) ) ** (1/4) #* 1.3


        n_lowest = 10000
        lowest_indices = (PZ+PR+SZ+SR+ST).argsort()[0:n_lowest]
        # lowest_indices = df['Total_misfit'].values.argsort()[0:n_lowest]

        PZ = PZ[PZ.argsort()[0:n_lowest]] #PZ[lowest_indices]
        PR = PR[PR.argsort()[0:n_lowest]] #PR[lowest_indices]
        SZ = SZ[SZ.argsort()[0:n_lowest]] #SZ[lowest_indices]
        SR = SR[SR.argsort()[0:n_lowest]] #SR[lowest_indices]
        ST = ST[ST.argsort()[0:n_lowest]] #ST[lowest_indices]
        AMP = AMP[AMP.argsort()[0:n_lowest]]

        plt.figure(figsize=(15,5))
        ax1 = plt.subplot(121)
        ax1.plot(PZ, label = 'PZ, Mean:%.2f Min%.2f' %(np.mean(PZ),PZ.min()))
        ax1.plot(ST, label = 'ST, Mean:%.2f Min%.2f' %(np.mean(ST),ST.min()))
        ax1.plot(SZ, label = 'SZ, Mean:%.2f Min%.2f' %(np.mean(SZ),SZ.min()))
        ax1.plot(SR, label = 'SR, Mean:%.2f Min%.2f' %(np.mean(SR),SR.min()))
        ax1.plot(PR, label = 'PR, Mean:%.2f Min%.2f' %(np.mean(PR),PR.min()))
        ax1.plot(AMP, label = 'AMP, Mean:%.2f Min%.2f' %(np.mean(AMP),AMP.min()))
        ax1.set_ylabel("Misfit", fontsize=25)
        ax1.set_xlabel("Sample", fontsize=25)
        ax1.tick_params(axis='x', labelsize=20)
        ax1.tick_params(axis='y', labelsize=20)
        ax1.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        plt.legend(fontsize=16, loc='lower right')
        plt.tight_layout()
        plt.yscale('log')
        plt.ylim((pow(20, -1), pow(10, 1)))
        ax2 = plt.subplot(122)
        ax2.hist(PZ, bins= 80, alpha=0.8,label = 'PZ, Mean:%.2f Min%.2f' %(np.mean(PZ),PZ.min()))
        ax2.hist(ST, bins= 80, alpha=0.8,label = 'ST, Mean:%.2f Min%.2f' %(np.mean(ST),ST.min()))
        ax2.hist(SZ, bins= 80, alpha=0.8,label = 'SZ, Mean:%.2f Min%.2f' %(np.mean(SZ),SZ.min()))
        ax2.hist(SR, bins= 80, alpha=0.8,label = 'SR, Mean:%.2f Min%.2f' %(np.mean(SR),SR.min()))
        ax2.hist(PR, bins= 80, alpha=0.8,label = 'PR, Mean:%.2f Min%.2f' %(np.mean(PR),PR.min()))
        ax1.hist(AMP,bins= 80, alpha=0.8, label = 'AMP, Mean:%.2f Min%.2f' %(np.mean(AMP),AMP.min()))
        ax2.set_xlabel("Misfit", fontsize=25)
        ax2.set_ylabel("Frequency", fontsize=25)
        ax2.tick_params(axis='x', labelsize=20)
        ax2.tick_params(axis='y', labelsize=20)
        ax2.ticklabel_format(style="sci", axis='x', scilimits=(-2, 2))
        plt.xscale('log')
        plt.xlim((pow(20, -1), pow(10, 1)))
        plt.ylim(0,500)
        # plt.legend(fontsize=16, loc='upper left')
        plt.tight_layout()
        # plt.show()
        plt.savefig(dir + '/misfits.pdf')

        #
        # ax2 = plt.subplot(512)
        # ax2.plot(PR, label = 'PR, Mean:%.2f Min%.2f' %(np.mean(PR),PR.min()))
        # ax2.set_ylabel("Misfit", fontsize=25)
        # ax2.tick_params(axis='x', labelsize=20)
        # ax2.tick_params(axis='y', labelsize=20)
        # ax2.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        # plt.legend(fontsize=20)
        # plt.tight_layout()
        # plt.yscale('log')
        # plt.ylim((pow(20, -1), pow(10, 1)))
        #
        #
        # ax3 = plt.subplot(513)
        # ax3.plot(SZ, label = 'SZ, Mean:%.2f Min%.2f' %(np.mean(SZ),SZ.min()))
        # ax3.set_ylabel("Misfit", fontsize=25)
        # ax3.tick_params(axis='x', labelsize=20)
        # ax3.tick_params(axis='y', labelsize=20)
        # ax3.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        # plt.legend(fontsize=20)
        # plt.tight_layout()
        # plt.yscale('log')
        # plt.ylim((pow(20, -1), pow(10, 1)))
        #
        # ax4 = plt.subplot(514)
        # ax4.plot(SR , label = 'SR, Mean:%.2f Min%.2f' %(np.mean(SR),SR.min()))
        # ax4.set_ylabel("Mean", fontsize=25)
        # ax4.tick_params(axis='x', labelsize=20)
        # ax4.tick_params(axis='y', labelsize=20)
        # ax4.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        # plt.legend(fontsize=20)
        # plt.tight_layout()
        # plt.yscale('log')
        # plt.ylim((pow(20, -1), pow(10, 1)))
        #
        # ax5 = plt.subplot(515)
        # ax5.plot( ST , label = 'ST, Mean:%.2f Min%.2f' %(np.mean(ST),ST.min()))
        # ax5.set_ylabel("Mean", fontsize=25)
        # ax5.set_xlabel("Sample", fontsize=25)
        # ax5.tick_params(axis='x', labelsize=20)
        # ax5.tick_params(axis='y', labelsize=20)
        # ax5.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        # plt.legend(fontsize=20)
        # plt.tight_layout()
        # plt.yscale('log')
        # plt.ylim((pow(20, -1), pow(10, 1)))
        # # plt.show()
        # if show == True:
        #     plt.show()
        # else:
        #     plt.savefig(dir + '/misfits.pdf')
        #     plt.close()

    def Get_cross_plots(self, filepath, savename, directory, skiprows, burnin,column_names, show=True):
        dir = directory + '/%s' % (savename)
        if not os.path.exists(dir):
            os.makedirs(dir)

        if filepath.endswith('.yaml') == True:
            with open(filepath, 'r') as stream:
                data_file = yaml.load(stream)
                stream.close()
                data = data_file['data']
                # parameters = data_file['parameters']
        else:
            # parameters = open(filepath, "r").readlines()[:33]
            data = np.loadtxt(filepath, delimiter=',', skiprows=skiprows)

        df = pd.DataFrame(data,
                          columns=column_names)
        #### PLOT:
        PZ = (df['p_z'].values)
        PR = (df['p_r'].values / 0.1)
        SZ = (df['s_z'].values / 0.14)
        SR = (df['s_r'].values / 0.35)
        ST = (df['s_t'].values)

        n_lowest = 10000
        Misfit = (PZ+PR+SZ+SR+ST)
        lowest_indices = Misfit.argsort()[0:n_lowest]
        # lowest_indices = df['Total_misfit'].values.argsort()[0:n_lowest]

        # PZ = PZ[PZ.argsort()[0:n_lowest]] #PZ[lowest_indices]
        # PR = PR[PR.argsort()[0:n_lowest]] #PR[lowest_indices]
        # SZ = SZ[SZ.argsort()[0:n_lowest]] #SZ[lowest_indices]
        # SR = SR[SR.argsort()[0:n_lowest]] #SR[lowest_indices]
        # ST = ST[ST.argsort()[0:n_lowest]] #ST[lowest_indices]

        PZ = PZ #PZ[lowest_indices]
        PR = PR #PR[lowest_indices]
        SZ = SZ #SZ[lowest_indices]
        SR = SR #SR[lowest_indices]
        ST = ST #ST[lowest_indices]
        TOTAL = Misfit

        Data = np.array([PZ,PR,SZ,SR,ST,TOTAL])
        dimensions_list = [0, 1, 2, 3, 4, 5]

        self.marginal_grid(savename,directory,Data.T,dimensions_list,bins=80,show = False)

    def full_moment_traces(self, filepath, savename, directory, skiprows, column_names, real_v, burnin):
        dir = directory + '/%s' % (savename.strip('.yaml'))
        if not os.path.exists(dir):
            os.makedirs(dir)

        if filepath.endswith('.yaml') == True:
            with open(filepath, 'r') as stream:
                data_file = yaml.load(stream)
                stream.close()
                data = data_file['data']
                # parameters = data_file['parameters']
        else:
            # parameters = open(filepath, "r").readlines()[:33]
            data = np.loadtxt(filepath, delimiter=',', skiprows=skiprows)

        length = len(data[0]) - (len(column_names) + 3)
        L_length = int(data[0][-1])
        R_length = int(data[0][-2])

        for i in range(R_length):
            column_names = np.append(column_names, "R_%i" % (i + 1))
        for i in range(L_length):
            column_names = np.append(column_names, "L_%i" % (i + 1))
        column_names = np.append(column_names, "Accepted")
        column_names = np.append(column_names, "Rayleigh_length")
        column_names = np.append(column_names, "Love_length")

        df = pd.DataFrame(data, columns=column_names)
        strike, dip, rake = aux_plane(real_v[2], real_v[3], real_v[4])

        M_true = self.convert(real_v[2], real_v[3], real_v[4], real_v[5])
        M_aux = self.convert(strike, dip, rake, real_v[5])

        moment = self.convert(df['Strike'][burnin:], df['Dip'][burnin:], df['Rake'][burnin:], df['M0'][burnin:])
        # moment = self.convert(df['Strike'][burnin:],df['Dip'][burnin:],df['Rake'][burnin:],1e+16)
        moment_names = np.array(["Mxx", "Myy", "Mzz", "Mxy", "Mxz", "Myz"])

        # beachball(moment, size=200, linewidth=2, facecolor='b', outfile=dir + '/beach.pdf')
        # beachball(M_true, size=200, linewidth=2, facecolor='b', outfile=dir + '/beach_true.pdf')
        Trace = np.mean(moment[0][burnin:]) + np.mean(moment[1][burnin:]) + np.mean(moment[2][burnin:])
        print("Mxx %.2E" % np.mean(moment[0]))
        print("Myy %.2E" % np.mean(moment[1]))
        print("Mzz %.2E" % np.mean(moment[2]))
        print("Trace %.2f" % Trace)
        column = 0

        params = {'legend.fontsize': 'x-large',
                  'figure.figsize': (15, 15),
                  'axes.labelsize': 20,
                  'axes.titlesize': 'x-large',
                  'xtick.labelsize': 14,
                  'ytick.labelsize': 13}
        pylab.rcParams.update(params)
        #

        fig = plt.figure(figsize=(25, 4))
        for i, v in enumerate(moment):
            ax1 = plt.subplot2grid((1, len(moment)), (0, column))
            plt.hist(v[burnin:], bins=100)
            ymin, ymax = ax1.get_ylim()
            xmin, xmax = ax1.get_xlim()
            plt.vlines(M_aux[i], ymin=ymin / 1.2, ymax=ymax, colors='k', linewidth=3, label='True model')
            # plt.text(0,ymax,"SD: %.2E" % np.std(v),fontsize =12)
            # plt.vlines(v[0], ymin=ymin, ymax=ymax, colors='r', linewidth=2)
            # plt.vlines(M_aux[i], ymin=ymin, ymax=ymax, colors='g', linewidth=2)
            # ax1.set_title("Density %s" % moment_names[i], color='b', fontsize=20)
            ax1.title.set_text("  SD: %.2E" % np.std(v))
            ax1.title.set_fontsize(18)
            ax1.set_xlabel("%s[Nm]" % moment_names[i], fontsize=18, color='b')
            if i == 0:
                ax1.set_ylabel("Posterior marginal", fontsize=20)
            ax1.tick_params(axis='x', labelsize=20)
            ax1.tick_params(axis='y', labelsize=20)
            ax1.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2), fontsize=18)
            # ax1.set_xlim(-0.5e17, 0.5e17)
            ax1.set_xlim(xmin, xmax)
            # plt.legend()
            plt.tight_layout()
            column += 1

        plt.savefig(dir + '/Full_moment_trace.pdf')

    def convert(self, strike, dip, rake, M0):
        rdip = np.deg2rad(strike)
        rstr = np.deg2rad(dip)
        rrake = np.deg2rad(rake)

        nx = -np.sin(rdip) * np.sin(rstr)
        ny = np.sin(rdip) * np.cos(rstr)
        nz = -np.cos(rdip)

        dx = np.cos(rrake) * np.cos(rstr) + np.cos(rdip) * np.sin(rrake) * np.sin(rstr)
        dy = np.cos(rrake) * np.sin(rstr) - np.cos(rdip) * np.sin(rrake) * np.cos(rstr)
        dz = -np.sin(rdip) * np.sin(rrake)

        # dx =  np.cos(rrake)*np.cos(rstr)+np.cos(rdip)*np.sin(rrake)*np.sin(rstr)
        # dy =  -np.cos(rrake)*np.sin(rstr)+np.cos(rdip)*np.sin(rrake)*np.cos(rstr)
        # dz = np.sin(rdip) *np.sin(rrake)

        Mxx = M0 * 2 * dx * nx
        Myy = M0 * 2 * dy * ny
        Mzz = M0 * 2 * dz * nz
        Mxy = M0 * dx * ny + dy * nx
        Mxz = M0 * dx * nz + dz * nx
        Myz = M0 * dy * nz + dz * ny

        moment = np.array([Mxx, Myy, Mzz, Mxy, Mxz, Myz])
        return moment

    def trace(self, filepath, savename, directory, skiprows, column_names,  burnin, lowest,real_v = None,):
        dir = directory + '/%s' % (savename.strip('.yaml'))
        if not os.path.exists(dir):
            os.makedirs(dir)

        if filepath.endswith('.yaml') == True:
            with open(filepath, 'r') as stream:
                data_file = yaml.load(stream)
                stream.close()
                data = data_file['data']
        else:
            data = np.loadtxt(filepath, delimiter=',', skiprows=skiprows)


        params = {'legend.fontsize': 'x-large',
                  'figure.figsize': (15, 15),
                  'axes.labelsize': 25,
                  'axes.titlesize': 'x-large',
                  'xtick.labelsize': 25,
                  'ytick.labelsize': 25}
        pylab.rcParams.update(params)
        #
        df = pd.DataFrame(data,
                          columns=column_names)
        df_select = df[["Epi", "Depth", "Strike", "Dip", "Rake"]]
        if real_v is not None:
            strike, dip, rake = aux_plane(real_v[2], real_v[3], real_v[4])
            M0_true = real_v[5]

        fig = plt.figure(figsize=(25, 6))
        # row = 0
        #
        # n_lowest = 200
        # lowest_indices = df[lowest].values.argsort()#[0:n_lowest]


        # PZ = (df['p_z'].values)
        # PR = (df['p_r'].values) * 0.01# / 0.1)
        # SZ = (df['s_z'].values)# / 0.14)
        # SR = (df['s_r'].values)# / 0.35)
        # ST = (df['s_t'].values)
        #
        # Misfit = (PZ)
        # lowest_indices = Misfit.argsort()
        #
        #
        # strike, dip, rake = aux_plane(df_select['Strike'][lowest_indices[0]], df_select['Dip'][lowest_indices[0]],
        #                               df_select['Rake'][lowest_indices[0]])
        # print(df_select['Strike'][lowest_indices[0]], df_select['Dip'][lowest_indices[0]], df_select['Rake'][lowest_indices[0]])
        # print(strike,dip,rake)

        #### TAYAK:
        # PZ = (df['p_z'].values)
        # PR = df['p_r'].values* 0.017
        # SZ = df['s_z'].values * 0.019
        # SR = df['s_r'].values* 0.0556
        # ST = (df['s_t'].values)  * 2.33
        # # AMP = df['Xi_amp'].values * 100

        ## TAYAK:
        # indices = np.where(Misfit < 0.13)
        # indices = np.where(Misfit < 0.27)
        # indices = np.where(Misfit < 1.25)

        ## DWTHot
        # PZ = (df['p_z'].values)
        # PR = (df['p_r'].values) * 0.1
        # SZ = (df['s_z'].values) * 0.055
        # SR = (df['s_r'].values) * 0.21
        # ST = (df['s_t'].values) * 0.2
        # AMP =(( (np.abs(np.log10(df['Npz'].values ) - np.log10(df['Nst'].values ))) / np.log(2)) ) ** (1/4)
        #
        # Misfit = PZ + ST + PR + SZ +SR
        #
        # ## DWTHot
        # # indices = np.where(Misfit < 0.13)
        # indices = np.where(Misfit < 5.95)
        # misfits = Misfit[indices]
        #
        # strike = df['Strike'].values[indices]
        # dip = df['Dip'].values[indices]
        # rake = df['Rake'].values[indices]

        ## EH45Tcold
        PZ = (df['p_z'].values)
        PR = (df['p_r'].values) * 0.5
        SZ = (df['s_z'].values) * 0.16
        SR = (df['s_r'].values) * 0.15
        ST = (df['s_t'].values) * 0.94
        AMP =(( (np.abs(np.log10(df['Npz'].values ) - np.log10(df['Nst'].values ))) / np.log(2)) ) ** (1/4) * 1.3

        Misfit = PZ +ST + PR + SZ + SR #+ AMP

        # indices = np.where(Misfit < 0.5)
        # indices = np.where(Misfit < 0.7)
        indices = np.where(Misfit < 8.75)
        misfits = Misfit[indices]

        strike = df['Strike'].values[indices]
        dip = df['Dip'].values[indices]
        rake = df['Rake'].values[indices]



        lowest_indice = misfits.argsort()[0]
        strike_fault= strike[lowest_indice]
        dip_fault= dip[lowest_indice]
        rake_fault= rake[lowest_indice]
        strike_aux, dip_aux, rake_aux = aux_plane(strike_fault,dip_fault,rake_fault)


        ax1 = plt.subplot2grid((1, 3), (0, 0))
        bin = int(360 * (len(strike)**(1/3) / (3.49 * np.std(strike))))
        binn = int(360 * (len(df['Strike'].values)**(1/3) / (3.49 * np.std(df['Strike'].values))))
        # plt.hist(df['Strike'][burnin:], weights=(df['Total_misfit'][burnin:].max() -  df['Total_misfit'][burnin:]) ,density = True,bins= 80, alpha=0.8)
        # plt.hist(df_select['Strike'][lowest_indices], bins= 80, alpha=0.8, label='Lowest Misfit')
        # plt.hist(df_select['Strike'][lowest_indices],weights=np.linspace(1,0,len(df_select['Strike'][lowest_indices])), bins= 80, alpha=0.8, label='Lowest Misfit')

        # n, x, _ = plt.hist(df['Strike'].values,weights=Misfit.max() - Misfit,density = True, bins= binn, alpha=0.4, color = 'grey')
        n, x, _ = plt.hist(strike, weights=Misfit.max() - misfits,density=True, bins=20, alpha=0.8, color='grey',
                           label='Lowest Misfit')
        # bin_centers = 0.5 * (x[1:] + x[:-1])
        # plt.plot(bin_centers, n)  ## using bin_centers rather than edges

        ymin, ymax = ax1.get_ylim()
        plt.vlines(strike_fault, ymin=ymin, ymax=ymax, colors='g', linewidth=3, label='Fault plane')
        plt.vlines(strike_aux, ymin=ymin, ymax=ymax, colors='k', linewidth=3, label='Auxiliary plane')

        # plt.vlines(df_select['Strike'][lowest_indices[0]], ymin=ymin, ymax=ymax, colors='g', linewidth=3, label='Fault plane')
        # plt.vlines(strike, ymin=ymin, ymax=ymax, colors='k', linewidth=3, label='Auxiliary plane')
        if real_v is not None:
            plt.vlines(real_v[2], ymin=ymin, ymax=ymax, colors='g', linewidth=3, label='Auxiliary plane')
            plt.vlines(strike, ymin=ymin, ymax=ymax, colors='k', linewidth=3, label='Fault plane')
        ax1.set_title("Density Strike", color='b', fontsize=25)
        ax1.set_xlabel("N=%i" % (len(df_select['Strike'])), fontsize=25)
        ax1.set_ylabel("Posterior marginal", fontsize=25)
        ax1.tick_params(axis='x', labelsize=20)
        ax1.tick_params(axis='y', labelsize=20)
        ax1.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        ax1.set_xlim(0, 360)
        # plt.legend( fontsize=20)
        plt.tight_layout()

        ax2 = plt.subplot2grid((1, 3), (0, 1))
        bin = int(90 * (len(dip) ** (1 / 3) / (3.49 * np.std(dip))))
        binn = int(90 * (len(df['Dip'].values) ** (1 / 3) / (3.49 * np.std(df['Dip'].values))))
        # plt.hist(df['Dip'][burnin:], weights=df['Total_misfit'][burnin:].max() -  df['Total_misfit'][burnin:],density = True,bins= 80, alpha=0.8)
        # plt.hist(df_select['Dip'][lowest_indices], bins=80, alpha=0.8, label='Lowest Misfit')
        # plt.hist(df_select['Dip'][lowest_indices],
        #          weights=np.linspace(1, 0, len(df_select['Dip'][lowest_indices])), bins=80, alpha=0.8,
        #          label='Lowest Misfit')
        # plt.hist(df['Dip'].values,weights=Misfit.max() - Misfit,bins=30, density = True,alpha=0.4, color = 'grey')
        plt.hist(dip,weights=Misfit.max() - misfits, bins=30, density = True,alpha=0.8, color = 'grey')

        ymin, ymax = ax2.get_ylim()
        plt.vlines(dip_fault, ymin=ymin, ymax=ymax, colors='g', linewidth=3, label='Fault plane')
        plt.vlines(dip_aux, ymin=ymin, ymax=ymax, colors='k', linewidth=3, label='Auxiliary plane')
        # plt.vlines(df_select['Dip'][lowest_indices[0]], ymin=ymin, ymax=ymax, colors='g', linewidth=3, label='Fault plane')
        # plt.vlines(dip, ymin=ymin, ymax=ymax, colors='k', linewidth=3, label='Auxiliary plane')
        if real_v is not None:
            plt.vlines(real_v[3], ymin=ymin, ymax=ymax, colors='g', linewidth=3)
            plt.vlines(dip, ymin=ymin, ymax=ymax, colors='k', linewidth=3, label='True model')
        # plt.vlines(df_select['Dip'][1], ymin=ymin, ymax=ymax, colors='r', linewidth=3,label='Start model')
        ax2.set_title("Density Dip", color='b', fontsize=25)
        ax2.set_xlabel("N=%i" % (len(df_select['Dip'])), fontsize=25)
        ax2.xaxis.set_ticks(np.arange(0, 90, 10))
        ax2.tick_params(axis='x', labelsize=20)
        ax2.tick_params(axis='y', labelsize=20)
        ax2.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        ax2.set_xlim(0, 90)
        # plt.legend( fontsize=20)
        plt.tight_layout()

        ax3 = plt.subplot2grid((1, 3), (0, 2))
        bin = int(360 * (len(rake) ** (1 / 3) / (3.49 * np.std(rake))))
        binn = int(360 * (len(df['Rake'].values) ** (1 / 3) / (3.49 * np.std(df['Rake'].values))))
        # plt.hist(df['Rake'][burnin:], weights=df['Total_misfit'][burnin:].max() -  df['Total_misfit'][burnin:],density = True,bins= 80, alpha=0.8)
        # plt.hist(df_select['Rake'][lowest_indices], bins=80, alpha=0.8, label='Lowest Misfit')
        # plt.hist(df_select['Rake'][lowest_indices],
        #          weights=np.linspace(1, 0, len(df_select['Rake'][lowest_indices])), bins=80, alpha=0.8,
        #          label='Lowest Misfit')
        #
        # plt.hist(df['Rake'].values,weights=Misfit.max() - Misfit, bins=binn, density = True,alpha=0.4,color = 'grey',
        #          label='All values')
        plt.hist(rake,weights=Misfit.max() - misfits, bins=30, density = True,alpha=0.8,color = 'grey',
                 label='Clustered Groups')
        ymin, ymax = ax3.get_ylim()
        plt.vlines(rake_fault, ymin=ymin, ymax=ymax, colors='g', linewidth=3, label='Fault plane')
        plt.vlines(rake_aux, ymin=ymin, ymax=ymax, colors='k', linewidth=3, label='Auxiliary plane')
        # plt.vlines(df_select['Rake'][lowest_indices[0]], ymin=ymin, ymax=ymax, colors='g', linewidth=3, label='Fault plane')
        # plt.vlines(rake, ymin=ymin, ymax=ymax, colors='k', linewidth=3, label='Auxiliary plane')
        if real_v is not None:
            plt.vlines(real_v[4], ymin=ymin, ymax=ymax, colors='g', linewidth=3, label='Auxiliary plane')
            plt.vlines(rake, ymin=ymin, ymax=ymax, colors='k', linewidth=3, label='Fault plane')
            plt.legend(loc='lower left', fontsize=18)
            # plt.legend(loc='center left', bbox_to_anchor=(1.1, 0.5))
        ax3.set_title("Density Rake", color='b', fontsize=25)
        ax3.set_xlabel("N=%i" % (len(df_select['Rake'])), fontsize=25)
        ax3.tick_params(axis='x', labelsize=20)
        ax3.tick_params(axis='y', labelsize=20)
        ax3.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        ax3.set_xlim(-180, 180)
        plt.tight_layout()
        plt.legend(loc='upper right', fontsize=20)
        # plt.show()
        plt.savefig(dir + '/Trace_fault.pdf')

        beachball([strike_fault, dip_fault, rake_fault], size=200, linewidth=2, facecolor='grey',
                  outfile=dir + '/beach_true.pdf')

        fig = self.plot_BBB(strike, dip, rake,color = 'grey')
        plt.savefig(dir + '/BBB.pdf')



        #
        # fig = plt.figure(figsize=(25, 6))
        #
        # ax4 = plt.subplot2grid((1, 3), (0, 0))
        # plt.hist(df_select['Depth'][burnin:] / 1000, bins=50, alpha=0.8)
        # ymin, ymax = ax4.get_ylim()
        # xmin, xmax = ax4.get_xlim()
        # if real_v is not None:
        #     plt.vlines(real_v[1], ymin=ymin, ymax=ymax, colors='k', linewidth=3, label='True model')
        #
        # ax4.set_title("Density Depth", color='b', fontsize=20)
        # # ax4.set_xlabel("N=%i" % (len(df_select['Depth'])), fontsize=20)
        # ax4.set_xlabel("Depth [km]", fontsize=20)
        # ax4.set_ylabel("Posterior marginal", fontsize=20)
        # ax4.tick_params(axis='x', labelsize=20)
        # ax4.tick_params(axis='y', labelsize=20)
        # ax4.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        # ax4.ticklabel_format(style="sci", axis='x', scilimits=(-2, 2))
        #
        # # ax4.set_xlim(real_v[1]-10000, real_v[1]+10000)
        # ax4.set_xlim(xmin, xmax)
        # # plt.legend( fontsize=20)
        # plt.tight_layout()
        #
        # ax5 = plt.subplot2grid((1, 3), (0, 1))
        # plt.hist(df_select['Epi'][burnin:], bins=50, alpha=0.8)
        # # plt.hist(df['Epi'][burnin:],bins=100,alpha=0.5,color = 'b',label = 'MAAK')
        # # plt.hist(df2['Epi'][burnin:],bins=100,alpha=0.5,color = 'r',label = 'EH45TcoldCrust_1b')
        # if real_v is not None:
        #     ymin, ymax = ax5.get_ylim()
        #
        # # plt.vlines(df_select['Epi'][0], ymin=ymin, ymax=ymax, colors='r', linewidth=3, label='Start model')
        #     plt.vlines(real_v[0], ymin=ymin, ymax=ymax, colors='k', linewidth=3,label = 'True model')
        # ax5.set_title("Density Epicentral Distance", color='b', fontsize=20)
        # ax5.set_xlabel("N=%i" % (len(df_select['Epi'])), fontsize=20)
        # ax5.tick_params(axis='x', labelsize=20)
        # ax5.tick_params(axis='y', labelsize=20)
        # ax5.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        # # plt.legend(fontsize=20)
        # plt.tight_layout()
        #
        # ax6 = plt.subplot2grid((1, 3), (0, 2))
        # Mw = 2.0 / 3.0 * (np.log10(df['M0'][burnin:]) - 9.1)
        # plt.hist(Mw, bins=np.arange(3, 5, 0.05), alpha=0.8)
        # # plt.hist(df['M0'][burnin:],bins=100,alpha=0.5,color = 'b',label = 'MAAK')
        # # plt.hist(df2['M0'][burnin:],bins=100,alpha=0.5,color = 'r',label = 'EH45TcoldCrust_1b')
        # # if real_v is not None:
        # #     ymin, ymax = ax6.get_ylim()
        # #     plt.vlines(M0_true, ymin=ymin, ymax=ymax, colors='k', linewidth=3, label= 'True model')
        # # plt.vlines(df['M0'][1], ymin=ymin, ymax=ymax, colors='r', linewidth=3,label = 'Start model')
        # ax6.set_title("Density Moment Magnitude", color='b', fontsize=20)
        # ax6.set_xlabel("N=%i" % (len(df['M0'])), fontsize=18)
        # ax6.tick_params(axis='x', labelsize=18)
        # ax6.tick_params(axis='y', labelsize=18)
        # ax6.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        # # ax6.set_xscale('log')
        # # ax6.set_xlim(0e18, 0.5e18 )
        # # plt.legend(loc='center left', bbox_to_anchor=(1.1, 0.5))
        # plt.legend(loc='upper right', fontsize=20)
        # plt.tight_layout()
        #
        # # plt.show()
        # plt.savefig(dir + '/Trace_position.pdf')

        # n_lowest = 10000
        # # pos = np.argmin(df['Total_misfit'].values)
        # lowest_indices = df[lowest].values.argsort()[0:n_lowest]
        # lowest_misfits = df[lowest].values#[lowest_indices]
        # lowest_strike = df['Strike'].values#[lowest_indices]
        # lowest_dip = df['Dip'].values#[lowest_indices]
        # lowest_rake = df['Rake'].values#[lowest_indices]
        # lowest_depth = df['Depth'].values#[lowest_indices]
        # lowest_P_shift = df['Shift_P'].values#[lowest_indices]
        # lowest_S_shift = df['Shift_S'].values#[lowest_indices]
        #
        # ind_Negative_Pshift = np.where(lowest_P_shift < 0)
        # ind_Positve_Pshift = np.where(lowest_P_shift > 0)
        #
        # fig = plt.figure(figsize=(20, 12))
        # row = 0
        #
        # ax1 = plt.subplot(231)
        # plt.plot(lowest_strike[ind_Negative_Pshift],lowest_misfits[ind_Negative_Pshift],'bo', label = 'Negative P - shift')
        # plt.plot(lowest_strike[ind_Positve_Pshift],lowest_misfits[ind_Positve_Pshift],'ro', label = 'Positive P - shift')
        # # plt.plot(df['Strike'].values,df['Total_misfit'].values,'bo')
        # ymin, ymax = ax1.get_ylim()
        # if real_v is not None:
        #     plt.vlines(real_v[2], ymin=ymin, ymax=ymax, colors='g', linewidth=3, label='Auxiliary plane')
        #     plt.vlines(strike, ymin=ymin, ymax=ymax, colors='k', linewidth=3, label='Fault plane')
        # ax1.set_title("Strike", color='b', fontsize=25)
        # ax1.set_xlabel("N=%i" % (len(df_select['Strike'])), fontsize=25)
        # ax1.set_ylabel("Misfit %s" % lowest, fontsize=25)
        # ax1.tick_params(axis='x', labelsize=20)
        # ax1.tick_params(axis='y', labelsize=20)
        # ax1.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        # ax1.set_xlim(0, 360)
        # plt.legend( fontsize=20)
        # plt.tight_layout()
        #
        # ax2 = plt.subplot(232)
        # plt.plot(lowest_dip[ind_Negative_Pshift],lowest_misfits[ind_Negative_Pshift],'bo', label = 'Negative P - shift')
        # plt.plot(lowest_dip[ind_Positve_Pshift],lowest_misfits[ind_Positve_Pshift],'ro', label = 'Positive P - shift')
        # # plt.plot(df['Dip'].values,df['Total_misfit'].values,'bo')
        # ymin, ymax = ax2.get_ylim()
        # if real_v is not None:
        #     plt.vlines(real_v[3], ymin=ymin, ymax=ymax, colors='g', linewidth=3)
        #     plt.vlines(dip, ymin=ymin, ymax=ymax, colors='k', linewidth=3, label='True model')
        # # plt.vlines(df_select['Dip'][1], ymin=ymin, ymax=ymax, colors='r', linewidth=3,label='Start model')
        # ax2.set_title("Dip", color='b', fontsize=25)
        # ax2.set_xlabel("N=%i" % (len(df_select['Dip'])), fontsize=25)
        # ax2.xaxis.set_ticks(np.arange(0, 90, 10))
        # ax2.tick_params(axis='x', labelsize=20)
        # ax2.tick_params(axis='y', labelsize=20)
        # ax2.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        # ax2.set_xlim(0, 90)
        # # plt.legend( fontsize=20)
        # plt.tight_layout()
        #
        # ax3 = plt.subplot(233)
        # plt.plot(lowest_rake[ind_Negative_Pshift], lowest_misfits[ind_Negative_Pshift], 'bo',
        #          label='Negative P - shift')
        # plt.plot(lowest_rake[ind_Positve_Pshift], lowest_misfits[ind_Positve_Pshift], 'ro',
        #          label='Positive P - shift')
        # # plt.plot(df['Rake'].values,df['Total_misfit'].values,'bo')
        # ymin, ymax = ax3.get_ylim()
        # if real_v is not None:
        #     plt.vlines(real_v[4], ymin=ymin, ymax=ymax, colors='g', linewidth=3, label='Auxiliary plane')
        #     plt.vlines(rake, ymin=ymin, ymax=ymax, colors='k', linewidth=3, label='Fault plane')
        #     plt.legend(loc='upper left', fontsize=20)
        #     # plt.legend(loc='center left', bbox_to_anchor=(1.1, 0.5))
        # ax3.set_title("Rake", color='b', fontsize=25)
        # ax3.set_xlabel("N=%i" % (len(df_select['Rake'])), fontsize=25)
        # ax3.tick_params(axis='x', labelsize=20)
        # ax3.tick_params(axis='y', labelsize=20)
        # ax3.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        # ax3.set_xlim(-180, 180)
        # plt.tight_layout()
        #
        #
        # ax4 = plt.subplot(234)
        # plt.plot(lowest_depth,lowest_misfits,'bo')
        # # plt.plot(df['Depth'].values,df['Total_misfit'].values,'bo')
        # ymin, ymax = ax4.get_ylim()
        # if real_v is not None:
        #     plt.vlines(real_v[1], ymin=ymin, ymax=ymax, colors='g', linewidth=3, label='Depth')
        #     plt.legend(loc='upper left', fontsize=20)
        #     # plt.legend(loc='center left', bbox_to_anchor=(1.1, 0.5))
        # ax4.set_title("Depth", color='b', fontsize=25)
        # ax4.set_xlabel("N=%i" % (len(df_select['Rake'])), fontsize=25)
        # ax4.tick_params(axis='x', labelsize=20)
        # ax4.tick_params(axis='y', labelsize=20)
        # ax4.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        # ax4.set_xlim(20000, 100000)
        # plt.tight_layout()
        #
        # ax5 = plt.subplot(235)
        # plt.plot(lowest_P_shift[ind_Negative_Pshift], lowest_misfits[ind_Negative_Pshift], 'bo',
        #          label='Negative P - shift')
        # plt.plot(lowest_P_shift[ind_Positve_Pshift], lowest_misfits[ind_Positve_Pshift], 'ro',
        #          label='Positive P - shift')
        # # plt.plot(df['Shift_P'].values,df['Total_misfit'].values,'bo')
        # ax5.set_title("P_Shift", color='b', fontsize=25)
        # ax5.set_xlabel("N=%i" % (len(df_select['Rake'])), fontsize=25)
        # ax5.tick_params(axis='x', labelsize=20)
        # ax5.tick_params(axis='y', labelsize=20)
        # ax5.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        # ax5.set_xlim(-60, 60)
        # plt.tight_layout()
        #
        # ax6 = plt.subplot(236)
        # plt.plot(lowest_S_shift[ind_Negative_Pshift], lowest_misfits[ind_Negative_Pshift], 'bo',
        #          label='Negative P - shift')
        # plt.plot(lowest_S_shift[ind_Positve_Pshift], lowest_misfits[ind_Positve_Pshift], 'ro',
        #          label='Positive P - shift')
        # # plt.plot(df['Shift_S'].values,df['Total_misfit'].values,'bo')
        # ax6.set_title("S_Shift", color='b', fontsize=25)
        # ax6.set_xlabel("N=%i" % (len(df_select['Rake'])), fontsize=25)
        # ax6.tick_params(axis='x', labelsize=20)
        # ax6.tick_params(axis='y', labelsize=20)
        # ax6.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        # ax6.set_xlim(-60, 60)
        # plt.tight_layout()
        #
        # # plt.show()
        # plt.savefig(dir + '/Misfit_vs_Parameters.pdf')

        #
        # fig = plt.figure(figsize=(10, 6))
        # row = 0
        #
        # ax1 = plt.subplot2grid((3,1), (0, 0))
        # plt.plot(df_select['Strike'])
        # ax1.set_title("Strike", color='b', fontsize=25)
        # ax1.set_xlabel("N=%i" % (len(df_select['Strike'])), fontsize=25)
        # ax1.set_ylabel("Misfit", fontsize=25)
        # ax1.tick_params(axis='x', labelsize=20)
        # ax1.tick_params(axis='y', labelsize=20)
        # ax1.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        # ax1.set_ylim(0, 360)
        # # plt.legend( fontsize=20)
        # plt.tight_layout()
        #
        # ax2 = plt.subplot2grid((3,1), (1,0))
        # plt.plot(df_select['Dip'])
        #
        # ax2.set_title("Dip", color='b', fontsize=25)
        # ax2.set_xlabel("N=%i" % (len(df_select['Dip'])), fontsize=25)
        # ax2.xaxis.set_ticks(np.arange(0, 90, 10))
        # ax2.tick_params(axis='x', labelsize=20)
        # ax2.tick_params(axis='y', labelsize=20)
        # ax2.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        # ax2.set_ylim(0, 90)
        # # plt.legend( fontsize=20)
        # plt.tight_layout()
        #
        # ax3 = plt.subplot2grid((3, 1), (2,0))
        # plt.plot(df_select['Rake'])
        # ax3.set_title("Rake", color='b', fontsize=25)
        # ax3.set_xlabel("N=%i" % (len(df_select['Rake'])), fontsize=25)
        # ax3.tick_params(axis='x', labelsize=20)
        # ax3.tick_params(axis='y', labelsize=20)
        # ax3.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        # ax3.set_ylim(-180, 180)
        # plt.tight_layout()
        # # plt.show()
        # plt.savefig(dir + '/Trace_plot.pdf')

    def get_beachballs(self, strike, dip, rake, M0, savepath):
        rdip = np.deg2rad(dip)
        rstr = np.deg2rad(strike)
        rrake = np.deg2rad(rake)

        nx = -np.sin(rdip) * np.sin(rstr)
        ny = np.sin(rdip) * np.cos(rstr)
        nz = -np.cos(rdip)

        dx = np.cos(rrake) * np.cos(rstr) + np.cos(rdip) * np.sin(rrake) * np.sin(rstr)
        dy = np.cos(rrake) * np.sin(rstr) - np.cos(rdip) * np.sin(rrake) * np.cos(rstr)
        dz = -np.sin(rdip) * np.sin(rrake)

        # dx =  np.cos(rrake)*np.cos(rstr)+np.cos(rdip)*np.sin(rrake)*np.sin(rstr)
        # dy =  -np.cos(rrake)*np.sin(rstr)+np.cos(rdip)*np.sin(rrake)*np.cos(rstr)
        # dz = np.sin(rdip) *np.sin(rrake)

        Mxx = M0 * 2 * dx * nx
        Myy = M0 * 2 * dy * ny
        Mzz = M0 * 2 * dz * nz
        Mxy = M0 * dx * ny + dy * nx
        Mxz = M0 * dx * nz + dz * nx
        Myz = M0 * dy * nz + dz * ny

        moment = np.array([Mxx, Myy, Mzz, Mxy, Mxz, Myz])
        stations = {'names': ['S01', 'S02', 'S03', 'S04'],
                    'azimuth': np.array([120., 5., 250., 75.]),
                    'takeoff_angle': np.array([30., 60., 45., 10.]),
                    'polarity': np.array([0.8, 0.5, 0.7, -0.9])}
        # MTplot(np.array([[1], [0], [-1], [0], [0], [0]]), 'beachball',  stations=stations, fault_plane=True)
        beachball(moment, size=200, linewidth=2, facecolor='b', outfile=savepath)

    def get_waveforms(self, filepath, savename, directory, skiprows, column_names, real_v, burnin, prior):
        dir = directory + '/%s' % (savename.strip('.yaml'))
        if not os.path.exists(dir):
            os.makedirs(dir)

        if filepath.endswith('.yaml') == True:
            with open(filepath, 'r') as stream:
                data_file = yaml.load(stream)
                stream.close()
                data = data_file['data']
                # parameters = data_file['parameters']
        else:
            # parameters = open(filepath, "r").readlines()[:33]
            data = np.loadtxt(filepath, delimiter=',', skiprows=0)

        # length = len(data[0]) - (len(column_names) + 3)
        # L_length = int(data[0][-1])
        # R_length = int(data[0][-2])

        # for i in range(R_length):
        #     column_names = np.append(column_names, "R_%i" % (i + 1))
        # for i in range(L_length):
        #     column_names = np.append(column_names, "L_%i" % (i + 1))
        # column_names = np.append(column_names, "Accepted")
        # column_names = np.append(column_names, "Rayleigh_length")
        # column_names = np.append(column_names, "Love_length")

        df = pd.DataFrame(data, columns=column_names)

        fig_bb, ax_bb = plt.subplots(1, 1, figsize=(4, 4))

        ax_bb.set_xticks([])
        ax_bb.set_yticks([])
        ax_bb.axis('off')
        img = None
        buf = io.BytesIO()

        epi = df['Epi']
        depth = df['Depth']
        strike = df['Strike']
        dip = df['Dip']
        rake = df['Rake']
        M0 = df['M0']

        P_shift = df['Shift_P']
        S_shift = df['Shift_S']

        from Get_Seismogram import Get_Seismogram
        from Cut_windows import Cut_windows
        seis = Get_Seismogram(prior)

        stream = obspy.read(prior['Observed'])

        BW_obs = Cut_windows(prior['VELOC_taup'])
        BW_obs.Get_bw_windows(stream, prior['Real_epi'], prior['Real_depth'], prior['or_time'], prior['npts'])
        # tt_P = obspy.UTCDateTime(2019, 1, 3, 15, 9, 54.9)
        # tt_S = obspy.UTCDateTime(2019, 1, 3, 15, 18, 34.6)
        # BW_obs.Get_bw_windows_MANUAL(stream, tt_P, tt_S, prior['or_time'], npts=prior['npts'])

        BW_syn = Cut_windows(prior['VELOC_taup'])

        fig = plt.figure(figsize=(10, 10))
        delta = BW_obs.P_stream.traces[0].meta.delta
        p_time_array = np.arange(len(BW_obs.P_stream.traces[0].data)) * delta
        s_time_array = np.arange(len(BW_obs.S_stream.traces[0].data)) * delta

        start_P = int((BW_obs.start_P.timestamp - prior['or_time'].timestamp - 50) / delta)
        end_P = int((BW_obs.start_P.timestamp - prior['or_time'].timestamp + 130) / delta)

        start_S = int((BW_obs.start_S.timestamp - prior['or_time'].timestamp - 50) / delta)
        end_S = int((BW_obs.start_S.timestamp - prior['or_time'].timestamp + 100) / delta)

        ax1 = plt.subplot2grid((5, 1), (0, 0))
        ax2 = plt.subplot2grid((5, 1), (1, 0))
        ax3 = plt.subplot2grid((5, 1), (2, 0))
        ax4 = plt.subplot2grid((5, 1), (3, 0))
        ax5 = plt.subplot2grid((5, 1), (4, 0))

        # a = np.where(np.logical_and(epi>=prior['Real_epi'] , epi <= prior['Real_epi']))
        # a = np.where(np.logical_and(epi>=88 -0.5  , epi <= 88 + 0.5))

        # for i in range(len(a[0])):
        for i in np.arange(len(epi) - 100, len(epi), 1):
            dict = geo.Geodesic(a=prior['radius'], f=prior['f']).ArcDirect(lat1=prior['la_r'], lon1=prior['lo_r'],
                                                                           azi1=prior['baz'], a12=epi[i], outmask=1929)

            st_syn = seis.get_seis_manual(la_s=dict['lat2'], lo_s=dict['lon2'], depth=depth[i],
                                          strike=strike[i], dip=dip[i], rake=rake[i],
                                          time=prior['or_time'], M0=M0[i])

            BW_syn.Get_bw_windows(st_syn, epi[i], depth[i], prior['or_time'], prior['npts'])

            P_shift_array = self.shift(BW_syn.P_stream.traces[0].data, -int(P_shift[i]))
            #
            # ax1.plot(p_time_array[start_P:end_P], BW_syn.P_stream.traces[0].data[start_P:end_P], 'g',
            #          label='Synthetic', linewidth = 0.1)

            ax1.plot(p_time_array[start_P:end_P], P_shift_array[start_P:end_P], 'r',
                     label='Synthetic', linewidth=0.1)
            # ax1.plot( P_shift_array, 'r',  label='Synthetic', linewidth = 0.1)
            ax1.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
            ax1.tick_params(axis='x', labelsize=18)
            ax1.tick_params(axis='y', labelsize=18)

            # plt.tight_layout()
            # plt.legend(loc='lower left', fontsize=15)

            P_shift_array = self.shift(BW_syn.P_stream.traces[1].data, -int(P_shift[i]))
            ax2.plot(p_time_array[start_P:end_P], P_shift_array[start_P:end_P], 'r',
                     label='Synthetic', linewidth=0.1)
            # ax2.plot(P_shift_array, 'r',label='Synthetic', linewidth = 0.1)
            ax2.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
            ax2.tick_params(axis='x', labelsize=18)
            ax2.tick_params(axis='y', labelsize=18)
            # plt.tight_layout()

            S_shift_array = self.shift(BW_syn.S_stream.traces[0].data, -int(S_shift[i]))
            ax3.plot(s_time_array[start_S:end_S], S_shift_array[start_S:end_S], 'r', linewidth=0.1)
            # ax3.plot(S_shift_array, 'r', linewidth = 0.1)
            ax3.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
            ax3.tick_params(axis='x', labelsize=18)
            ax3.tick_params(axis='y', labelsize=18)
            # plt.tight_layout()

            S_shift_array = self.shift(BW_syn.S_stream.traces[1].data, -int(S_shift[i]))
            ax4.plot(s_time_array[start_S:end_S], S_shift_array[start_S:end_S], 'r', linewidth=0.1)
            # ax4.plot(S_shift_array, 'r', linewidth = 0.1)
            ax4.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
            ax4.tick_params(axis='x', labelsize=18)
            ax4.tick_params(axis='y', labelsize=18)
            # plt.tight_layout()

            S_shift_array = self.shift(BW_syn.S_stream.traces[2].data, -int(S_shift[i]))
            ax5.plot(s_time_array[start_S:end_S], S_shift_array[start_S:end_S], 'r', linewidth=0.1)
            # ax5.plot( S_shift_array, 'r', linewidth = 0.1)
            ax5.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
            ax5.tick_params(axis='x', labelsize=18)
            ax5.tick_params(axis='y', labelsize=18)
            ax5.set_xlabel(BW_obs.start_P.strftime('From P arrival: %Y-%m-%dT%H:%M:%S + [sec]'), fontsize=18)

        ax1.plot(p_time_array[start_P:end_P], BW_obs.P_stream.traces[0].data[start_P:end_P], 'b')
        # ax1.plot( BW_obs.P_stream.traces[0].data, 'b', label='Observed', linewidth = 0.1)
        ymin, ymax = ax1.get_ylim()
        xmin, xmax = ax1.get_xlim()
        ax1.text(xmax - 5, ymax / 1.7, "P-Z", fontsize=20, color='b')

        ax2.plot(p_time_array[start_P:end_P], BW_obs.P_stream.traces[1].data[start_P:end_P], 'b')
        # ax2.plot(BW_obs.P_stream.traces[1].data, 'b', linewidth = 0.1)
        ymin, ymax = ax2.get_ylim()
        xmin, xmax = ax2.get_xlim()
        ax2.text(xmax - 5, ymax / 1.7, "P-R", fontsize=20, color='b')

        ax3.plot(s_time_array[start_S:end_S], BW_obs.S_stream.traces[0].data[start_S:end_S], 'b')
        # ax3.plot( BW_obs.S_stream.traces[0].data, 'b', linewidth = 0.1)
        ymin, ymax = ax3.get_ylim()
        xmin, xmax = ax3.get_xlim()
        ax3.text(xmax - 10, ymax / 1.5, "S-Z", fontsize=20, color='b')

        ax4.plot(s_time_array[start_S:end_S], BW_obs.S_stream.traces[1].data[start_S:end_S], 'b')
        # ax4.plot(BW_obs.S_stream.traces[1].data, 'b', linewidth = 0.1)
        ymin, ymax = ax4.get_ylim()
        xmin, xmax = ax4.get_xlim()
        ax4.text(xmax - 10, ymax / 1.5, "S-R", fontsize=20, color='b')

        ax5.plot(s_time_array[start_S:end_S], BW_obs.S_stream.traces[2].data[start_S:end_S], 'b')
        # ax5.plot(BW_obs.S_stream.traces[2].data, 'b', linewidth = 0.1)
        ymin, ymax = ax5.get_ylim()
        xmin, xmax = ax5.get_xlim()
        ax5.text(xmax - 10, ymax / 1.7, "S-T", fontsize=20, color='b')

        plt.tight_layout()

        # plt.show()
        plt.savefig(dir + '/Waveforms.pdf')
        # plt.show()
        plt.close()

    def get_BBB(self, filepath, savename, directory, skiprows, column_names, burnin, lowest,real_v = None):
        dir = directory + '/%s' % (savename.strip('.yaml'))
        if not os.path.exists(dir):
            os.makedirs(dir)
        if real_v is not None:
            beachball([real_v[2], real_v[3], real_v[4]], size=200, linewidth=2, facecolor='b',
                      outfile=dir + '/beach_true.pdf')

        if filepath.endswith('.yaml') == True:
            with open(filepath, 'r') as stream:
                data_file = yaml.load(stream)
                stream.close()
                data = data_file['data']
        else:
            data = np.loadtxt(filepath, delimiter=',', skiprows=skiprows)


        df = pd.DataFrame(data, columns=column_names)
        n_lowest = 100
        lowest_indices = df[lowest].values.argsort()[0:n_lowest]
        strike = df['Strike'].values[lowest_indices]
        dip = df['Dip'].values[lowest_indices]
        rake = df['Rake'].values[lowest_indices]
        depth = df['Depth'].values[lowest_indices]
        epi = df['Epi'].values[lowest_indices]

        ## TAYAK: EVENT_235_1
        # PZ = (df['p_z'].values)
        # PR = (df['p_r'].values) * 0.017
        # SZ = (df['s_z'].values) * 0.019
        # SR = (df['s_r'].values) * 0.0556
        # ST = (df['s_t'].values)  * 2.33
        # AMP = (((np.abs(np.log10(df['Npz'].values) - np.log10(df['Nst'].values))) / np.log(2))) ** (1 / 4) / 2
        # Misfit = PZ + ST + PR + SR + SZ + AMP
        # # indices = np.where(Misfit < 0.13)
        # # indices = np.where(Misfit < 0.27)
        # # indices = np.where(Misfit < 1.25)
        # indices = np.where(Misfit < 1.55)
        # #
        # strike = df['Strike'].values[indices]
        # dip = df['Dip'].values[indices]
        # rake = df['Rake'].values[indices]
        #### PZ
        # G1_S = np.hstack((strike[(strike < 70) ] , strike[(strike > 349) ] ))
        # G1_D = np.hstack(( dip[(strike < 70) ] ,  dip[(strike > 349) ] ))
        # G1_R = np.hstack(( rake[(strike < 70) ] , rake[(strike > 349) ] ))
        #
        # G21_S = strike[( ( (70 < strike) & (strike < 250)) & ((-180 < rake ) & (rake <-100)) )]
        # G21_D = dip[( ( (70 < strike) & (strike < 250)) & ((-180 < rake ) & (rake <-100)) )]
        # G21_R = rake[( ( (70 < strike) & (strike < 250)) & ((-180 < rake ) & (rake <-100)) )]
        #
        #
        # G22_S = strike[( ( (70 < strike) & (strike < 250)) & ((-100 < rake ) & (rake <0)) )]
        # G22_D = dip[( ( (70 < strike) & (strike < 250)) & ((-100 < rake ) & (rake <0)) )]
        # G22_R = rake[( ( (70 < strike) & (strike < 250)) & ((-100 < rake ) & (rake <0)) )]
        #
        # G3_S = strike[(250< strike) & (strike < 349) ]
        # G3_D = dip[(250< strike) & (strike < 349) ]
        # G3_R = rake[(250< strike) & (strike < 349) ]

        ####PZ+ST
        # G1_S = strike[(strike < 250) ]
        # G1_D = dip[(strike < 250) ]
        # G1_R = rake[(strike < 250) ]
        #
        # G2_S = strike[(strike > 250) ]
        # G2_D = dip[(strike > 250) ]
        # G2_R = rake[(strike > 250) ]


        ####TOTAL
        # G1_S = np.hstack(( strike[(strike < 100) ] , strike[(strike>255) ] ))
        # G1_D = np.hstack(( dip[(strike < 100)  ], dip[(strike>255)] ))
        # G1_R = np.hstack(( rake[(strike < 100) ],  rake[(strike>255)] ))
        #
        # G2_S = strike[((strike > 100)  & (strike<255) )]
        # G2_D = dip[((strike > 100)  & (strike<255) )]
        # G2_R = rake[((strike > 100)  & (strike<255) )]


        ## DWTHot_Shift_1
        # PZ = (df['p_z'].values)
        # PR = (df['p_r'].values) * 0.1
        # SZ = (df['s_z'].values) * 0.055
        # SR = (df['s_r'].values) * 0.21
        # ST = (df['s_t'].values) * 0.2
        # AMP =(( (np.abs(np.log10(df['Npz'].values ) - np.log10(df['Nst'].values ))) / np.log(2)) ) ** (1/4)
        #
        # Misfit = PZ + ST + PR + SZ + SR #+ AMP
        #
        # # indices = np.where(Misfit < 0.13)
        # indices = np.where(Misfit < 5.95)
        #
        # misfits = Misfit[indices]
        #
        # strike = df['Strike'].values[indices]
        # dip = df['Dip'].values[indices]
        # rake = df['Rake'].values[indices]

        ###PZ
        # G1_S = strike[(rake < -100) ]
        # G1_D =  dip[(rake < -100) ]
        # G1_R = rake[(rake < -100) ]
        # M1 =  misfits[(rake < -100) ]
        #
        # G2_S =  strike[((rake > -100) & (rake < 0))]
        # G2_D =  dip[((rake > -100) & (rake < 0))]
        # G2_R =  rake[((rake > -100) & (rake < 0))]
        # M2   =  misfits[((rake > -100) & (rake < 0))]
        #
        # G3_S =  strike[((rake > 0) & (rake < 100))]
        # G3_D =  dip[((rake > 0) & (rake < 100))]
        # G3_R =  rake[((rake > 0) & (rake < 100))]
        # M3   =  misfits[((rake > 0) & (rake < 100))]
        #
        # G4_S = strike[(rake > 100) ]
        # G4_D = dip[(rake > 100) ]
        # G4_R = rake[(rake > 100) ]
        # M4   = misfits[(rake > 100) ]

        ## TOTAL

        # G1_S = strike[(strike > 225) ]
        # G1_D =  dip[(strike > 225) ]
        # G1_R = rake[(strike > 225) ]
        # M1 =  misfits[(strike > 225) ]
        #
        # G2_S =  strike[(strike < 225) ]
        # G2_D =  dip[(strike < 225) ]
        # G2_R =  rake[(strike < 225) ]
        # M2   =  misfits[(strike < 225) ]


        ## EH45Tcold
        PZ = (df['p_z'].values)
        PR = (df['p_r'].values) * 0.5
        SZ = (df['s_z'].values) * 0.16
        SR = (df['s_r'].values) * 0.15
        ST = (df['s_t'].values) * 0.94
        AMP =(( (np.abs(np.log10(df['Npz'].values ) - np.log10(df['Nst'].values ))) / np.log(2)) ) ** (1/4) * 1.3

        Misfit = PZ + ST + PR + SZ + SR #+ AMP

        # indices = np.where(Misfit < 0.5)
        # indices = np.where(Misfit < 0.7)
        indices = np.where(Misfit < 8.75)


        misfits = Misfit[indices]

        strike = df['Strike'].values[indices]
        dip = df['Dip'].values[indices]
        rake = df['Rake'].values[indices]

        ## PZ
        # G1_S =  strike[(strike < 160) ]
        # G1_D =  dip[(strike < 160) ]
        # G1_R =  rake[(strike < 160) ]
        # M1   =  misfits[(strike < 160) ]
        #
        # G2_S =  strike[(strike > 160) ]
        # G2_D =  dip[(strike > 160) ]
        # G2_R =  rake[(strike > 160) ]
        # M2   =  misfits[(strike > 160) ]

        ## PZST

        # G1_S =  strike[(strike < 200) ]
        # G1_D =  dip[(strike < 200) ]
        # G1_R =  rake[(strike < 200) ]
        # M1   =  misfits[(strike < 200) ]
        #
        # G2_S =  strike[(strike > 200) ]
        # G2_D =  dip[(strike > 200) ]
        # G2_R =  rake[(strike > 200) ]
        # M2   =  misfits[(strike > 200) ]


        ## TOTAL


        G1_S =  strike[(strike < 125) ]
        G1_D =  dip[(strike < 125) ]
        G1_R =  rake[(strike < 125) ]
        M1   =  misfits[(strike < 125) ]

        G2_S =  strike[(strike > 125) ]
        G2_D =  dip[(strike > 125) ]
        G2_R =  rake[(strike > 125) ]
        M2   =  misfits[(strike > 125) ]



        strike = G2_S
        dip = G2_D
        rake = G2_R

        with open(filepath) as f:
            content = f.readlines()
        azimuth = float(content[54].strip('\n').split(':')[-1])

        veloc_model = '/home/nienke/Documents/Master/Data/Database/TAYAK.npz'

        from obspy.taup import TauPyModel
        model = TauPyModel(model=veloc_model)
        inc_angles = []
        azimuths = []
        phase_names = []
        for i in range(0,1):
            tt_P = model.get_travel_times(source_depth_in_km=depth[i] / 1000, distance_in_degree=epi[i], phase_list=['P'])
            inc_angles.append(tt_P[0].takeoff_angle)
            azimuths.append(azimuth)
            phase_names.append('P')
            tt_S = model.get_travel_times(source_depth_in_km=depth[i] / 1000, distance_in_degree=epi[i], phase_list=['S'])
            inc_angles.append(tt_S[0].takeoff_angle)
            azimuths.append(azimuth)
            phase_names.append('S')




        fig = self.plot_BBB(strike, dip, rake,
                            azimuths=azimuths, inc_angles=inc_angles,
                            phase_names=phase_names,color = 'orange')
        plt.savefig(dir + '/BBB.pdf')




        # fig_bb, ax_bb = plt.subplots(1, 1, figsize=(4, 4))
        #
        # ax_bb.set_xticks([])
        # ax_bb.set_yticks([])
        # ax_bb.axis('off')
        # img = None
        # buf = io.BytesIO()
        #
        # n_lowest = 10
        # lowest_indices = df['Total_misfit'].values.argsort()[0:n_lowest]
        #
        # for v,i in enumerate(lowest_indices):
        # # for i in range(0, len(strike)):
        #
        #     b = beach(fm=[strike[i], dip[i], rake[i]],
        #               width=200, linewidth=0, facecolor='b',
        #               xy=(0, 0), axes=ax_bb, alpha=1, zorder=i)
        #     ax_bb.add_collection(b)
        #     ax_bb.set_xlim((-0.1, 0.1))
        #     ax_bb.set_ylim((-0.1, 0.1))
        #
        #     p = Circle((0., 0,), 0.065, linewidth=2, edgecolor='k', facecolor='b', zorder=5)
        #     ax_bb.add_patch(p)
        #
        #     buf.seek(0)
        #     fig_bb.savefig(buf, format='png')
        #     buf.seek(0)
        #     if img is None:
        #         img = mpimg.imread(buf)
        #     else:
        #         img += mpimg.imread(buf)
        # plt.close(fig_bb)
        # fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        # ax.imshow(img / np.max(img.flatten()))
        # ax.set_xticks([])
        # ax.set_yticks([])
        # ax.axis('off')
        # plt.savefig(dir + '/BBB.pdf')
        ## plt.show()

    def Cluster(self, filepath, savename, directory, skiprows, column_names, burnin, real_v = None):
        dir = directory + '/%s' % (savename.strip('.yaml'))
        if not os.path.exists(dir):
            os.makedirs(dir)

        data = np.loadtxt(filepath, delimiter=',', skiprows=skiprows)

        df = pd.DataFrame(data, columns=column_names)
        n_lowest = 10000
        lowest_indices = df['Total_misfit'].values.argsort()[0:n_lowest]

        PZ = (df['p_z'].values)
        PR = (df['p_r'].values / 0.1) * 0.017
        SZ = (df['s_z'].values / 0.14) * 0.019
        SR = (df['s_r'].values / 0.35) * 0.0556
        ST = (df['s_t'].values) * 2.33

        misfits = PZ +ST +PR +SZ+SR
        strike = df['Strike'].values
        dip = df['Dip'].values
        rake = df['Rake'].values


        # df = pd.read_csv('https://raw.githubusercontent.com/plotly/datasets/master/alpha_shape.csv')
        # df.head()
        #
        # scatter = dict(
        #     mode="markers",
        #     name="y",
        #     type="scatter3d",
        #     x=df['x'], y=df['y'], z=df['z'],
        #     marker=dict(size=2, color="rgb(23, 190, 207)")
        # )
        # clusters = dict(
        #     alphahull=7,
        #     name="y",
        #     opacity=0.1,
        #     type="mesh3d",
        #     x=df['x'], y=df['y'], z=df['z']
        # )
        # layout = dict(
        #     title='3d point clustering',
        #     scene=dict(
        #         xaxis=dict(zeroline=False),
        #         yaxis=dict(zeroline=False),
        #         zaxis=dict(zeroline=False),
        #     )
        # )
        # fig = dict(data=[scatter, clusters], layout=layout)
        # Use py.iplot() for IPython notebook
        # py.iplot(fig, filename='3d point clustering')




        fig = plt.figure(figsize=(25, 6))
        ax1 = plt.subplot(131)

        plt.scatter(strike, dip, c=misfits,cmap = 'rainbow')
        # plt.scatter(strike, dip, c=misfits, s=10000.0 / rake)
        # plt.text(260, 5, "Size = Rake", fontsize=25, bbox=dict(boxstyle="round",
        #                                                        ec=(1., 0.5, 0.5),
        #                                                        fc=(1., 0.8, 0.8),
        #                                                        ))
        plt.xlabel("Strike")
        plt.ylabel("Dip")
        plt.tight_layout()
        # cbar = plt.colorbar()
        # cbar.set_label('Misfit', rotation=270, labelpad=20)
        ax2 = plt.subplot(132)
        plt.scatter(strike, rake, c=misfits,cmap = 'rainbow')
        plt.xlabel("Strike")
        plt.ylabel("Rake")
        # cbar = plt.colorbar()
        # cbar.set_label('Misfit', rotation=270, labelpad=20)
        ax3 = plt.subplot(133)
        plt.scatter(dip, rake, c=misfits,cmap = 'rainbow')
        cbar = plt.colorbar()
        cbar.set_label('Misfit', rotation=270, labelpad=20)
        plt.xlabel("Dip")
        plt.ylabel("Rake")

        plt.tight_layout()
        plt.savefig(dir + '/Clustered.pdf')

        fig = plt.figure(figsize=(15,15))
        ax1 = plt.subplot(111)
        plt.scatter(strike, rake, c=misfits, s=10000.0 / dip,cmap = 'rainbow')
        # plt.text(260, 5, "Size = Rake", fontsize=25, bbox=dict(boxstyle="round",
        #                                                        ec=(1., 0.5, 0.5),
        #                                                        fc=(1., 0.8, 0.8),
        #                                                        ))
        plt.xlabel("Strike")
        plt.ylabel("Rake")
        plt.tight_layout()
        plt.savefig(dir + '/strike_dip_rake_misfit.pdf')

    def plot_BBB(self, strikes, dips, rakes, azimuths=None, inc_angles=None,
                 phase_names=None,
                 color='blue'):
        fig_bb = plt.figure(figsize=(5, 5), dpi=200)
        ax_bb = fig_bb.add_axes([0., 0., 1., 1.])

        ax_bb.set_xticks([])
        ax_bb.set_yticks([])
        ax_bb.axis('off')
        img = None
        buf = io.BytesIO()
        i = 0
        for strike, dip, rake in zip(strikes, dips, rakes):
            i += 1
            b = beach(fm=[strike, dip, rake],
                      width=990, linewidth=0, facecolor=color,
                      xy=(0, 0), axes=ax_bb, alpha=1, zorder=i)
            ax_bb.add_collection(b)
            ax_bb.set_xlim((-1, 1))
            ax_bb.set_ylim((-1, 1))

            buf.seek(0)
            fig_bb.savefig(buf, format='png', dpi=200)
            buf.seek(0)
            if img is None:
                img = mpimg.imread(buf)
            else:
                img += mpimg.imread(buf)

        plt.close(fig_bb)
        # fig, ax = plt.subplots(1, 1, figsize=(5, 5), dpi=250)
        fig = plt.figure(figsize=(5, 5), dpi=200)
        ax = fig.add_axes([0., 0., 1., 1.], label='BBB')
        ax.imshow(img / np.max(img.flatten()))
        ax_2 = fig.add_axes([0., 0., 1., 1.], label='Circle_ray')
        ax_2.set_xlim((-1, 1))
        ax_2.set_ylim((-1, 1))
        p = Circle((0., 0,), 0.99, linewidth=2, edgecolor='k',
                   zorder=0, fill=False)
        ax_2.add_patch(p)
        if azimuths is not None and inc_angles is not None:
            for a, i, phase in zip(azimuths, inc_angles, phase_names):
                x = np.sin(np.deg2rad(a)) * i / 90.
                y = np.cos(np.deg2rad(a)) * i / 90.
                p = Circle((x, y), 0.015, linewidth=2, edgecolor='k',
                           zorder=0, facecolor='k', fill=True)
                ax_2.add_patch(p)
                ax_2.text(x + 0.01, y + 0.01, s=phase, fontsize=24)

        # p = Circle((0., 0,), 0.99, linewidth=2, edgecolor='k',
        #           zorder=5)
        # ax_2.add_patch(p)

        for a in [ax, ax_2]:
            a.set_xticks([])
            a.set_yticks([])
            a.axis('off')
        # plt.tight_layout()
        # plt.savefig('result_%5.1f_km.pdf' % depth)
        return fig

    def event_plot(self, savename,directory,la_receiver, lo_receiver, la_source, lo_source):
        dir = directory + '/%s' % (savename.strip('.yaml'))
        la_r = la_receiver
        lo_r = lo_receiver
        la_s = la_source
        lo_s = lo_source


        mars_dir = '/home/nienke/Documents/Master/Data/mars_pictures/Mars_lightgray.jpg'

        fig = plt.figure()

        m= Basemap(projection='moll',lon_0=round(0.0))

        # draw parallels and meridians.
        par = np.arange(-90, 90, 30)
        label_par = np.full(len(par), True, dtype=bool)
        meridians = np.arange(-180, 180, 30)
        label_meri = np.full(len(meridians), True, dtype=bool)

        m.drawmeridians(np.arange(-180, 180, 30),labels=label_meri)
        m.drawparallels(np.arange(-90, 90, 30),label=label_par)


        m.warpimage(mars_dir)
        mstatlon, mstatlat = m(lo_r, la_r)
        m.plot(mstatlon, mstatlat, 'k^', markersize=8)

        EQlon, EQlat = m(lo_s, la_s)
        m.plot(EQlon, EQlat, 'r*', markersize=10, zorder=10)
        plt.savefig(dir + '/Location_Event.pdf')

    def write_to_dict(self, list_of_parameters):
        parameters = {
            'velocity_model': list_of_parameters[0].strip('\n'),
            'MCMC': list_of_parameters[1].strip('\r\n'),
            'misfit': list_of_parameters[2].strip('\r\n'),
            'noise': eval(list_of_parameters[3]),
            'sdr': eval(list_of_parameters[4].strip('\n')),
            'plot_modus': eval(list_of_parameters[5].strip('\n')),
            'alpha': np.float(list_of_parameters[6]),
            'beta': np.float(list_of_parameters[7]),
            'azimuth': np.float(list_of_parameters[8]),
            'components': np.asarray(list_of_parameters[9].strip('\r\n').split(',')),
            'la_r': np.float(list_of_parameters[10]),
            'lo_r': np.float(list_of_parameters[11]),
            'filter': list_of_parameters[12].strip('\r\n'),
            'definition': list_of_parameters[13].strip('\r\n'),
            'kind': list_of_parameters[14].strip('\r\n'),
            'network': list_of_parameters[15].strip('\r\n'),
            'sample_number': np.int(list_of_parameters[16]),
            'var_est': np.float(list_of_parameters[17]),
            'epi_range_max': np.int(list_of_parameters[18]),
            'epi_range_min': np.int(list_of_parameters[19]),
            'epi_spread': np.int(list_of_parameters[20]),
            'depth_range_max': np.int(list_of_parameters[21]),
            'depth_range_min': np.int(list_of_parameters[22]),
            'depth_spread': np.int(list_of_parameters[23]),
            'strike_range_max': np.int(list_of_parameters[24]),
            'strike_range_min': np.int(list_of_parameters[25]),
            'strike_spread': np.int(list_of_parameters[26]),
            'dip_range_max': np.int(list_of_parameters[27]),
            'dip_range_min': np.int(list_of_parameters[28]),
            'dip_spread': np.int(list_of_parameters[29]),
            'rake_range_max': np.int(list_of_parameters[30]),
            'rake_range_min': np.int(list_of_parameters[31]),
            'rake_spread': np.int(list_of_parameters[32]),
            'time_at_rec': np.asarray(list_of_parameters[33].strip('\r\n').split(','), dtype=int)}
        return parameters

    def shift(self, np_array, time_shift):
        new_array = np.zeros_like(np_array)
        if time_shift < 0:
            new_array[-time_shift:] = np_array[:time_shift]
        elif time_shift == 0:
            new_array[:] = np_array[:]
        else:
            new_array[:-time_shift] = np_array[time_shift:]
        return new_array

    def marginal_grid(
            self, savename, directory,
            samples: np.ndarray,
            dimensions_list,
            bins: int = 25,
            show: bool = True,
            colormap_2d=plt.get_cmap("Greys"),
            color_1d="black"
    ):
        number_of_plots = len(dimensions_list)
        import matplotlib.gridspec as _gridspec
        dir = directory + '/%s' % (savename.strip('.yaml'))
        plt.figure(figsize=(12, 12))
        gs1 = _gridspec.GridSpec(number_of_plots, number_of_plots)
        gs1.update(wspace=0.05, hspace=0.05)  # set the spacing between axes.

        # Get extent of every set
        dim_range = []
        for i_dim in range(number_of_plots):
            min = samples[:,dimensions_list[i_dim]].min()
            max = samples[:,dimensions_list[i_dim]].max()
            dim_range.append((min, max))

        for i_plot in range(number_of_plots):
            # print(i_plot, i_plot) # grid indices for diagonal
            axis = plt.subplot(gs1[i_plot + (number_of_plots) * i_plot])

            # Modify axes
            if i_plot != number_of_plots - 1:
                axis.set_xticklabels([])
                axis.tick_params(axis="x", which="both", bottom=False, top=False)
            else:
                labels = ['PZ','PR','SZ','SR','ST','Total']
                # axis.set_xlabel(f"dimension {dimensions_list[number_of_plots - 1]}")
                axis.set_xlabel(labels[number_of_plots - 1])

            axis.set_yticklabels([])
            axis.tick_params(axis="y", which="both", left=False, right=False)
            if i_plot == 0:
                axis.set_ylabel("PZ")

            # Plot histogram on diagonal
            axis.plot(
                samples[:,dimensions_list[i_plot]]
            )

            for j_plot in range(i_plot):
                # print(i_plot, j_plot) # grid indices for lower left
                axis = plt.subplot(gs1[j_plot + (number_of_plots) * i_plot])

                # Modify axes
                if i_plot != number_of_plots - 1:
                    axis.set_xticklabels([])
                    axis.tick_params(axis="x", which="both", bottom=False, top=False)
                else:
                    labels = ['PZ', 'PR', 'SZ', 'SR', 'ST', 'Total']
                    axis.set_xlabel(labels[j_plot])
                    # axis.set_xlabel(f"dimension {dimensions_list[j_plot]}")

                if j_plot != 0:
                    axis.set_yticklabels([])
                    axis.tick_params(axis="y", which="both", left=False, right=False)
                else:
                    # axis.set_ylabel(f"dimension {dimensions_list[i_plot]}")
                    labels = ['PZ', 'PR', 'SZ', 'SR', 'ST', 'Total']
                    axis.set_ylabel(labels[i_plot])

                # Plot 2d marginals
                # axis.hist2d(
                #     samples[:,dimensions_list[j_plot]],
                #     samples[:,dimensions_list[i_plot]],
                #     bins=bins,
                #     range=[dim_range[j_plot], dim_range[i_plot]],
                #     cmap=colormap_2d,
                # )

                axis.scatter(
                    samples[:,dimensions_list[j_plot]],
                    samples[:,dimensions_list[i_plot]],
                )

                # print(i_plot, j_plot) # grid indices for lower left
                axis = plt.subplot(gs1[i_plot + (number_of_plots) * j_plot])

                axis.set_xticklabels([])
                axis.tick_params(axis="x", which="both", bottom=False, top=False)
                axis.set_yticklabels([])
                axis.tick_params(axis="y", which="both", left=False, right=False)

                axis.axis("off")

                # correlation = np.corrcoef(
                #     samples[:,dimensions_list[j_plot]], samples[:,dimensions_list[i_plot]]
                # )[1][0]
                # axis.text(
                #     0.5,
                #     0.5,
                #     f"{correlation:.2f}",
                #     horizontalalignment="center",
                #     verticalalignment="center",
                #     fontsize=40 * np.abs(correlation),
                #     transform=axis.transAxes,
                # )

        if show:
            plt.show()
        else:
            plt.savefig(dir + '/Marginals.pdf')



if __name__ == '__main__':
    main()
