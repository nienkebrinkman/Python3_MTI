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
    ## Post - Processing [processing the results from inversion]
    result = Post_processing_sdr()

    strike, dip, rake = aux_plane(238, 80, 143)

    directory = '/home/nienke/Documents/Master/Data/Mars/S0235b/waveforms/Output/'
    path_to_file = '/home/nienke/Documents/Master/Data/Mars/S0235b/waveforms/Output/TAYAK_NEW_5.txt'
    path_to_file_BBB = '/home/nienke/Documents/Master/Data/Mars/S0235b/waveforms/Output/TAYAK_BBB_NEW_5.txt'
    # path_to_file ='/home/nienke/Documents/Applied_geophysics/Thesis/anaconda/Earth/python3/Events/new-trials/7_BBB.txt'

    # par = Get_Parameters()
    # REAL = par.get_unkown()
    # PRIOR = par.get_prior()
    # strike, dip, rake = aux_plane( REAL['strike'], REAL['dip'], REAL['rake'])
    # real_v = np.array([REAL['epi'], REAL['depth_s'],strike,dip,rake,1e+16])
    # real_v = np.array([35.2068855191, 16848.1405882, 99.88000245897199, 77.02994881296266, 68.80551508147109,3.16227766016798e+16])
    # real_v=np.array([90.9362911326,17961.4407815,268.68760185363493, 89.27813406577755,59.759221608124186,2818382931267215.0]) # 4.3 Event
    # real_v=np.array([64.8000651008,14567.2191587,189.90144092357255,77.55305609924525,79.99955781783979,7943282347240532.0])
    # real_v=np.array([80.32,35000,333,61,83,1.016e+19]) # EARTH Japan 6.6
    # real_v=np.array([80.25,7700,186,36,69,4.965e+17]) # Earth Vigrinia 5.8
    # real_v=np.array([49.62,1000,211,85,-171,7.300e+17 ]) # Earth Eastern Gulf Of Aden 6.0
    # real_v=np.array([88.52638906313514,10000,28,22,-107,4.4668359215096166e+17 ]) # -- 2 --
    # real_v=np.array([93.30019067657031,10000,87,44,-96,1.7782794100389228e+18 ]) # -- 3 --
    # real_v=np.array([54.01252893905707,10000,261,85,-175,1.2589254117941714e+18 ]) # -- 4 --
    # real_v=np.array([54.01252893905707,10000,262,83,-175,6.309573444801892e+17 ]) # -- 5 --
    # real_v=np.array([95.18311630916064,18000,150,17,-76,3.548133892335731e+18]) # -- 7 --

    magnitude = 4.46  # Magnitude of the earthquake
    exp = mtm.magnitude_to_moment(magnitude)  # convert the mag to moment in [Nm] # 5495408738576270.0
    m_pp = -214996686134000.0
    m_rp = 247303763622000.0
    m_rr = 3531870796970000.0
    m_rt = -3487091494290000.0
    m_tp = 1008979372750000.0
    m_tt = -3316874110840000.0

    m_dd = m_rr
    m_nn = m_tt
    m_ee = m_pp
    m_nd = m_rt
    m_ed = -m_rp
    m_ne = -m_tp
    # init pyrocko moment tensor
    m = mtm.MomentTensor(
        mnn=m_nn * exp,
        mee=m_ee * exp,
        mdd=m_dd * exp,
        mne=m_ne * exp,
        mnd=m_nd * exp,
        med=m_ed * exp)

    # gives out both nodal planes:
    (s1, d1, r1), (s2, d2, r2) = m.both_strike_dip_rake()

    print('strike1=%s, dip1=%s, rake1=%s' % (
        s1, d1, r1))  # strike1=244.6022883548614, dip1=73.52631851233593, rake1=67.50622085517007
    print('strike2=%s, dip2=%s, rake2=%s' % (
        s2, d2, r2))  # strike2=120.19814045454422, dip2=27.62588224193955, rake2=142.29811781206305
    # real_v = np.array([88.4756, 38438, s2, d2, r2, exp])  # MSS event 5.0
    # beachball([s1,d1,r1], size=200, linewidth=2, facecolor='b')

    real_v = None # Real Marsquakes

    savename = 'Trials'
    show = False  # Choose True for direct show, choose False for saving
    skiprows = 26  # 67
    # column_names= ["Epi", "Depth", "Strike", "Dip", "Rake", "M0","Total_misfit","S_z","S_r","S_t","P_z","P_r","BW_misfit","Rtot","Ltot"]
    column_names = ["Epi", "Depth", "Strike", "Dip", "Rake", "M0", "Total_misfit", "p_z", "p_r", "s_z", "s_r", "s_t",
                    'bw_tot', 'Shift_S', 'Shift_P', 'accept']
    # column_names= ["Epi", "Depth", "Strike", "Dip", "Rake","M0","Total_misfit"]
    # path_to_file, save = result.convert_txt_folder_to_yaml(path, savename)
    # result.get_stereonets(filepath=path_to_file, savename=savename, directory=directory, show=show)

    # result.plot_streams(stream_filepath=path_to_stream,filepath=path_to_file,savename=savename, directory=directory,skiprows=skiprows ,column_names=column_names)
    burnin = 2000

    File = np.loadtxt(path_to_file, delimiter=',', skiprows=skiprows)
    strike, dip, rake = aux_plane(np.mean(File[:,2]), np.mean(File[:,3]), np.mean(File[:,4]))
    # real_v = np.array([None,None, strike, dip, rake, None])
    real_v = None

    # result.trace(filepath=path_to_file, savename=savename, directory=directory, skiprows=skiprows,
    #              column_names=column_names,burnin=burnin ,real_v=real_v)
    # result.get_BBB(filepath=path_to_file_BBB, savename=savename, directory=directory, skiprows=skiprows,
    #                column_names=column_names,  burnin=burnin,real_v=real_v)
    #
    # result.marginal_grid( savename=savename, directory=directory,samples=File[:,:], dimensions_list = [0,1,2,3,4,5], show = False)
    result.get_convergence(filepath=path_to_file, savename = savename, directory = directory, skiprows = skiprows, column_names = column_names, show=False)

    # result.event_plot(savename = savename,directory = directory,la_receiver = 4.5, lo_receiver = 136 , la_source = 10.99, lo_source = 160.95)


class Post_processing_sdr:
    def convert_txt_file_to_yaml(self, filepath):
        save_path = filepath.replace('.txt', '.yaml')
        if os.path.isfile(save_path) == True:
            print("File is already converted to .yaml file, located in \n %s" % save_path)
            return save_path, save_path
        with open(filepath) as infile:
            dat = np.genfromtxt(infile, delimiter=',', skip_header=34, skip_footer=1)
            # par = open(fname, "r").readlines()[:34]
            # parameters = self.write_to_dict(par)
        # data_file = {'data': dat, 'parameters': parameters}
        data_file = {'data': dat}
        with open(save_path, 'w') as yaml_file:
            yaml.dump(data_file, yaml_file, default_flow_style=False)
        yaml_file.close()
        savename = save_path.split('/')[-1].strip('.yaml')
        return save_path, savename

    def convert_txt_folder_to_yaml(self, dir_of_txt_files, savename):
        filenames = glob.glob("%s/*.txt" % dir_of_txt_files)
        save_path = dir_of_txt_files + '/%s.yaml' % savename
        if os.path.isfile(save_path) == True:
            return save_path, savename
        for i, fname in enumerate(filenames):

            with open(fname) as infile:
                if i == 0:
                    # dat = np.genfromtxt(infile, delimiter=',',skip_footer=1)
                    dat = np.genfromtxt(infile, delimiter=',', skip_header=35, skip_footer=1)
                    # par = open(fname, "r").readlines()[:34]
                    # parameters = self.write_to_dict(par)

                else:
                    dat = np.vstack((dat, np.genfromtxt(infile, delimiter=',', skip_header=35, skip_footer=1)))
                    # dat = np.vstack((dat, np.genfromtxt(infile, delimiter=',',skip_footer=1)))
            # data_file = {'data': dat, 'parameters': parameters}
            data_file = {'data': dat}
            with open(save_path, 'w') as yaml_file:
                yaml.dump(data_file, yaml_file, default_flow_style=False)
            yaml_file.close()
            return save_path, savename

    def get_accepted_samples(self, filepath, savename, directory, skiprows, column_names):
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
        length = len(data[0]) - (len(column_names) + 3)
        R_length = int(data[0][-2])
        L_length = int(data[0][-1])

        for i in range(R_length):
            column_names = np.append(column_names, "R_%i" % (i + 1))
        for i in range(L_length):
            column_names = np.append(column_names, "L_%i" % (i + 1))
        column_names = np.append(column_names, "Accepted")
        column_names = np.append(column_names, "Rayleigh_length")
        column_names = np.append(column_names, "Love_length")

        df = pd.DataFrame(data,
                          columns=column_names)
        total = len(df["Accepted"])
        accepted = np.sum(df["Accepted"] == 1)
        savepath = dir + '/Accept_%i_outof_%i.txt' % (accepted, total)
        with open(savepath, 'w') as save_file:
            save_file.write("%i,%i\n\r" % (total, accepted))
        save_file.close()
        print("Total amount of samples = %i" % len(df["Accepted"]))
        print("Total amount of accepted samples = %i" % np.sum(df["Accepted"] == 1))

    def get_convergence(self, filepath, savename, directory, skiprows, column_names, show=True):
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
        plt.figure(1)
        ax = plt.subplot(111)
        ax.plot(np.arange(0, len(df['Total_misfit'])), df['Total_misfit'])
        # ax.plot(np.arange(0, len(df['Total_misfit'])), df['p_z'],c='r')
        # ax.plot(np.arange(0, len(df['Total_misfit'])), df['s_t'],c='g')
        plt.yscale('log')
        plt.xlabel('Iteration')
        plt.ylabel('-Log(likelihood)')
        # ax.invert_yaxis()
        ax.xaxis.tick_top()

        plt.tight_layout()
        if show == True:
            plt.show()
        else:
            plt.savefig(dir + '/Convergence.pdf')
            plt.close()


        amount_of_samples = len(df['Total_misfit'])
        loop_length = int(amount_of_samples/10 -1)
        epi_mean    = np.zeros(loop_length)
        depth_mean  = np.zeros(loop_length)
        strike_mean = np.zeros(loop_length)
        dip_mean    = np.zeros(loop_length)
        rake_mean   = np.zeros(loop_length)
        M0_mean     = np.zeros(loop_length)


        for i in range(loop_length):
            epi_mean[i] = np.mean(df['Epi'].values[0:int(i * 10 + 1)])
            depth_mean[i] = np.mean(df['Depth'].values[0:int(i * 10 + 1)])
            strike_mean[i] = np.mean(df['Strike'].values[0:int(i * 10 + 1)])
            dip_mean[i] = np.mean(df['Dip'].values[0:int(i * 10 + 1)])
            rake_mean[i] = np.mean(df['Rake'].values[0:int(i * 10 + 1)])
            M0_mean[i] = np.mean(df['M0'].values[0:int(i * 10 + 1)])

        plt.close()
        plt.figure(figsize=(20,10))
        ax1 = plt.subplot(231)
        ax1.plot(np.arange(loop_length) * 10, epi_mean , label = 'Epicentral distance')
        ax1.set_ylabel("Mean", fontsize=25)
        ax1.set_xlabel("Sample", fontsize=25)
        ax1.tick_params(axis='x', labelsize=20)
        ax1.tick_params(axis='y', labelsize=20)
        ax1.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        plt.legend(fontsize=20)
        plt.tight_layout()
        # plt.yscale('log')


        ax2 = plt.subplot(232)
        ax2.plot(np.arange(loop_length) * 10, depth_mean , label = 'Depth')
        ax2.set_ylabel("Mean", fontsize=25)
        ax2.set_xlabel("Sample", fontsize=25)
        ax2.tick_params(axis='x', labelsize=20)
        ax2.tick_params(axis='y', labelsize=20)
        ax2.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        plt.legend(fontsize=20)
        plt.tight_layout()
        # plt.yscale('log')


        ax3 = plt.subplot(233)
        ax3.plot(np.arange(loop_length) * 10, M0_mean , label = 'M0')
        ax3.set_ylabel("Mean", fontsize=25)
        ax3.set_xlabel("Sample", fontsize=25)
        ax3.tick_params(axis='x', labelsize=20)
        ax3.tick_params(axis='y', labelsize=20)
        ax3.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        plt.legend(fontsize=20)
        plt.tight_layout()
        # plt.yscale('log')

        ax4 = plt.subplot(234)
        ax4.plot(np.arange(loop_length) * 10, strike_mean , label = 'Strike')
        ax4.set_ylabel("Mean", fontsize=25)
        ax4.set_xlabel("Sample", fontsize=25)
        ax4.tick_params(axis='x', labelsize=20)
        ax4.tick_params(axis='y', labelsize=20)
        ax4.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        plt.legend(fontsize=20)
        plt.tight_layout()
        # plt.yscale('log')

        ax5 = plt.subplot(235)
        ax5.plot(np.arange(loop_length) * 10, dip_mean , label = 'Dip')
        ax5.set_ylabel("Mean", fontsize=25)
        ax5.set_xlabel("Sample", fontsize=25)
        ax5.tick_params(axis='x', labelsize=20)
        ax5.tick_params(axis='y', labelsize=20)
        ax5.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        plt.legend(fontsize=20)
        plt.tight_layout()
        # plt.yscale('log')

        ax6 = plt.subplot(236)
        ax6.plot(np.arange(loop_length) * 10, rake_mean , label = 'Rake')
        ax6.set_ylabel("Mean", fontsize=25)
        ax6.set_xlabel("Sample", fontsize=25)
        ax6.tick_params(axis='x', labelsize=20)
        ax6.tick_params(axis='y', labelsize=20)
        ax6.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        plt.legend(fontsize=20)
        plt.tight_layout()
        # plt.yscale('log')


        if show == True:

            plt.show()
        else:
            plt.savefig(dir + '/Parameter_convergence.pdf')
            plt.close()


        plt.figure(figsize=(20,10))
        ax1 = plt.subplot(511)
        ax1.plot(df['p_z'], label = 'Epicentral distance')
        ax1.set_ylabel("Mean", fontsize=25)
        ax1.set_xlabel("Sample", fontsize=25)
        ax1.tick_params(axis='x', labelsize=20)
        ax1.tick_params(axis='y', labelsize=20)
        ax1.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        plt.legend(fontsize=20)
        plt.tight_layout()
        # plt.yscale('log')


        ax2 = plt.subplot(512)
        ax2.plot( df['p_r'] , label = 'Depth')
        ax2.set_ylabel("Mean", fontsize=25)
        ax2.set_xlabel("Sample", fontsize=25)
        ax2.tick_params(axis='x', labelsize=20)
        ax2.tick_params(axis='y', labelsize=20)
        ax2.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        plt.legend(fontsize=20)
        plt.tight_layout()
        # plt.yscale('log')


        ax3 = plt.subplot(513)
        ax3.plot( df['s_z']  , label = 'M0')
        ax3.set_ylabel("Mean", fontsize=25)
        ax3.set_xlabel("Sample", fontsize=25)
        ax3.tick_params(axis='x', labelsize=20)
        ax3.tick_params(axis='y', labelsize=20)
        ax3.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        plt.legend(fontsize=20)
        plt.tight_layout()
        # plt.yscale('log')

        ax4 = plt.subplot(514)
        ax4.plot( df['s_r'] , label = 'Strike')
        ax4.set_ylabel("Mean", fontsize=25)
        ax4.set_xlabel("Sample", fontsize=25)
        ax4.tick_params(axis='x', labelsize=20)
        ax4.tick_params(axis='y', labelsize=20)
        ax4.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        plt.legend(fontsize=20)
        plt.tight_layout()
        # plt.yscale('log')

        ax5 = plt.subplot(515)
        ax5.plot( df['s_t'] , label = 'Dip')
        ax5.set_ylabel("Mean", fontsize=25)
        ax5.set_xlabel("Sample", fontsize=25)
        ax5.tick_params(axis='x', labelsize=20)
        ax5.tick_params(axis='y', labelsize=20)
        ax5.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        plt.legend(fontsize=20)
        plt.tight_layout()
        # plt.yscale('log')
        plt.show()
        if show == True:
            plt.show()
        else:
            plt.savefig(dir + '/misfits.pdf')
            plt.close()
    def combine_all(self, filepath, savename, directory, skiprows, column_names, real_v):
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

        par = Get_Paramters()
        REAL = par.get_unkown()

        df = pd.DataFrame(data,
                          columns=column_names)

        df_select = df[["Epi", "Depth", "Strike", "Dip", "Rake", "M0"]]
        # fig = plt.figure(figsize=(20,10))
        fig = plt.figure()
        ax1 = plt.subplot2grid((4, 1), (0, 0))
        ax1.axhline(y=real_v[0], linewidth=0.3, linestyle=':', color='b')
        ax1.plot(df_select['Epi'], label="Epi", c='b')
        ax1.tick_params("y", colors='b')
        ax1.set_ylabel('Epi[degree]', color='b')
        ax2 = ax1.twinx()
        ax2.plot(df_select['Depth'], label="Depth", c='r')
        ax2.axhline(y=real_v[1], linewidth=0.3, linestyle=':', color='r')
        ax2.tick_params('y', colors='r')
        ax2.set_ylabel('Depth [m]', color='r')
        plt.tight_layout()
        ax3 = plt.subplot2grid((4, 1), (1, 0))
        lns1 = ax3.plot(df_select['Strike'], label="Strike", color="b")
        ax3.axhline(y=real_v[2], linewidth=0.3, linestyle=':', color='b')
        ax3.axhline(y=real_v[3], linewidth=0.3, linestyle=':', color='g')
        ax3.axhline(y=real_v[4], linewidth=0.3, linestyle=':', color='r')
        lns2 = ax3.plot(df_select['Dip'], label="Dip", color='g')
        lns3 = ax3.plot(df_select['Rake'], label="Rake", color='r')
        # ax3.set_ylabel('[degree]', color='k')
        plt.tight_layout()
        axes = ax3.twinx()
        lns4 = axes.plot(df_select['M0'], label="M0", c='m')
        # ax2.axhline(y=REAL['depth_s'], linewidth=0.3, linestyle=':', color='r')
        axes.tick_params('y', colors='m')
        # axes.set_ylabel('M0', color='r')
        plt.yscale('log')
        lns = lns1 + lns2 + lns3 + lns4
        labs = [l.get_label() for l in lns]
        ax3.legend(lns, labs, loc='center left', bbox_to_anchor=(1.1, 0.5))
        plt.tight_layout()
        R_select = df.filter(like='R_')
        L_select = df.filter(like='L_')
        df_select_xi = df[["Total_misfit", "S_z", "S_r", "S_t", "P_z", "P_r", "BW_misfit", "Rtot", "Ltot"]]
        ax4 = plt.subplot2grid((4, 1), (2, 0))
        ax4.plot(df_select_xi['Total_misfit'], label="Total_misfit", c='r')
        # ax4.plot(df_select_xi['S_z'], label = "S_z")
        # ax4.plot(df_select_xi['S_r'], label = "S_r")
        # ax4.plot(df_select_xi['S_t'], label = "S_t")
        # ax4.plot(df_select_xi['P_z'], label = "P_z")
        # ax4.plot(df_select_xi['P_r'], label = "P_r")
        ax4.plot(df_select_xi['BW_misfit'], label="BW_tot")
        plt.yscale('log')
        # plt.ylim((pow(10, 0), pow(10, 3)))
        plt.legend(loc='center left', bbox_to_anchor=(1.1, 0.5))
        plt.tight_layout()
        ax5 = plt.subplot2grid((4, 1), (3, 0))
        ax5.plot(df_select_xi['Rtot'], label="Rtot")
        ax5.plot(df_select_xi['Ltot'], label="Ltot")
        # for i in R_select:
        #     ax5.plot(R_select[i])
        # for i in L_select:
        #     ax5.plot(L_select[i])
        plt.yscale('log')
        # plt.ylim((pow(10, 0), pow(10, 4)))
        # ax6 = ax5.twinx()
        # # for i in L_select:
        # #     ax6.plot(L_select[i])
        # l2 = ax6.plot(df_select_xi['Ltot'], label = "Ltot")
        # plt.yscale('log')
        # lns2 = l1 + l2
        # labs = [l.get_label() for l in lns2]
        # ax5.legend(lns, labs, loc='center left', bbox_to_anchor=(1.1, 0.5))
        plt.legend(loc='center left', bbox_to_anchor=(1.1, 0.5))
        # plt.ylim((pow(10, 0), pow(10, 4)))
        plt.tight_layout()
        # plt.show()
        plt.savefig(dir + '/combined_all_par.pdf')
        plt.close()

    def Misfits(self, filepath, savename, directory, skiprows, column_names):
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

        par = Get_Paramters()
        REAL = par.get_unkown()

        df = pd.DataFrame(data,
                          columns=column_names)

        fig = plt.figure(figsize=(20, 5))
        # fig = plt.figure()
        R_select = df.filter(like='R_')
        L_select = df.filter(like='L_')
        df_select_xi = df[["Total_misfit", "S_z", "S_r", "S_t", "P_z", "P_r", "BW_misfit", "Rtot", "Ltot"]]
        ax4 = plt.subplot2grid((1, 1), (0, 0))
        ax4.plot(df_select_xi['Total_misfit'], label="Total_misfit", c='r', linewidth=1.1)
        plt.legend(fontsize=20)
        ax4.plot(df_select_xi['Total_misfit'][0:200], c='k', linewidth=2)
        ymin, ymax = ax4.get_ylim()
        ax4.vlines(200, ymin=ymin, ymax=ymax, colors='k', linewidth=2, linestyles=':')
        # ax4.plot(df_select_xi['S_z'], label = "S_z")
        # ax4.plot(df_select_xi['S_r'], label = "S_r")
        # ax4.plot(df_select_xi['S_t'], label = "S_t")
        # ax4.plot(df_select_xi['P_z'], label = "P_z")
        # ax4.plot(df_select_xi['P_r'], label = "P_r")
        # ax4.plot(df_select_xi['BW_misfit'], label = "BW_tot",linewidth = 0.3)
        # ax4.plot(df_select_xi['Rtot'], label="Rtot",linewidth = 0.3)
        # ax4.plot(df_select_xi['Ltot'], label="Ltot",linewidth = 0.3)
        ax4.set_xlabel('N = %i' % len(df_select_xi['Total_misfit']), color='b', fontsize=20)
        ax4.set_ylabel('Misfit value', color='k', fontsize=20)
        ax4.tick_params(axis='x', labelsize=18)
        ax4.tick_params(axis='y', labelsize=18)
        plt.yscale('log')

        # plt.xlim((0, 10000))
        # plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
        plt.tight_layout()

        plt.savefig(dir + '/misfits.pdf')
        plt.close()

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

    def trace(self, filepath, savename, directory, skiprows, column_names,  burnin, real_v = None,):
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
        row = 0



        ax1 = plt.subplot2grid((1, 3), (0, 0))
        bin = int(360 * (len(df_select['Strike'][burnin:])**(1/3) / (3.49 * np.std(df_select['Strike'][burnin:]))))
        plt.hist(df_select['Strike'][burnin:], bins= 80, alpha=0.8)

        ymin, ymax = ax1.get_ylim()
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
        bin = int(90 * (len(df_select['Dip'][burnin:]) ** (1 / 3) / (3.49 * np.std(df_select['Dip'][burnin:]))))
        plt.hist(df_select['Dip'][burnin:], bins=80, alpha=0.8)

        ymin, ymax = ax2.get_ylim()
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
        bin = int(360 * (len(df_select['Rake'][burnin:]) ** (1 / 3) / (3.49 * np.std(df_select['Rake'][burnin:]))))
        plt.hist(df_select['Rake'][burnin:], bins=80, alpha=0.8)
        ymin, ymax = ax3.get_ylim()
        if real_v is not None:
            plt.vlines(real_v[4], ymin=ymin, ymax=ymax, colors='g', linewidth=3, label='Auxiliary plane')
            plt.vlines(rake, ymin=ymin, ymax=ymax, colors='k', linewidth=3, label='Fault plane')
            plt.legend(loc='upper left', fontsize=20)
            # plt.legend(loc='center left', bbox_to_anchor=(1.1, 0.5))
        ax3.set_title("Density Rake", color='b', fontsize=25)
        ax3.set_xlabel("N=%i" % (len(df_select['Rake'])), fontsize=25)
        ax3.tick_params(axis='x', labelsize=20)
        ax3.tick_params(axis='y', labelsize=20)
        ax3.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        ax3.set_xlim(-180, 180)
        plt.tight_layout()

        # plt.show()
        plt.savefig(dir + '/Trace_fault.pdf')

        fig = plt.figure(figsize=(25, 6))

        ax4 = plt.subplot2grid((1, 3), (0, 0))
        plt.hist(df_select['Depth'][burnin:] / 1000, bins=50, alpha=0.8)
        ymin, ymax = ax4.get_ylim()
        xmin, xmax = ax4.get_xlim()
        if real_v is not None:
            plt.vlines(real_v[1], ymin=ymin, ymax=ymax, colors='k', linewidth=3, label='True model')

        ax4.set_title("Density Depth", color='b', fontsize=20)
        # ax4.set_xlabel("N=%i" % (len(df_select['Depth'])), fontsize=20)
        ax4.set_xlabel("Depth [km]", fontsize=20)
        ax4.set_ylabel("Posterior marginal", fontsize=20)
        ax4.tick_params(axis='x', labelsize=20)
        ax4.tick_params(axis='y', labelsize=20)
        ax4.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        ax4.ticklabel_format(style="sci", axis='x', scilimits=(-2, 2))

        # ax4.set_xlim(real_v[1]-10000, real_v[1]+10000)
        ax4.set_xlim(xmin, xmax)
        # plt.legend( fontsize=20)
        plt.tight_layout()

        ax5 = plt.subplot2grid((1, 3), (0, 1))
        plt.hist(df_select['Epi'][burnin:], bins=50, alpha=0.8)
        # plt.hist(df['Epi'][burnin:],bins=100,alpha=0.5,color = 'b',label = 'MAAK')
        # plt.hist(df2['Epi'][burnin:],bins=100,alpha=0.5,color = 'r',label = 'EH45TcoldCrust_1b')
        if real_v is not None:
            ymin, ymax = ax5.get_ylim()

        # plt.vlines(df_select['Epi'][0], ymin=ymin, ymax=ymax, colors='r', linewidth=3, label='Start model')
            plt.vlines(real_v[0], ymin=ymin, ymax=ymax, colors='k', linewidth=3,label = 'True model')
        ax5.set_title("Density Epicentral Distance", color='b', fontsize=20)
        ax5.set_xlabel("N=%i" % (len(df_select['Epi'])), fontsize=20)
        ax5.tick_params(axis='x', labelsize=20)
        ax5.tick_params(axis='y', labelsize=20)
        ax5.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        # plt.legend(fontsize=20)
        plt.tight_layout()

        ax6 = plt.subplot2grid((1, 3), (0, 2))
        Mw = 2.0 / 3.0 * (np.log10(df['M0'][burnin:]) - 9.1)
        plt.hist(Mw, bins=np.arange(3, 5, 0.05), alpha=0.8)
        # plt.hist(df['M0'][burnin:],bins=100,alpha=0.5,color = 'b',label = 'MAAK')
        # plt.hist(df2['M0'][burnin:],bins=100,alpha=0.5,color = 'r',label = 'EH45TcoldCrust_1b')
        # if real_v is not None:
        #     ymin, ymax = ax6.get_ylim()
        #     plt.vlines(M0_true, ymin=ymin, ymax=ymax, colors='k', linewidth=3, label= 'True model')
        # plt.vlines(df['M0'][1], ymin=ymin, ymax=ymax, colors='r', linewidth=3,label = 'Start model')
        ax6.set_title("Density Moment Magnitude", color='b', fontsize=20)
        ax6.set_xlabel("N=%i" % (len(df['M0'])), fontsize=18)
        ax6.tick_params(axis='x', labelsize=18)
        ax6.tick_params(axis='y', labelsize=18)
        ax6.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        # ax6.set_xscale('log')
        # ax6.set_xlim(0e18, 0.5e18 )
        # plt.legend(loc='center left', bbox_to_anchor=(1.1, 0.5))
        plt.legend(loc='upper right', fontsize=20)
        plt.tight_layout()

        # plt.show()
        plt.savefig(dir + '/Trace_position.pdf')

        n_lowest = 1000
        # pos = np.argmin(df['Total_misfit'].values)
        lowest_indices = df['Total_misfit'].values.argsort()[0:n_lowest]
        lowest_misfits = df['Total_misfit'].values[lowest_indices]
        lowest_strike = df['Strike'].values[lowest_indices]
        lowest_dip = df['Dip'].values[lowest_indices]
        lowest_rake = df['Rake'].values[lowest_indices]

        fig = plt.figure(figsize=(25, 6))
        row = 0

        ax1 = plt.subplot2grid((1, 3), (0, 0))
        plt.plot(lowest_strike,lowest_misfits,'bo')

        ymin, ymax = ax1.get_ylim()
        if real_v is not None:
            plt.vlines(real_v[2], ymin=ymin, ymax=ymax, colors='g', linewidth=3, label='Auxiliary plane')
            plt.vlines(strike, ymin=ymin, ymax=ymax, colors='k', linewidth=3, label='Fault plane')
        ax1.set_title("Strike", color='b', fontsize=25)
        ax1.set_xlabel("N=%i" % (len(df_select['Strike'])), fontsize=25)
        ax1.set_ylabel("Misfit", fontsize=25)
        ax1.tick_params(axis='x', labelsize=20)
        ax1.tick_params(axis='y', labelsize=20)
        ax1.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        ax1.set_xlim(0, 360)
        # plt.legend( fontsize=20)
        plt.tight_layout()

        ax2 = plt.subplot2grid((1, 3), (0, 1))
        plt.plot(lowest_dip,lowest_misfits,'bo')

        ymin, ymax = ax2.get_ylim()
        if real_v is not None:
            plt.vlines(real_v[3], ymin=ymin, ymax=ymax, colors='g', linewidth=3)
            plt.vlines(dip, ymin=ymin, ymax=ymax, colors='k', linewidth=3, label='True model')
        # plt.vlines(df_select['Dip'][1], ymin=ymin, ymax=ymax, colors='r', linewidth=3,label='Start model')
        ax2.set_title("Dip", color='b', fontsize=25)
        ax2.set_xlabel("N=%i" % (len(df_select['Dip'])), fontsize=25)
        ax2.xaxis.set_ticks(np.arange(0, 90, 10))
        ax2.tick_params(axis='x', labelsize=20)
        ax2.tick_params(axis='y', labelsize=20)
        ax2.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        ax2.set_xlim(0, 90)
        # plt.legend( fontsize=20)
        plt.tight_layout()

        ax3 = plt.subplot2grid((1, 3), (0, 2))
        plt.plot(lowest_rake,lowest_misfits,'bo')
        ymin, ymax = ax3.get_ylim()
        if real_v is not None:
            plt.vlines(real_v[4], ymin=ymin, ymax=ymax, colors='g', linewidth=3, label='Auxiliary plane')
            plt.vlines(rake, ymin=ymin, ymax=ymax, colors='k', linewidth=3, label='Fault plane')
            plt.legend(loc='upper left', fontsize=20)
            # plt.legend(loc='center left', bbox_to_anchor=(1.1, 0.5))
        ax3.set_title("Rake", color='b', fontsize=25)
        ax3.set_xlabel("N=%i" % (len(df_select['Rake'])), fontsize=25)
        ax3.tick_params(axis='x', labelsize=20)
        ax3.tick_params(axis='y', labelsize=20)
        ax3.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        ax3.set_xlim(-180, 180)
        plt.tight_layout()
        # plt.show()
        plt.savefig(dir + '/Misfit_vs_Moment.pdf')


        fig = plt.figure(figsize=(10, 6))
        row = 0

        ax1 = plt.subplot2grid((3,1), (0, 0))
        plt.plot(df_select['Strike'])
        ax1.set_title("Strike", color='b', fontsize=25)
        ax1.set_xlabel("N=%i" % (len(df_select['Strike'])), fontsize=25)
        ax1.set_ylabel("Misfit", fontsize=25)
        ax1.tick_params(axis='x', labelsize=20)
        ax1.tick_params(axis='y', labelsize=20)
        ax1.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        ax1.set_ylim(0, 360)
        # plt.legend( fontsize=20)
        plt.tight_layout()

        ax2 = plt.subplot2grid((3,1), (1,0))
        plt.plot(df_select['Dip'])

        ax2.set_title("Dip", color='b', fontsize=25)
        ax2.set_xlabel("N=%i" % (len(df_select['Dip'])), fontsize=25)
        ax2.xaxis.set_ticks(np.arange(0, 90, 10))
        ax2.tick_params(axis='x', labelsize=20)
        ax2.tick_params(axis='y', labelsize=20)
        ax2.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        ax2.set_ylim(0, 90)
        # plt.legend( fontsize=20)
        plt.tight_layout()

        ax3 = plt.subplot2grid((3, 1), (2,0))
        plt.plot(df_select['Rake'])
        ax3.set_title("Rake", color='b', fontsize=25)
        ax3.set_xlabel("N=%i" % (len(df_select['Rake'])), fontsize=25)
        ax3.tick_params(axis='x', labelsize=20)
        ax3.tick_params(axis='y', labelsize=20)
        ax3.ticklabel_format(style="sci", axis='y', scilimits=(-2, 2))
        ax3.set_ylim(-180, 180)
        plt.tight_layout()
        # plt.show()
        plt.savefig(dir + '/Trace_plot.pdf')





    def trace_density(self, filepath, savename, directory, skiprows, column_names, real_v, burnin):
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
        strike, dip, rake = aux_plane(real_v[2], real_v[3], real_v[4])
        fig = plt.figure(figsize=(20, 20))
        row = 0

        for i in df_select:
            ax1 = plt.subplot2grid((5, 2), (row, 0))
            ax1.plot(df_select[i], label=i)
            # ax1.axhline(y=REAL[df_select], linewidth=0.3, linestyle=':')
            ax1.set_title("Trace %s" % i, color='b', fontsize=20)
            ax1.set_xlabel("Iteration")
            # ax1.set_ylabel("Epicentral ")
            plt.tight_layout()
            ax2 = plt.subplot2grid((5, 2), (row, 1))
            plt.hist(df_select[i][burnin:], bins=100)
            ymin, ymax = ax2.get_ylim()
            plt.vlines(df_select[i][1], ymin=ymin, ymax=ymax, colors='r', linewidth=3)
            plt.vlines(real_v[row], ymin=ymin, ymax=ymax, colors='k', linewidth=3)
            if i == 'Strike':
                plt.vlines(strike, ymin=ymin, ymax=ymax, colors='g', linewidth=3)
            if i == 'Dip':
                plt.vlines(dip, ymin=ymin, ymax=ymax, colors='g', linewidth=3)
            if i == 'Rake':
                plt.vlines(rake, ymin=ymin, ymax=ymax, colors='g', linewidth=3)

            # y, binEdges = np.histogram(df_select[i], bins=100)
            # bincenters = 0.5 * (binEdges[1:] + binEdges[:-1])
            # pylab.plot(bincenters, y, '-',label = "%s" % i)
            ax2.set_title("Density %s" % i, color='b', fontsize=20)
            ax2.set_xlabel("N=%i" % (len(df_select[i])))
            plt.tight_layout()
            row += 1
        plt.savefig(dir + '/Trace_density.pdf')
        plt.close()

    def get_pdf(self, filepath, savename, directory, skiprows, column_names, show=True):
        dir = directory + '/%s' % (savename.strip('.yaml'))
        if not os.path.exists(dir):
            os.makedirs(dir)

        dir_pdf = dir + '/PDF'
        if not os.path.exists(dir_pdf):
            os.makedirs(dir_pdf)

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
        R_length = int(data[0][-2])
        L_length = int(data[0][-1])

        for i in range(R_length):
            column_names = np.append(column_names, "R_%i" % (i + 1))
        for i in range(L_length):
            column_names = np.append(column_names, "L_%i" % (i + 1))
        column_names = np.append(column_names, "Accepted")
        column_names = np.append(column_names, "Rayleigh_length")
        column_names = np.append(column_names, "Love_length")

        df = pd.DataFrame(data,
                          columns=column_names)
        df_select = df[["Epi", "Depth", "Strike", "Dip", "Rake"]]
        pdf = Plots()
        for i, value in df_select.iteritems():
            pdf.marginal_1D(data=value, name=i, amount_bins=20, directory=dir_pdf, show=show)
        for i in itertools.combinations(df_select, 2):
            pdf.marginal_2D(df_select[i[0]], i[0], df_select[i[1]], i[1], amount_bins=20, directory=dir_pdf, show=show)

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

    def get_BBB(self, filepath, savename, directory, skiprows, column_names, burnin, real_v = None):
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
            data = np.loadtxt(filepath, delimiter=',', skiprows=0)


        df = pd.DataFrame(data, columns=column_names)

        fig_bb, ax_bb = plt.subplots(1, 1, figsize=(4, 4))

        ax_bb.set_xticks([])
        ax_bb.set_yticks([])
        ax_bb.axis('off')
        img = None
        buf = io.BytesIO()

        strike = df['Strike']
        dip = df['Dip']
        rake = df['Rake']

        for i in range(0, len(strike)):

            b = beach(fm=[strike[i], dip[i], rake[i]],
                      width=200, linewidth=0, facecolor='b',
                      xy=(0, 0), axes=ax_bb, alpha=1, zorder=i)
            ax_bb.add_collection(b)
            ax_bb.set_xlim((-0.1, 0.1))
            ax_bb.set_ylim((-0.1, 0.1))

            p = Circle((0., 0,), 0.065, linewidth=2, edgecolor='k', facecolor='b', zorder=5)
            ax_bb.add_patch(p)

            buf.seek(0)
            fig_bb.savefig(buf, format='png')
            buf.seek(0)
            if img is None:
                img = mpimg.imread(buf)
            else:
                img += mpimg.imread(buf)
        plt.close(fig_bb)
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        ax.imshow(img / np.max(img.flatten()))
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('off')
        plt.savefig(dir + '/BBB.pdf')
        # plt.show()

        # beachball([strike,dip,rake], size=200, linewidth=2, facecolor='b', outfile=dir + '/beach_true_aux.pdf')

    def get_stereonets(self, filepath, savename, directory, skiprows, column_names, real_v, burnin):
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
        #
        for i, v in enumerate(df['Strike'][burnin:]):
            if v < 200:
                df['Strike'][i + burnin] = v + 360

        for i, v in enumerate(df['Rake'][burnin:]):
            if v < -100:
                df['Rake'][i + burnin] = -v

        mean_strike = np.mean(df['Strike'][burnin:])
        mean_dip = np.mean(df['Dip'][burnin:])
        mean_rake = np.mean(df['Rake'][burnin:])
        strike, dip, rake = aux_plane(real_v[2], real_v[3], real_v[4])

        print("True solution strike: %.4f" % strike)
        print("True solution dip: %.4f" % dip)
        print("True solution rake: %.4f" % rake)
        # print("SD M0: %.4E" %np.std(df['M0'][burnin:]))
        # print("mean M0: %.4E" %np.mean(df['M0'][burnin:]))
        df_select = df[["Epi", "Depth", "Strike", "Dip", "Rake"]]
        for i in df_select:
            print("SD %s: %.4f" % (i, np.std(df_select[i][burnin:])))
            print("mean %s: %.4f" % (i, np.mean(df_select[i][burnin:])))

        # moment = self.convert(df['Strike'][burnin:], df['Dip'][burnin:], df['Rake'][burnin:], df['M0'][burnin:])
        # M = np.array([np.mean(moment[0]),np.mean(moment[1]),np.mean(moment[2]),np.mean(moment[3]),np.mean(moment[4]),np.mean(moment[5])])
        # M_true = self.convert(real_v[2],real_v[3],real_v[4],3.16227766016798e+16)
        # M_true = self.convert(real_v[2],real_v[3],real_v[4],real_v[5])
        # M_true = self.convert(strike,dip,rake,real_v[5])

        # beachball([np.mean(df['Strike'][burnin:]), np.mean(df['Dip'][burnin:]), np.mean(df['Rake'][burnin:])], size=200, linewidth=2, facecolor='b', outfile=dir + '/beach.pdf')
        beachball([real_v[2], real_v[3], real_v[4]], size=200, linewidth=2, facecolor='b',
                  outfile=dir + '/beach_true.pdf')
        beachball([strike, dip, rake], size=200, linewidth=2, facecolor='b', outfile=dir + '/beach_true_aux.pdf')

        # SD_strike= np.std(df['Strike'][burnin:])
        # SD_dip=np.std(df['Dip'][burnin:])
        # SD_rake=np.std(df['Rake'][burnin:])
        #
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='stereonet')
        # ax.rake(mean_strike-SD_strike, mean_dip-SD_dip, mean_rake-SD_rake, c='b', zorder=5, markersize=4)
        # ax.rake(mean_strike+SD_strike, mean_dip+SD_dip, mean_rake+SD_rake, c='r', zorder=5, markersize=4)
        # ax.plane(mean_strike-SD_strike, mean_dip-SD_dip, c='b', linewidth=0.3)
        # ax.plane(mean_strike+SD_strike, mean_dip+SD_dip, c='r', linewidth=0.3)
        # ax.rake(mean_strike, mean_dip, mean_rake, c='k', zorder=5)
        # ax.plane(mean_strike, mean_dip, c='k', linewidth=0.3,linestyle=':')
        # ax.rake(strike, dip, rake, c='k', zorder=5)
        # ax.plane(strike, dip, c='k', linewidth=0.3)
        # ax.rake(real_v[2], real_v[3], real_v[4], c='k', zorder=5)
        # ax.plane(real_v[2], real_v[3], c='k', linewidth=0.3)
        # # ax.grid()
        # plt.savefig(dir+'/stereonet.pdf')

    def Seaborn_plots(self, filepath, savename, directory, show, skiprows, column_names):
        dir = directory + '/%s' % (savename.strip('.yaml'))
        if not os.path.exists(dir):
            os.makedirs(dir)

        dir_seaborn = dir + '/Seaborn'
        if not os.path.exists(dir_seaborn):
            os.makedirs(dir_seaborn)

        if filepath.endswith('.yaml') == True:
            with open(filepath, 'r') as stream:
                data_file = yaml.load(stream)
                stream.close()
                data = data_file['data']
                # parameters = data_file['parameters']
        else:
            # parameters = open(filepath, "r").readlines()[:33]
            # file = '/home/nienke/Documents/Applied_geophysics/Thesis/anaconda/Final/Random_Result/Exploring.txt'
            data = np.loadtxt(filepath, delimiter=',', skiprows=skiprows)
        length = len(data[0]) - (len(column_names) + 3)
        R_length = int(data[0][-2])
        L_length = int(data[0][-1])

        for i in range(R_length):
            column_names = np.append(column_names, "R_%i" % (i + 1))
        for i in range(L_length):
            column_names = np.append(column_names, "L_%i" % (i + 1))
        column_names = np.append(column_names, "Accepted")
        column_names = np.append(column_names, "Rayleigh_length")
        column_names = np.append(column_names, "Love_length")

        df = pd.DataFrame(data,
                          columns=column_names)
        df_select = df[["Epi", "Depth", "Strike", "Dip", "Rake"]]
        params = {'legend.fontsize': 'x-large',
                  'figure.figsize': (15, 15),
                  'axes.labelsize': 25,
                  'axes.titlesize': 'x-large',
                  'xtick.labelsize': 25,
                  'ytick.labelsize': 25}
        pylab.rcParams.update(params)
        plt.figure()
        # autocorrelation_plot(df_select['Strike'],linewidth = 1,label = "Strike")
        # autocorrelation_plot(df_select['Dip'],linewidth = 1, label = "Dip")
        # autocorrelation_plot(df_select['Rake'],linewidth = 1, label = "Rake")
        # autocorrelation_plot(df_select['Epi'],linewidth = 1, label = "Epicentral_distance")
        # autocorrelation_plot(df_select['Depth'],linewidth = 1, label = "Depth")
        # plt.legend()
        # plt.savefig(directory + '/autocorrelation.pdf')
        # plt.show()
        # params = {'legend.fontsize': 'x-large',
        #           'figure.figsize': (15, 15),
        #           'axes.labelsize': 25,
        #           'axes.titlesize': 'x-large',
        #           'xtick.labelsize': 25,
        #           'ytick.labelsize': 25}
        # pylab.rcParams.update(params)
        # # plt.figure()
        scatter_matrix(df_select, diagonal='kde')
        plt.savefig(directory + '/correlations.pdf')
        # #

        #
        #

        # plot = Plots()
        # for i in itertools.combinations(df_select, 2):
        #     plot.Kernel_density(data=df_select, name_x=i[0], name_y=i[1], directory=dir_seaborn,
        #                         savename=savename.strip(".yaml"), show=show)
        #     plot.hist(data=df, name_x=i[0], name_y=i[1], directory=dir_seaborn, savename=savename.strip(".yaml"),
        #               show=show)

        ## Pair Grid approximation
        # plot.Pair_Grid(data=df_select,directory=dir_seaborn,savename=savename,show=show)

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


    def get_seismogram_plots(self, directory, sdr=False):
        if sdr == True:
            dir = directory + '/proc_sdr'
        else:
            dir = directory + '/proc'
        filenames = glob.glob("%s/*.yaml" % dir)
        for file in filenames:
            with open(file, 'r') as stream:
                data = yaml.load(stream)
                stream.close
            f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
            for i, v in data.iteritems():
                ax1.plot(v['trace_z'], alpha=0.2)
                ax2.plot(v['trace_r'], alpha=0.2)
                ax3.plot(v['trace_t'], alpha=0.2)
            plt.show()

    def plot_streams(self, stream_filepath, filepath, savename, directory, skiprows, column_names):
        dir = directory + '/%s' % (savename.strip('.yaml'))
        if not os.path.exists(dir):
            os.makedirs(dir)

        dir_pdf = dir + '/stream_plots'
        if not os.path.exists(dir_pdf):
            os.makedirs(dir_pdf)

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
        R_length = int(data[0][-2])
        L_length = int(data[0][-1])

        for i in range(R_length):
            column_names = np.append(column_names, "R_%i" % (i + 1))
        for i in range(L_length):
            column_names = np.append(column_names, "L_%i" % (i + 1))
        column_names = np.append(column_names, "Accepted")
        column_names = np.append(column_names, "Rayleigh_length")
        column_names = np.append(column_names, "Love_length")

        st = obspy.read(stream_filepath)
        for i in st.traces:
            if "Z" in i.meta.channel:
                a = 1

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
        plt.figure(figsize=(8, 8))
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
                axis.set_xlabel(f"dimension {dimensions_list[number_of_plots - 1]}")

            axis.set_yticklabels([])
            axis.tick_params(axis="y", which="both", left=False, right=False)
            if i_plot == 0:
                axis.set_ylabel("relative density")

            # Plot histogram on diagonal
            axis.hist(
                samples[:,dimensions_list[i_plot]],
                bins=bins,
                density=False,
                range=dim_range[i_plot],
                color=color_1d,
            )

            for j_plot in range(i_plot):
                # print(i_plot, j_plot) # grid indices for lower left
                axis = plt.subplot(gs1[j_plot + (number_of_plots) * i_plot])

                # Modify axes
                if i_plot != number_of_plots - 1:
                    axis.set_xticklabels([])
                    axis.tick_params(axis="x", which="both", bottom=False, top=False)
                else:
                    axis.set_xlabel(f"dimension {dimensions_list[j_plot]}")

                if j_plot != 0:
                    axis.set_yticklabels([])
                    axis.tick_params(axis="y", which="both", left=False, right=False)
                else:
                    axis.set_ylabel(f"dimension {dimensions_list[i_plot]}")

                # Plot 2d marginals
                axis.hist2d(
                    samples[:,dimensions_list[j_plot]],
                    samples[:,dimensions_list[i_plot]],
                    bins=bins,
                    range=[dim_range[j_plot], dim_range[i_plot]],
                    cmap=colormap_2d,
                )

                # print(i_plot, j_plot) # grid indices for lower left
                axis = plt.subplot(gs1[i_plot + (number_of_plots) * j_plot])

                axis.set_xticklabels([])
                axis.tick_params(axis="x", which="both", bottom=False, top=False)
                axis.set_yticklabels([])
                axis.tick_params(axis="y", which="both", left=False, right=False)

                axis.axis("off")

                correlation = np.corrcoef(
                    samples[:,dimensions_list[j_plot]], samples[:,dimensions_list[i_plot]]
                )[1][0]
                axis.text(
                    0.5,
                    0.5,
                    f"{correlation:.2f}",
                    horizontalalignment="center",
                    verticalalignment="center",
                    fontsize=40 * np.abs(correlation),
                    transform=axis.transAxes,
                )

        if show:
            plt.show()
        else:
            plt.savefig(dir + '/Marginals.pdf')



if __name__ == '__main__':
    main()
