#!/usr/bin/python
# -*- coding: utf-8 -*-
# 2021-8-13 22:44:15 update this code
from pylab import *
import numpy as np
import os
import glob

#
# # spinless080804_N12_0816_contour_load_merge1115_v4_newdata1116.py
# def colorbar_mesh_small(mappable, title):
#     from mpl_toolkits.axes_grid1 import make_axes_locatable
#     import matplotlib.pyplot as plt
#     from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#     last_axes = plt.gca()
#     ax = mappable.axes
#     fig = ax.figure
#     divider = make_axes_locatable(ax)
#     cax = divider.append_axes("right", size="5%", pad=0.05)
#
#     label_size = 22
#     cbar = fig.colorbar(mappable, cax=cax, orientation='vertical')
#     cbar.ax.tick_params(labelsize=label_size - 4)
#     cbar.ax.set_title(title, fontsize=label_size - 4)
#     plt.sca(last_axes)
#
#     return cbar


def reading_data(filename, str_to_find):
    try:
        f_in = open(filename, 'r')
        word_shift = len(str_to_find.split()) - 2
        for line in f_in:
            if str_to_find in line:
                data_ave = line.split()[2 + word_shift]
                data_err = line.split()[4 + word_shift]
                break

        f_in.close()
        return data_ave, data_err
    except:
        return np.nan, np.nan

def reading_mu(filename, str_to_find):
    try:
        f_in = open(filename, 'r')
        word_shift = 0#len(str_to_find.split()) - 2
        for line in f_in:
            if str_to_find in line:
                mu = line.split()[2]
                mu = float(mu)
                break

        f_in.close()
        return mu
    except:
        return np.nan


##################################################################
def get_ave(data, num_realizations):
    ave_data = np.nanmean(data)
    err_data = np.std(data) / np.sqrt(num_realizations)
    # err_data = np.nan
    return ave_data, err_data


##################################################################

######################## Plotting  parameters #############################################
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True

list_of_ls = [80, 160, 200, 240]
# list_of_ls = [200, 200]
# list_of_mus = [0, -0.5]
list_of_mus = [0.5, 0, -0.5, -1]
# list_of_mus = [0,0]
nrows = len(list_of_mus) + 1
ncols = 3*2
fig, axs = plt.subplots(nrows, ncols, sharex=True, figsize=(8 * ncols, 3 * nrows))

txt_size = 20
label_size = 20
########################################################################################################################
Lx = 8
Ly = 4
Lx = 10
Ly = 4
# l = 240
U = - 7

# load = True
# load = False
file_name = 'ladder'
path = './' + file_name + '/'

# list_of_Us = data[0, :]
list_of_colors = ['crimson', 'blue', 'forestgreen', 'black', 'purple', 'orange', 'darkblue', 'rosybrown', 'green', 'dimgrey', 'limegreen']

# get_dat = False
get_dat = True
if get_dat:

    for i_mu, mu in enumerate(list_of_mus):

        for i_l, l in enumerate(list_of_ls):
            #
            # list_of_files = glob.glob('./' + file_name + '_L_%d_l_%d_U_%.4f_mu_%.4f.dat' % (Lx, l, U, mu))
            # num_files = len(list_of_files)
            # if len(list_of_files) > 0:
            #     print('dat file exsit!')
            #     continue
            # else:

            # if load:
            #     # txt_file_name = file_name + '_L_%d_l_%d_U_%.4f_mu_%.4f.dat' % (Lx, l, U, mu)
            #     # # np.savetxt(txt_file_name, np.c_[list_of_hs, ave_sign, err_sign, ave_sign_sigma, err_sign_sigma, ave_dens, err_dens])
            #     #
            #     # data = np.loadtxt(txt_file_name, unpack=True)
            #     # list_of_hs = data[:, 0]
            #     # ave_sign = data[:, 1]
            #     # err_sign = data[:, 2]
            #     # ave_sign_sigma = data[:, 3]
            #     # err_sign_sigma = data[:, 4]
            #     # ave_dens = data[:, 5]
            #     # err_dens = data[:, 6]
            #     continue

            print(' ')
            #rLx8Ly2l200U-7.00mu0.00h0.00r000
            # bn      rLx8Ly2l200U-7.00mu0.00h0.00r000.out
            #./ladder/rLx8rLy2l200U-7.00mu0.00h*r000.out
            list_of_files = glob.glob(path + 'rLx%dLy%dl%dU%.2fmu%.2fh*r000.out' % (Lx, Ly, l, U, mu))  # rL16l200U-7.00mu0.50h0.00r000.out
            # print(path + 'rLx%dLy%dl%dU%.2fmu%.2fh*r000.out' % (Lx, Ly, l, U, mu))
            list_of_hs = np.full((len(list_of_files)), np.nan)
            if len(list_of_files) == 0:
                continue
            # else:
                # print('exist')

            for i_file, filename in enumerate(list_of_files):
                mu_up = reading_mu(filename, 'mu_up :')
                mu_dn = reading_mu(filename, 'mu_dn :')
                list_of_hs[i_file] = (mu_dn - mu_up) / 2
                # print(i_file, filename, list_of_hs[i_file])

            print(list_of_hs)
            # list_of_hs = np.sort(list_of_hs)
            # print(list_of_hs)

            ave_dens = np.full(len(list_of_hs), np.nan)
            err_dens = np.full(len(list_of_hs), np.nan)
            ave_sign = np.full(len(list_of_hs), np.nan)
            err_sign = np.full(len(list_of_hs), np.nan)

            ave_sign_up = np.full(len(list_of_hs), np.nan)
            err_sign_up = np.full(len(list_of_hs), np.nan)
            ave_sign_dn = np.full(len(list_of_hs), np.nan)
            err_sign_dn = np.full(len(list_of_hs), np.nan)
            ave_sign_sigma = np.full(len(list_of_hs), np.nan)
            err_sign_sigma = np.full(len(list_of_hs), np.nan)

            #####

            ave_dens = np.full(len(list_of_hs), np.nan)
            err_dens = np.full(len(list_of_hs), np.nan)

            ave_dens_up = np.full(len(list_of_hs), np.nan)
            err_dens_up = np.full(len(list_of_hs), np.nan)
            ave_dens_dn = np.full(len(list_of_hs), np.nan)
            err_dens_dn = np.full(len(list_of_hs), np.nan)
            ave_dens_sigma = np.full(len(list_of_hs), np.nan)
            err_dens_sigma = np.full(len(list_of_hs), np.nan)
            for i_h, h in enumerate(list_of_hs):
                list_of_files = glob.glob(path + 'rLx%dLy%dl%dU%.2fmu%.2fh%.2fr*.out' % (Lx, Ly, l, U, mu, h))  # rL16l200U-7.00mu0.50h0.00r000.out
                num_files = len(list_of_files)
                if len(list_of_files) == 0:
                    continue
                # else:
                #     print('load')

                r_ave_dens = np.full(len(list_of_files), np.nan)
                r_err_dens = np.full(len(list_of_files), np.nan)
                r_ave_sign = np.full(len(list_of_files), np.nan)
                r_err_sign = np.full(len(list_of_files), np.nan)

                r_ave_sign_up = np.full(len(list_of_files), np.nan)
                r_err_sign_up = np.full(len(list_of_files), np.nan)
                r_ave_sign_dn = np.full(len(list_of_files), np.nan)
                r_err_sign_dn = np.full(len(list_of_files), np.nan)
                r_ave_sign_sigma = np.full(len(list_of_files), np.nan)
                r_err_sign_sigma = np.full(len(list_of_files), np.nan)


                r_ave_dens = np.full(len(list_of_files), np.nan)
                r_err_dens = np.full(len(list_of_files), np.nan)

                r_ave_dens_up = np.full(len(list_of_files), np.nan)
                r_err_dens_up = np.full(len(list_of_files), np.nan)
                r_ave_dens_dn = np.full(len(list_of_files), np.nan)
                r_err_dens_dn = np.full(len(list_of_files), np.nan)
                r_ave_dens_sigma = np.full(len(list_of_files), np.nan)
                r_err_dens_sigma = np.full(len(list_of_files), np.nan)
                for i_file, filename in enumerate(list_of_files):
                    r_ave_dens[i_file], r_err_dens[i_file] = reading_data(filename, 'Density :')
                    r_ave_sign[i_file], r_err_sign[i_file] = reading_data(filename, 'Avg sign :')
                    r_ave_sign_up[i_file], r_err_sign_up[i_file] = reading_data(filename, 'Avg up sign :')
                    r_ave_sign_dn[i_file], r_err_sign_dn[i_file] = reading_data(filename, 'Avg dn sign :')
                    r_ave_sign_sigma[i_file] = (r_ave_sign_up[i_file] + r_ave_sign_dn[i_file]) / 2


                    r_ave_dens[i_file], r_err_dens[i_file] = reading_data(filename, 'Density :')
                    r_ave_dens_up[i_file], r_err_dens_up[i_file] = reading_data(filename, 'Up spin occupancy :')
                    r_ave_dens_dn[i_file], r_err_dens_dn[i_file] = reading_data(filename, 'Down spin occupancy :')
                    r_ave_dens_sigma[i_file] = (r_ave_dens_up[i_file] + r_ave_dens_dn[i_file]) / 2

                ave_dens[i_h], err_dens[i_h] = get_ave(r_ave_dens, num_files)
                ave_sign[i_h], err_sign[i_h] = get_ave(r_ave_sign, num_files)
                ave_sign_up[i_h], err_sign_up[i_h] = get_ave(r_ave_sign_up, num_files)
                ave_sign_dn[i_h], err_sign_up[i_h] = get_ave(r_ave_sign_dn, num_files)
                ave_sign_sigma[i_h] = (ave_sign_up[i_h] + ave_sign_dn[i_h]) / 2
                err_sign_sigma[i_h] = (err_sign_up[i_h] + err_sign_dn[i_h]) / 2
                ave_dens[i_h], err_dens[i_h] = get_ave(r_ave_dens, num_files)
                ave_dens_up[i_h], err_dens_up[i_h] = get_ave(r_ave_dens_up, num_files)
                ave_dens_dn[i_h], err_dens_up[i_h] = get_ave(r_ave_dens_dn, num_files)
                ave_dens_sigma[i_h] = (ave_dens_up[i_h] + ave_dens_dn[i_h]) / 2
                # err_dens_sigma[i_h] = (err_dens_up[i_h] + err_dens_dn[i_h]) / 2
                # ave_HD[i_h], err_HD[i_h] = get_ave(r_ave_HD, num_files)
                print('l = %d   mu = %.4f   h = %.4f   ave_sign = %.4f   ave_sign_sigma = %.4f   ' % (l, mu, h, ave_sign[i_h], ave_sign_sigma[i_h]))

            # print((Lx, l, U, mu))
            print(ave_dens_up[0:3])

            txt_file_name = './dat/chain_Lx_%d_Ly_%d_l_%d_U_%.4f_mu_%.4f.dat' % (Lx, Ly, l, U, mu)
            # txt_file_name = file_name + '_L_%d_l_%d_U_%.4f_mu_%.4f.dat' % (Lx, l, U, mu)
            # np.savetxt(txt_file_name, np.c_[list_of_hs, ave_sign, err_sign, ave_sign_sigma, err_sign_sigma, ave_dens, err_dens])
            np.savetxt(txt_file_name, np.c_[list_of_hs, ave_sign, err_sign, ave_sign_up, err_sign_up, ave_sign_dn, err_sign_dn, ave_dens, err_dens, ave_dens_up, err_dens_up, ave_dens_dn, err_dens_dn])
            continue
            # ################################################################################################
            #
            #
            # # ave_signs[:, i_l] = ave_sign
            # # ave_sign_sigmas[:, i_l] = ave_sign_sigma
            # # ave_sign_ups[:, i_l] = ave_sign_up
            # # ave_sign_dns[:, i_l] = ave_sign_dn
            ########################################################################################################################
            # for i_l, mu in enumerate(list_of_mus):
            for i in range(3):
                if i == 0: z = ave_sign
                # if i == 1: z = ave_sign_sigmas
                if i == 1: z = ave_sign_up
                if i == 2: z = ave_sign_dn

                # ax = axs[i]
                ax = axs[i_mu][i]  # [i_l]

                # p1, = ax.plot(list_of_hs, z, color='crimson', marker='o', linestyle='--', label=r'$\mu/t=%.2f$' % mu)
                # p1, = ax.plot(list_of_hs, z, color='crimson', marker='o', linestyle='--', label=r'$\mu/t=%.2f$' % mu), label=r'$L_\tau=%d$'%l
                ax.plot(list_of_hs, z, color=list_of_colors[i_l], marker='o', linestyle='--', label=r'$L_\tau=%d$' % l)
    # 1
# 1
# if True:
for i_mu, mu in enumerate(list_of_mus):
    # for i_mu in range(len(list_of_mus)+1):
    #     if i_mu < len(list_of_mus):
    #         mu = list_of_mus[i_mu]
    #
    ax = axs[i_mu][0]

    # axs[0][i_mu].set_title(r'$L_{\tau}=%d $' % (list_of_ls[i_l]), fontsize=txt_size)
    # axs[0][i_mu].set_title(r'$\mu/t=%d $' % (mu - 3.5), fontsize=txt_size)

    for i_l, l in enumerate(list_of_ls):

        if True:
            list_of_files = glob.glob('./dat/chain_Lx_%d_Ly_%d_l_%d_U_%.4f_mu_%.4f.dat' % (Lx, Ly, l, U, mu))  # rL16l200U-7.00mu0.50h0.00r000.out

            if len(list_of_files) == 0:
                continue
            txt_file_name = './dat/chain_Lx_%d_Ly_%d_l_%d_U_%.4f_mu_%.4f.dat' % (Lx, Ly, l, U, mu)
            # np.savetxt(txt_file_name, np.c_[list_of_hs, ave_sign, err_sign, ave_sign_sigma, err_sign_sigma, ave_dens, err_dens])
            #np.savetxt(txt_file_name, np.c_[list_of_hs, ave_sign, err_sign, ave_sign_up, err_sign_up, ave_sign_dn, err_sign_up, ave_dens, err_dens])

            # data = np.loadtxt(txt_file_name, unpack=True)

            data = np.swapaxes(np.loadtxt(txt_file_name, unpack=True), 0, 1)
            list_of_hs = data[:, 0]
            ave_sign = data[:, 1]
            # print(ave_sign[0:10])
            err_sign = data[:, 2]
            ave_sign_up = data[:, 3]
            err_sign_up = data[:, 4]
            ave_sign_dn = data[:, 5]
            err_sign_dn = data[:, 6]
            ave_dens = data[:, 7]
            err_dens = data[:, 8]
            ave_dens_up = data[:, 9]
            err_dens_up = data[:, 10]
            ave_dens_dn = data[:, 11]
            err_dens_dn = data[:, 12]
            # list_of_hs = data[:, 0]
            # ave_sign = data[:, 1]
            # err_sign = data[:, 2]
            # # ave_sign_sigma = data[:, 3]
            # # err_sign_sigma = data[:, 4]
            # ave_dens = data[:, 5]
            # err_dens = data[:, 6]
            # continue



        for i in range(6):

            if i == 0: z = ave_sign
            if i == 1: z = ave_sign_up
            if i == 2: z = ave_sign_dn

            if i == 0: zerr = err_sign
            if i == 1: zerr = err_sign_up
            if i == 2: zerr = err_sign_dn

            if i == 3: z = ave_dens
            if i == 4: z = ave_dens_up
            if i == 5: z = ave_dens_dn

            if i == 3: zerr = err_dens
            if i == 4: zerr = err_dens_up
            if i == 5: zerr = err_dens_dn

            # if i == 1: z = ave_sign_sigmas
            # if i == 0: z = ave_sign
            # # if i == 1: z = ave_sign_sigmas
            # if i == 1: z = ave_sign_up
            # if i == 2: z = ave_sign_dn
            #
            # # ax = axs[i]
            ax = axs[i_mu][i]  # [i_l]

            # p1, = ax.plot(list_of_hs, z, color='crimson', marker='o', linestyle='--', label=r'$\mu/t=%.2f$' % mu)
            # p1, = ax.plot(list_of_hs, z, color='crimson', marker='o', linestyle='--', label=r'$\mu/t=%.2f$' % mu), label=r'$L_\tau=%d$'%l
            # ax.plot(list_of_hs, z, color=list_of_colors[i_l], marker='o', linestyle='--', label=r'$L_\tau=%d$' % l)

            # ax.errorbar(np.arange(0, 30, 1), abs(ave_Kc), yerr=err_Kc, linestyle=' 'marker='o',
            #             color='blue')
            ax.errorbar(list_of_hs, z, yerr=zerr, color=list_of_colors[i_l], marker='o', linestyle='--', capsize=4, label=r'$L_\tau=%d$' % l)
            #
            # list_of_colors = ['crimson', 'blue', 'forestgreen', 'purple', 'orange', 'darkblue', 'rosybrown', 'green',
            #                   'dimgrey', 'limegreen', 'black']
            #
            # # p1, = ax.plot(list_of_hs, z, color='crimson', marker='o', linestyle='--', label=r'$\mu/t=%.2f$' % mu)
            # # p1, = ax.plot(list_of_hs, z, color='crimson', marker='o', linestyle='--', label=r'$\mu/t=%.2f$' % mu), label=r'$L_\tau=%d$'%l
            # p1, = ax.plot(list_of_hs, z, color=list_of_colors[i_l], marker='o', linestyle='--', label=r'$L_\tau=%d$' % l)
            # p1, = ax.plot(list_of_hs, z[:, i_l], color='crimson', marker='o', linestyle='--', label=r'$\mu/t=%.2f$' % mu)
            # plot_pcolor = ax.pcolormesh(list_of_hs - dh / 2, list_of_mus - dmu / 2 - 3.5, np.swapaxes(z, 0, 1), cmap='RdBu',
            #                             edgecolors='face', vmin=0, vmax=1.0)  # , ticks=[0.46, 0.5])
            # plot_pcolor = ax.pcolormesh(list_of_hs-dh/2, list_of_mus-dmu/2 , np.swapaxes(z, 0, 1), cmap='RdBu',  edgecolors='face', vmin=0, vmax=1.0)#, ticks=[0.46, 0.5])
            # if i == 0: colorbar_mesh_small(plot_pcolor, r'${\langle S \rangle}$')
            # if i == 1: colorbar_mesh_small(plot_pcolor, r'${\langle S \rangle}_{\sigma}$')
            # if i == 1: colorbar_mesh_small(plot_pcolor, r'${\langle S_{\uparrow} \rangle}$')
            # if i == 2: colorbar_mesh_small(plot_pcolor, r'${\langle S_{\downarrow} \rangle}$')

            if i == 0: ax.set_ylabel(r'${\langle S \rangle}$', fontsize=txt_size)
            if i == 1: ax.set_ylabel(r'${\langle S_{\uparrow} \rangle}$', fontsize=txt_size)
            if i == 2: ax.set_ylabel(r'${\langle S_{\downarrow} \rangle}$', fontsize=txt_size)
            ax.tick_params(axis='both', which='major', labelsize=label_size, length=10)
            ax.minorticks_on()
            ax.set_xlabel(r'$h/t$', fontsize=txt_size)
            # ax.set_ylabel(r'$\mu/t$', fontsize=txt_size)
            # ax.set_title(r'$ L = %d$' % N, fontsize=txt_size)
            ax.set_xlim(0, 8)
            ax.set_ylim(-0.05, 1.05)
            ymin, ymax = ax.get_ylim()

            txt_size1 = txt_size + 5

            if abs((mu - 3.5) - (-3.5)) < 1e-6:
                Uc = 7.13
                ax.plot([Uc, Uc], [ymin, ymax], color='black', marker=' ', linestyle='--')  # , label=txt_label)
            if abs((mu - 3.5) - (-4)) < 1e-6:
                Uc = 1.516
                ax.plot([Uc, Uc], [ymin, ymax], color='black', marker=' ', linestyle='--')  # , label=txt_label)
                Uc = 6
                ax.plot([Uc, Uc], [ymin, ymax], color='black', marker=' ', linestyle='--')  # , label=txt_label)
                Uc = 7.62
                ax.plot([Uc, Uc], [ymin, ymax], color='black', marker=' ', linestyle='--')  # , label=txt_label)
            if abs((mu - 3.5) - (-4.5)) < 1e-6:
                Uc = 2.09
                ax.plot([Uc, Uc], [ymin, ymax], color='black', marker=' ', linestyle='--')  # , label=txt_label)
            #
            # Uc = 1.74
            # ax.plot([Uc, Uc], [ymin, ymax], color='black', marker=' ', linestyle='--')  # , label=txt_label)
            # Uc = 6.6
            # ax.plot([Uc, Uc], [ymin, ymax], color='black', marker=' ', linestyle='--')  # , label=txt_label)
            # # if i_l == 0:
            #     # ax.plot([1, 1.4 + 1], [-4.6, -5.0 - 1], color='blue', marker=' ')  # , label=txt_label)
            #
            #     # ax.text(0.08, 0.25, r"$V$", transform=ax.transAxes, fontsize=txt_size1, va='center', ha='center', color='w')
            dy = 0
            # if i == 1: dy = 0.5
            # ax.text(0.5, 0.25, r"$PP$", transform=ax.transAxes, fontsize=txt_size1, va='center', ha='center', color='w')
            ax.text(0.08, 0.2 + dy, r"$ED$", transform=ax.transAxes, fontsize=txt_size1, va='center', ha='center',
                    color='black')
            ax.text(0.5, 0.2 + dy, r"$PP$", transform=ax.transAxes, fontsize=txt_size1, va='center', ha='center', color='black')
            # ax.text(0.5, 0.2 + dy, r"$FP_1$", transform=ax.transAxes, fontsize=txt_size1, va='center', ha='center', color='black')
            ax.text(0.99, 0.2 + dy, r"$FP_2$", transform=ax.transAxes, fontsize=txt_size1, va='center', ha='center',
                    color='black')
            #
            # if i == 0: ax.set_ylabel(r'${\langle S \rangle}$', fontsize=txt_size)
            # if i == 1: ax.set_ylabel(r'${\langle S_{\uparrow} \rangle}$', fontsize=txt_size)
            # if i == 2: ax.set_ylabel(r'${\langle S_{\downarrow} \rangle}$', fontsize=txt_size)

axs[0][0].legend(loc='upper right', numpoints=1, fontsize=txt_size-6)

# for i in range(nrows-1):
for i_mu, mu in enumerate(list_of_mus):
    ax = axs[i_mu][0]
    ax.text(-0.2, 0.5, r'$\mu/t = %.1f$' % (mu - 3.5), color='black', bbox=dict(facecolor='green', edgecolor='green', boxstyle='round', alpha=0.2), transform=ax.transAxes, rotation=90, va='center', fontsize=txt_size - 6)

# assert False
txt_size1 = txt_size
txt_size2 = txt_size1 - 4
x1 = 0.08 + 0.02 + 0.02
x3 = 0.92
# i = 1
for i in range(3):
    ax = axs[len(list_of_mus)][i]
    y = 0.3
    ax.text(x1, 0.5+y, r'${\langle S\rangle}=1$', transform=ax.transAxes, fontsize=txt_size2, va='center',
            ha='center', color='black')
    ax.text(x1, 0.3+y, r'${\langle S_{\uparrow} \rangle}<1$', transform=ax.transAxes, fontsize=txt_size2,
            va='center', ha='center', color='black')
    ax.text(x1, 0.1+y, r'${\langle S_{\downarrow} \rangle}<1$', transform=ax.transAxes, fontsize=txt_size2,
            va='center', ha='center', color='black')
    ax.text(0.5, 0.5+y, r'${\langle S \rangle}<1$', transform=ax.transAxes, fontsize=txt_size2, va='center',
            ha='center', color='black')
    # ax.text(0.5, 0.3+y, r'${\langle S_{\uparrow} \rangle}=1$', transform=ax.transAxes, fontsize=txt_size2,
    #         va='center', ha='center', color='black')
    ax.text(0.5, 0.3+y, r'${\langle S_{\uparrow} \rangle}<1$', transform=ax.transAxes, fontsize=txt_size2,
            va='center', ha='center', color='black')

    ax.text(0.5, 0.1+y, r'${\langle S_{\downarrow} \rangle}<1$', transform=ax.transAxes, fontsize=txt_size2,
            va='center', ha='center', color='black')
    ax.text(x3, 0.5+y, r'${\langle S \rangle}=1$', transform=ax.transAxes, fontsize=txt_size2, va='center',
            ha='center', color='black')
    ax.text(x3, 0.3+y, r'${\langle S_{\uparrow} \rangle}=1$', transform=ax.transAxes, fontsize=txt_size2,
            va='center', ha='center', color='black')
    ax.text(x3, 0.1+y, r'${\langle S_{\downarrow} \rangle}=1$', transform=ax.transAxes, fontsize=txt_size2,
            va='center', ha='center', color='black')
    #########################################################

    if i == 0: ax.set_ylabel(r'${\langle S \rangle}$', fontsize=txt_size)
    if i == 1: ax.set_ylabel(r'${\langle S_{\uparrow} \rangle}$', fontsize=txt_size)
    if i == 2: ax.set_ylabel(r'${\langle S_{\downarrow} \rangle}$', fontsize=txt_size)
    ax.tick_params(axis='both', which='major', labelsize=label_size, length=10)
    ax.minorticks_on()
    ax.set_xlabel(r'$h/t$', fontsize=txt_size)
    # ax.set_ylabel(r'$\mu/t$', fontsize=txt_size)
    # ax.set_title(r'$ L = %d$' % N, fontsize=txt_size)
    ax.set_xlim(0, 8)
    ax.set_ylim(-0.05, 1.05)
    ymin, ymax = ax.get_ylim()

    txt_size1 = txt_size + 5
    #
    # Uc = 1.74
    # ax.plot([Uc, Uc], [ymin, ymax], color='black', marker=' ', linestyle='--')  # , label=txt_label)
    # Uc = 6.6
    # ax.plot([Uc, Uc], [ymin, ymax], color='black', marker=' ', linestyle='--')  # , label=txt_label)

    # if i_l == 0:
    #     # ax.plot([1, 1.4 + 1], [-4.6, -5.0 - 1], color='blue', marker=' ')  # , label=txt_label)
    #
    #     # ax.text(0.08, 0.25, r"$V$", transform=ax.transAxes, fontsize=txt_size1, va='center', ha='center', color='w')
    dy = 0
    # if i == 1: dy = 0.5
    # ax.text(0.5, 0.25, r"$PP$", transform=ax.transAxes, fontsize=txt_size1, va='center', ha='center', color='w')
    ax.text(0.08, 0.2 + dy, r"$ED$", transform=ax.transAxes, fontsize=txt_size1, va='center', ha='center',
            color='black')
    ax.text(0.5, 0.2 + dy, r"$PP$", transform=ax.transAxes, fontsize=txt_size1, va='center', ha='center', color='black')
    # ax.text(0.5, 0.2 + dy, r"$FP_1$", transform=ax.transAxes, fontsize=txt_size1, va='center', ha='center', color='black')
    ax.text(0.99, 0.2 + dy, r"$FP_2$", transform=ax.transAxes, fontsize=txt_size1, va='center', ha='center',
            color='black')
    #
    # if i == 0: ax.set_ylabel(r'${\langle S \rangle}$', fontsize=txt_size)
    # if i == 1: ax.set_ylabel(r'${\langle S_{\uparrow} \rangle}$', fontsize=txt_size)
    # if i == 2: ax.set_ylabel(r'${\langle S_{\downarrow} \rangle}$', fontsize=txt_size)

# if i == 0 and i_l == 0:
fig.subplots_adjust(hspace=0.0)
# axs[0][0].set_title(r'${2D}\ \ L_x=8\ \ L_y=2\ \ L_{\tau}=200 $', fontsize=txt_size)
# axs[1].set_title(r'$nwarm=800\ \ npass=2000\ \ (\rm{shift}\ \rm{ finished})$', fontsize=txt_size)
# axs[0][1].set_title(r'$nwarm=800\ \ npass=2000$', fontsize=txt_size)
# axs[0][1].set_title(r'$\mu/t=%.2f$'%(list_of_mus[1]-3.5), fontsize=txt_size)
# axs[0][2].set_title(r'$  (\rm{shift}\ \rm{ finished})$', fontsize=txt_size)
#
# fig.suptitle(r'${1D}\ \ L=%d\ \ \mu=%.2f \ \ $'  % (Lx, mu-3.5)+ r'${\rm{nwarm}}=800\ \ {\rm{npass}}=2000 \ \ 24 \ {\rm{samples}}$', fontsize = txt_size, y = 0.95)
# fig.suptitle(r'${1D}\ \ L=%d \ \ $'  % (Lx)+ r'${\rm{nwarm}}=800\ \ {\rm{npass}}=2000 \ \ 24 \ {\rm{samples}}$', fontsize = txt_size, y = 0.95)
fig.suptitle(r'${2D}\ \ L_x=%d \ \ L_y=%d \ \ $'  % (Lx, Ly)+ r'${\rm{nwarm}}=800\ \ {\rm{npass}}=2000 \ \ 24 \ {\rm{samples}}$', fontsize = txt_size, y=0.91)
#################### Outputting figure ###################################################################
# import glob
print('save')
output_path = './'
ii = 0
# # str_add = '_slice_Tvar_merge' + '_L_%d' % (Lx) + '_mu_%.2f' % (mu-3.5)
# str_add = '_slice_Tvar_merge' + '_L_%d' % (Lx) + '_loop_mu' #% (mu-3.5)
# while 1:
#     # file_name = file_name
#     # png_file_name = file_name + '_mesh_mushift_half_Lx_%d_Ly_%d_l_%d_U_%.4f' % (Lx, Ly, l, U) + '_%d.png' % ii
#     png_file_name = file_name + str_add + '_%d.png' % ii
#     list_of_files = glob.glob(output_path + png_file_name)
#     if len(list_of_files) > 0:
#         print('ii = %d   len(list_of_files_all) > 0' % ii)
#         ii = ii + 1
#     else:
#         break
#
# png_file_name = os.path.join(output_path, png_file_name)
# png_file_name = './ladder_Lx_%d_Ly_%d_5.png' % (Lx, Ly) # 8 2
# png_file_name = './ladder_Lx_%d_Ly_%d_0930.png' % (Lx, Ly) #
png_file_name = './ladder_Lx_%d_Ly_%d_1004.png' % (Lx, Ly) #
plt.savefig(png_file_name, facecolor='w', edgecolor='w', orientation='portrait', format='png', bbox_inches='tight')
#
# pdf_file_name = file_name + str_add + '_%d.pdf' % ii
# list_of_files = glob.glob(output_path + pdf_file_name)
# pdf_file_name = os.path.join(output_path, pdf_file_name)
# plt.savefig(pdf_file_name, facecolor='w', edgecolor='w', orientation='portrait', format="pdf", bbox_inches='tight')
# plt.show()

'''

            182677 hcco (5.4)ladder(cut) Lx=8 Ly=2 l=200 U=-7.00 mu=-4.0000 h=0.0000   R       1:55    1 node-19
            182678 hcco (5.4)ladder(cut) Lx=8 Ly=2 l=200 U=-7.00 mu=-4.0000 h=1.0000   R       1:55    1 node-18
            182679 hcco (5.4)ladder(cut) Lx=8 Ly=2 l=200 U=-7.00 mu=-4.0000 h=2.0000  PD       0:00    1 (Resources)
            182680 hcco (5.4)ladder(cut) Lx=8 Ly=2 l=200 U=-7.00 mu=-4.5000 h=0.0000  PD       0:00    1 (Priority)
            182681 hcco (5.4)ladder(cut) Lx=8 Ly=2 l=200 U=-7.00 mu=-4.5000 h=1.0000  PD       0:00    1 (Priority)
            182682 hcco (5.4)ladder(cut) Lx=8 Ly=2 l=200 U=-7.00 mu=-4.5000 h=2.0000  PD       0:00    1 (Priority)
            182683 hcco (5.4)ladder(cut) Lx=8 Ly=2 l=200 U=-7.00 mu=-5.0000 h=0.0000  PD       0:00    1 (Priority)
            182684 hcco (5.4)ladder(cut) Lx=8 Ly=2 l=200 U=-7.00 mu=-5.0000 h=1.0000  PD       0:00    1 (Priority)
            182685 hcco (5.4)ladder(cut) Lx=8 Ly=2 l=200 U=-7.00 mu=-5.0000 h=2.0000  PD       0:00    1 (Priority)


            182686 hcco                      (3)ne=15 eta=1.02 V=2.40 W=0.00 realz=0  PD       0:00    1 (Priority)
            182687 hcco                      (3)ne=15 eta=1.02 V=3.20 W=0.00 realz=0  PD       0:00    1 (Priority)
            182688 hcco                      (3)ne=15 eta=1.03 V=2.40 W=0.00 realz=0  PD       0:00    1 (Priority)
            182689 hcco                      (3)ne=15 eta=1.03 V=3.20 W=0.00 realz=0  PD       0:00    1 (Priority)
            182690 hcco                      (3)ne=15 eta=1.07 V=3.20 W=0.00 realz=0  PD       0:00    1 (Priority)
            182691 hcco                      (3)ne=15 eta=1.08 V=3.20 W=0.00 realz=0  PD       0:00    1 (Priority)
            182692 hcco                      (3)ne=15 eta=1.09 V=3.20 W=0.00 realz=0  PD       0:00    1 (Priority)
            182693 hcco                      (3)ne=15 eta=1.04 V=3.20 W=0.00 realz=0  PD       0:00    1 (Priority)
            182694 hcco                      (3)ne=15 eta=1.05 V=3.20 W=0.00 realz=0  PD       0:00    1 (Priority)
            182695 hcco                      (3)ne=15 eta=1.06 V=3.20 W=0.00 realz=0  PD       0:00    1 (Priority)

'''