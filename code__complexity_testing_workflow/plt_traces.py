import matplotlib.pyplot as plt
import numpy as np



def plt_six_traces(dd_rdx, plotting_window_0, plotting_window, sample_rate, simulation_location):

    times = (np.arange(len(dd_rdx[0,plotting_window_0:plotting_window_0+plotting_window])) + plotting_window_0)/sample_rate

    fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6, 1, figsize=(10, 30), sharex=True, sharey=True)

    ax1.plot(times, dd_rdx[0,plotting_window_0:plotting_window_0+plotting_window])
    ax1.set_title("Soma Signal")
    ax1.margins(0, .1)
    ax1.grid(alpha=.5, ls='--')

    ax2.plot(times, dd_rdx[1,plotting_window_0:plotting_window_0+plotting_window])
    ax2.set_title("den-1 Signal")
    ax2.margins(0, .1)
    ax2.grid(alpha=.5, ls='--')

    ax3.plot(times, dd_rdx[2,plotting_window_0:plotting_window_0+plotting_window])
    ax3.set_title("den-2 Signal")
    ax3.margins(0, .1)
    ax3.grid(alpha=.5, ls='--')

    ax4.plot(times, dd_rdx[3,plotting_window_0:plotting_window_0+plotting_window])
    ax4.set_title("den-3 Signal")
    ax4.margins(0, .1)
    ax4.grid(alpha=.5, ls='--')

    ax5.plot(times, dd_rdx[4,plotting_window_0:plotting_window_0+plotting_window])
    ax5.set_title("den-4 Signal")
    ax5.margins(0, .1)
    ax5.grid(alpha=.5, ls='--')

    ax6.plot(times, dd_rdx[4,plotting_window_0:plotting_window_0+plotting_window])
    ax6.set_title("den-5 Signal")
    ax6.margins(0, .1)
    ax6.grid(alpha=.5, ls='--')


    plt.tight_layout()

    plt.savefig( simulation_location + '/six_trace.png')



def plt_three_traces(dd, dd_lp, dd_hp, plotting_window_0, plotting_window, sample_rate, filter_f, simulation_location ):

    times = (np.arange(len(dd[1,plotting_window_0:plotting_window_0+plotting_window])) + plotting_window_0)/sample_rate

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(30, 10), sharex=True, sharey=True)
    ax1.plot(times, dd[0,plotting_window_0:plotting_window_0+plotting_window])
    ax1.set_title("Original Signal")
    ax1.margins(0, .1)
    ax1.grid(alpha=.5, ls='--')
    ax2.plot(times, dd_lp[0,plotting_window_0:plotting_window_0+plotting_window])
    ax2.set_title("Low-Pass Filter (" + str(filter_f*1000) + "Hz)")
    ax2.grid(alpha=.5, ls='--')
    ax3.plot(times, dd_hp[0,plotting_window_0:plotting_window_0+plotting_window])
    ax3.set_title("High-Pass Filter (" + str(filter_f*1000) + "Hz)")
    ax3.grid(alpha=.5, ls='--')
    plt.tight_layout()

    plt.savefig( simulation_location + '/hi_lo_read.png')
