#***********************************************************************!
#                       The Zeltron code project.                       !
#***********************************************************************!
# Copyright (C) 2012-2015. Authors: Benoit Cerutti & Greg Werner        !
#                                                                       !
# This program is free software: you can redistribute it and/or modify  !
# it under the terms of the GNU General Public License as published by  !
# the Free Software Foundation, either version 3 of the License, or     !
# (at your option) any later version.                                   !
#                                                                       !
# This program is distributed in the hope that it will be useful,       !
# but WITHOUT ANY WARRANTY; without even the implied warranty of        !
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !
# GNU General Public License for more details.                          !
#                                                                       !
# You should have received a copy of the GNU General Public License     !
# along with this program. If not, see <http://www.gnu.org/licenses/>.  !
#***********************************************************************!
#
# This python script will draw the total electron or ion energy spectrum:
# dN/du(u), where u is the 4 velocity, at a given time step it.
#
# To execute, type for instance:$ python plot_spectrum.py 0 ions &
# This command will draw the ions spectrum at time step 0.
#
#***********************************************************************!

import matplotlib
import numpy
import scipy
import math
import matplotlib.pyplot as plt
import sys
np = numpy

import plotutil as pu
from custom_cmaps import register_custom_cmaps

plt.style.use("~/repos/public_zeltron/plot3D/aa_ppt.mplstyle")
register_custom_cmaps(plt)
figsizeHere = np.array([1., 4.]) * np.array(plt.rcParams["figure.figsize"])

c = 299792458e2
e = 4.8032068e-10
me=9.1093897e-28

def plot_spectrum(*args, **kwargs):
    it = args[0]

    simdir   = "./../"
    params   = pu.get_sim_params_dict(simdir)

    #===========================================================================
    # Parameters of the simulation
    xmax = params["xmax [cm]"]
    xmin = params["xmin [cm]"]
    ymax = params["ymax [cm]"]
    ymin = params["ymin [cm]"]
    zmax = params["zmax [cm]"]
    zmin = params["zmin [cm]"]
    dt   = params["dt [s]"]
    dx   = params["dx [cm]"]
    dy   = params["dy [cm]"]
    dz   = params["dz [cm]"]

    rlc     = params["R_LC [cm]"]
    rnozzle = params["R_nozzle [cm]"]
    sigmabg = params["sigma"]
    B0      = params["B0 [G]"]

    #===========================================================================
    # Parse args

    spec = "both"
    if "species" in kwargs:
        spec = kwargs["species"]

    ymax_plot = None
    if "ymax" in kwargs.keys():
        ymax_plot = kwargs["ymax"]
    ymin_plot = None
    if "ymin" in kwargs.keys():
        ymin_plot = kwargs["ymin"]

    xmin_plot = None
    if "xmin" in kwargs.keys():
        xmin_plot = kwargs["xmin"]
    xmax_plot = None
    if "xmax" in kwargs.keys():
        xmax_plot = kwargs["xmax"]

    ynorm = 1.0
    if "ynorm" in kwargs.keys():
        ynorm = kwargs["ynorm"]
    
    compfac = 2
    if "compfac" in kwargs.keys():
        compfac = kwargs["compfac"]

    legend_loc = "best"
    if "legend_loc" in kwargs.keys():
        legend_loc = kwargs["legend_loc"]

    #===============================================================================
    # 4-velocity
    u=numpy.loadtxt(".././data/u.dat")
    u=u[0:len(u)-1]
    
    #===============================================================================
    # Spectrum
    if spec == "both":
        dNdud =pu.readSpectrumFromHdf5(".././data/spectrum_electrons_drift_"+it+".h5")
        dNdub =pu.readSpectrumFromHdf5(".././data/spectrum_electrons_bg_"+it+".h5")
        dNdud+=pu.readSpectrumFromHdf5(".././data/spectrum_ions_drift_"+it+".h5")
        dNdub+=pu.readSpectrumFromHdf5(".././data/spectrum_ions_bg_"+it+".h5")
    else:
        dNdud=pu.readSpectrumFromHdf5(".././data/spectrum_"+spec+"_drift_"+it+".h5")
        dNdub=pu.readSpectrumFromHdf5(".././data/spectrum_"+spec+"_bg_"+it+".h5")
    
    # Total spectrum
    dNdu=dNdud+dNdub

    #===============================================================================
    # Build figure

    xvals = u
    xstr = r"u"
    xlabel = r"$u = \gamma \beta$"

    compstr = r"%s^%g" % (xstr, compfac)
    if compfac == 1:
        compstr = r"%s" % xstr
    if compfac == 0:
        compstr = r""

    yvals_list = [
        xvals**compfac * dNdud,
        xvals**compfac * dNdub,
        xvals**compfac * dNdu ,
    ]
    yvals_labels = [
        "inj.",
        "bg.",
        "tot.",
    ]

    ynorms = [ynorm] * len(yvals_list)
    if ynorm == "auto":
        for i in range(len(yvals_list)):
            ynorms[i] = scipy.integrate.simpson(yvals_list[i], x=xvals)
            if ynorms[i] == 0.0:
                ynorms[i] == 1.0

    maxval = np.max(yvals_list[0] / ynorms[0])
    if len(yvals_list) > 1:
        for i, yvals in enumerate(yvals_list[1:]):
            maxval = max(maxval, np.max(yvals / ynorms[i]))

    if ymin_plot is None:
        ymin_plot = 1e-3 * maxval
    if ymax_plot is None:
        ymax_plot = 3.0 * maxval

    ylabel = r"$%s \frac{dN}{d%s}$" % (compstr, xstr)

    itpad = '{:>6}'.format(it)
    timestr = "timestep="+itpad+r"; $ct/R_{\rm LC}$=%.2f" \
        % (c*dt*int(it)/rlc)
    titlestr = timestr

    fig = plt.figure()
    ax = fig.gca()
    ax.semilogx()
    ax.semilogy()
    
    # ax.plot(u,dNdu,color='blue',lw=2)
    for i, yvals in enumerate(yvals_list):
        ax.plot(
            xvals, yvals / ynorms[i],
            label = yvals_labels[i]
        )

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(titlestr)
    ax.set_xlim(xmin_plot, xmax_plot)
    ax.set_ylim(ymin_plot, ymax_plot)
    if len(yvals_list) > 1:
        ax.legend(loc = legend_loc)
    
    #===============================================================================

    # plt.show()
    fname=".././data/plots/spectrum_%s.png" % (it)
    if "save" in kwargs.keys():
        fname = kwargs["save"]
    plt.savefig(fname, bbox_inches="tight")
    
    #===============================================================================

args, kwargs = pu.getArgsAndKwargs(sys.argv[1:])
if len(args) == 0:
    print("usage: python plot_spectrum.py timestep [keywords]")
    print("  species = [electrons|ions|both] default is 'both'")
    print("  ymin/ymax = the min/max abscissa axis values")
    print("  xmin/xmax = the min/max ordinate axis values")
    # print("  cmin/cmax = the min/max colorbar values")
    print("  save = the file name (including extension) to save the figure to")
    print("  verbose = True/False -- print runtime diagnostic information")
    print("  compfac = 2(default) -- compensate factor of distribution")
    print("  legend_loc = 'best'(default) -- where to place matplotlib legend.")
    print("  ynorm = 1.0(default) -- normalize all abscissa by this number.")
    print("    'auto' normalizes all spectra by their integral.")
    # print("  fieldnorm = normalize E/B fields by this multiple of B0")
    # print("  reduceRes = 1 (no reduction), 2, 3, 4, ...")
    # print("  redResPow2 = 0 (no reduction), 1, 2, ...")
    # print("     reduce resolution in each dimension by specified factor of 2")
    # print("     'smooth' smooths is applied each time (including if redResPow2 = 0)")
    # print("  smooth = 0 (no smoothing), 1, 2, 3, ...")
    # print("  norm = [bunif (default), bcone, bparab]")
    # print("  smooth = 0 (no smoothing), 1, 2, 3, ...")
    sys.exit(0)
plot_spectrum(*args, **kwargs)
