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

from matplotlib import pyplot as plt
import matplotlib.patheffects as pe
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import ImageGrid 
from scipy.signal import find_peaks
import numpy as np
import copy
import sys
import os

from matplotlib.ticker import AutoLocator, MaxNLocator

import plotutil as pu
from custom_cmaps import register_custom_cmaps

plt.style.use("~/repos/public_zeltron/plot3D/aa_ppt.mplstyle")
register_custom_cmaps(plt)
figsizeHere = np.array([1., 4.]) * np.array(plt.rcParams["figure.figsize"])

c = 299792458e2
e = 4.8032068e-10
me=9.1093897e-28

def plot_densities_slice(*args, **kwargs):
    it = args[0]

    simdir   = "./../"
    datadir  = os.path.join(simdir, "data")
    fielddir = os.path.join(datadir, "fields")
    currentdir = os.path.join(datadir, "currents")
    densitydir = os.path.join(datadir, "densities")
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
    # zmono   = params["z_monopole [cm]"]

    # Derived parameters
    omega = c / rlc
    J0 = B0 * omega / (2 * np.pi)
    ngj = B0 * omega / (2 * np.pi * c * e)

    xnorm = rlc
    ynorm = rlc
    znorm = rlc

    xlabel = r"$x/R_{\rm LC}$"
    ylabel = r"$y/R_{\rm LC}$"
    zlabel = r"$z/R_{\rm LC}$"

    fieldnorm = J0
    fieldnormstr = r"n_{\rm GJ}"

    Nx = int(np.round((xmax-xmin)/dx))
    Ny = int(np.round((ymax-ymin)/dy))
    Nz = int(np.round((zmax-zmin)/dz))

    # Coordinates
    z = pu.readArrayFromHdf5(os.path.join(fielddir, "By_%s.h5" % it), "axis0coords")
    y = pu.readArrayFromHdf5(os.path.join(fielddir, "By_%s.h5" % it), "axis1coords")
    x = pu.readArrayFromHdf5(os.path.join(fielddir, "Bx_%s.h5" % it), "axis2coords")

    #===========================================================================
    # Parse args

    ix = None
    iy = int(np.round(Ny/2))
    iz = None
    slicestr = "iy%d" % iy
    dslice = (slice(None), slice(iy,iy+1), slice(None))
    planelabel = r"%s = %.2f" % (ylabel, y[iy]/znorm)
    hlabel = zlabel
    hnorm  = znorm
    hmax   = zmax
    hmin   = zmin
    hcoords= z
    vlabel = xlabel
    vnorm  = xnorm
    vmin   = xmin
    vmax   = xmax
    vcoords= x
    need_transpose = False

    numixyz = 0
    if "ix" in kwargs.keys():
        ix = kwargs["ix"]
        iy = None
        iz = None
        slicestr = "ix%d" % ix
        dslice = (slice(None), slice(None), slice(ix,ix+1))
        planelabel = r"%s = %.2f" % (xlabel, x[ix]/znorm)
        hlabel = zlabel
        hnorm  = znorm
        hmax   = zmax
        hmin   = zmin
        hcoords= z
        vlabel = ylabel
        vnorm  = ynorm
        vmin   = ymin
        vmax   = ymax
        vcoords= y
        need_transpose = True
        numixyz += 1
    if "iy" in kwargs.keys():
        iy = kwargs["iy"]
        ix = None
        iz = None
        slicestr = "iy%d" % iy
        dslice = (slice(None), slice(iy,iy+1), slice(None))
        planelabel = r"%s = %.2f" % (ylabel, y[iy]/znorm)
        hlabel = zlabel
        hnorm  = znorm
        hmax   = zmax
        hmin   = zmin
        hcoords= z
        vlabel = xlabel
        vnorm  = xnorm
        vmin   = xmin
        vmax   = xmax
        vcoords= x
        need_transpose = True
        numixyz += 1
    if "iz" in kwargs.keys():
        iz = kwargs["iz"]
        ix = None
        iy = None
        slicestr = "iz%d" % iz
        dslice = (slice(iz,iz+1), slice(None), slice(None))
        planelabel = r"%s = %.2f" % (zlabel, z[iz]/znorm)
        hlabel = xlabel
        hnorm  = xnorm
        hmax   = xmax
        hmin   = xmin
        hcoords= x
        vlabel = ylabel
        vnorm  = ynorm
        vmin   = ymin
        vmax   = ymax
        vcoords= y
        need_transpose = False
        numixyz += 1

    if numixyz > 1:
        print("You may only specify one of 'ix', 'iy', or 'iz' args.")
        print("Got %d" % numixyz)
        print("Aborting...")
        sys.exit(0)

    hmax_plot = hmax/hnorm
    # xmax_plot = xmax/(sigc*rho0)
    if "hmax" in kwargs.keys():
        hmax_plot = kwargs["hmax"]
    hmin_plot = hmin/hnorm
    # xmin_plot = xmin/(sigc*rho0)
    if "hmin" in kwargs.keys():
        hmin_plot = kwargs["hmax"]

    vmin_plot = vmin/vnorm
    # ymin_plot = ymin/(sigc*rho0)
    if "vmin" in kwargs.keys():
        vmin_plot = kwargs["vmin"]
    vmax_plot = vmax/vnorm
    # ymax_plot = ymax/(sigc*rho0)
    if "vmax" in kwargs.keys():
        vmax_plot = kwargs["vmax"]

    cmap = "inferno"
    if "cmap" in kwargs.keys():
        cmap = kwargs["cmap"]
    cmin = 0.01
    if "cmin" in kwargs.keys():
        cmin = kwargs["cmin"]
    cmax = 10.0
    if "cmax" in kwargs.keys():
        cmax = kwargs["cmax"]

    #===========================================================================
    # Read simulation data
    nbge = pu.readArrayFromHdf5(os.path.join(densitydir, "mapxyz_electrons_bg_%s.h5" % it), "field")
    nbgi = pu.readArrayFromHdf5(os.path.join(densitydir, "mapxyz_ions_bg_%s.h5" % it), "field")
    nje  = pu.readArrayFromHdf5(os.path.join(densitydir, "mapxyz_electrons_drift_%s.h5" % it), "field")
    nji  = pu.readArrayFromHdf5(os.path.join(densitydir, "mapxyz_ions_drift_%s.h5" % it), "field")

    nbg = nbge  + nbgi
    nj  = nje + nji

    nj  = nj[dslice].squeeze()
    nbg = nbg[dslice].squeeze()

    if need_transpose:
        nj  = np.transpose(nj)
        nbg = np.transpose(nbg)

    #===========================================================================
    # Build figure
    nrows = 3
    ncols = 1
    cbar_location = "top"
    cbar_mode = "single"
    cbar_size = "7%"
    cbar_pad = 0.6

    pad = 0.05
    delta = pad * (hmax_plot - hmin_plot)
    texth = hmax_plot - delta
    textv = vmax_plot - delta
    textha = "right"
    textva = "top"
    lw = plt.rcParams["lines.linewidth"]
    path_effects = [
        pe.Stroke(linewidth = 1.5*lw, foreground = 'w', alpha=1.0),
        pe.Normal()
    ]

    itpad = '{:>6}'.format(it)
    timestr = "timestep="+itpad+r"; $ct/R_{\rm LC}$=%.2f" \
        % (c*dt*int(it)/rlc)
    titlestr = timestr + "; " + planelabel

    fig = plt.figure(1, figsize=figsizeHere)
    grid = ImageGrid(
        fig, 111,
        nrows_ncols = (nrows, ncols),
        axes_pad = 0.0,
        cbar_location = cbar_location,
        cbar_mode = cbar_mode,
        cbar_size = cbar_size,
        cbar_pad = cbar_pad,
    )

    norm = colors.LogNorm(vmin=cmin, vmax=cmax)

    #===========================================================================
    # nj
    ax = grid[0]
    ax.grid(False)
    ax.set_title(titlestr, loc = "left")
    pm = ax.pcolormesh(
        hcoords/hnorm, vcoords/vnorm, nj/fieldnorm,
        cmap = cmap,
        shading = "auto",
        norm = norm,
        # vmin = cmin, vmax = cmax
    )
    ax.set_xlabel(hlabel)
    ax.set_ylabel(vlabel)
    ax.set_xlim((hmin_plot, hmax_plot))
    ax.set_ylim((vmin_plot, vmax_plot))
    ax.text(
        texth, textv, r"$n_{\rm j}/%s$" % fieldnormstr,
        ha = textha, va = textva,
        path_effects = path_effects,
    )
    ax.set_aspect("equal")

    #===========================================================================
    # nbg
    ax = grid[1]
    ax.grid(False)
    pm = ax.pcolormesh(
        hcoords/hnorm, vcoords/vnorm, nbg/fieldnorm,
        cmap = cmap,
        shading = "auto",
        norm = norm,
        # vmin = cmin, vmax = cmax
    )
    ax.set_xlabel(hlabel)
    ax.set_ylabel(vlabel)
    ax.set_xlim((hmin_plot, hmax_plot))
    ax.set_ylim((vmin_plot, vmax_plot))
    ax.text(
        texth, textv, r"$n_{\rm bg}/%s$" % fieldnormstr,
        ha = textha, va = textva,
        path_effects = path_effects,
    )
    ax.set_aspect("equal")

    #===========================================================================
    # nbg + nj
    ax = grid[2]
    ax.grid(False)
    pm = ax.pcolormesh(
        hcoords/hnorm, vcoords/vnorm, (nbg + nj)/fieldnorm,
        cmap = cmap,
        shading = "auto",
        norm = norm,
        # vmin = cmin, vmax = cmax
    )
    ax.set_xlabel(hlabel)
    ax.set_ylabel(vlabel)
    ax.set_xlim((hmin_plot, hmax_plot))
    ax.set_ylim((vmin_plot, vmax_plot))
    ax.text(
        texth, textv, r"$n_{\rm tot}/%s$" % fieldnormstr,
        ha = textha, va = textva,
        path_effects = path_effects,
    )
    ax.set_aspect("equal")

    #===========================================================================
    # Do the colorbar
    ax.cax.grid(False)
    ax.cax.colorbar(pm)

    #===========================================================================
    fname=".././data/plots/densities_%s_%s.png" % (slicestr, it)
    if "save" in kwargs.keys():
        fname = kwargs["save"]
    plt.savefig(fname, bbox_inches="tight")

args, kwargs = pu.getArgsAndKwargs(sys.argv[1:])
if len(args) == 0:
    print("usage: python plot_densities_slice.py timestep [keywords]")
    print("  Plot E&B fields on one 2D plane")
    print("  ix=;iy=;iz=")
    print("  Plane is specied by one cell index in x, y, or z")
    print("  Default is iy=0")
    print("  hmin/hmax = the min/max abscissa axis values")
    print("  vmin/vmax = the min/max ordinate axis values")
    print("  cmin/cmax = the min/max colorbar values")
    print("  save = the file name (including extension) to save the figure to")
    print("  verbose = True/False -- print runtime diagnostic information")
    # print("  fieldnorm = normalize E/B fields by this multiple of B0")
    # print("  reduceRes = 1 (no reduction), 2, 3, 4, ...")
    # print("  redResPow2 = 0 (no reduction), 1, 2, ...")
    # print("     reduce resolution in each dimension by specified factor of 2")
    # print("     'smooth' smooths is applied each time (including if redResPow2 = 0)")
    # print("  smooth = 0 (no smoothing), 1, 2, 3, ...")
    # print("  norm = [bunif (default), bcone, bparab]")
    # print("  smooth = 0 (no smoothing), 1, 2, 3, ...")
    sys.exit(0)
plot_densities_slice(*args,**kwargs)
