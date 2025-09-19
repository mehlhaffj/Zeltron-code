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

def plot_fields_slice(*args, **kwargs):
    it = args[0]

    simdir   = "./../"
    datadir  = os.path.join(simdir, "data")
    fielddir = os.path.join(datadir, "fields")
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
    xnorm = rlc
    ynorm = rlc
    znorm = rlc

    xlabel = r"$x/R_{\rm LC}$"
    ylabel = r"$y/R_{\rm LC}$"
    zlabel = r"$z/R_{\rm LC}$"

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
        hmin_plot = kwargs["hmin"]

    vmin_plot = vmin/vnorm
    # ymin_plot = ymin/(sigc*rho0)
    if "vmin" in kwargs.keys():
        vmin_plot = kwargs["vmin"]
    vmax_plot = vmax/vnorm
    # ymax_plot = ymax/(sigc*rho0)
    if "vmax" in kwargs.keys():
        vmax_plot = kwargs["vmax"]

    cmap = "RdBu_r"
    if "cmap" in kwargs.keys():
        cmap = kwargs["cmap"]
    cmin = -1.0
    if "cmin" in kwargs.keys():
        cmin = kwargs["cmin"]
    cmax = 1.0
    if "cmax" in kwargs.keys():
        cmax = kwargs["cmax"]
    clog = False
    if "clog" in kwargs.keys():
        clog = kwargs["clog"]
    clinthresh = cmax / 10.0
    if "clinthresh" in kwargs.keys():
        clinthresh = kwargs["clinthresh"]
    clinscale = 0.5
    if "clinscale" in kwargs.keys():
        clinscale = kwargs["clinscale"]

    nrows = 2
    if "nrows" in kwargs.keys():
        nrows = kwargs["nrows"]
    ncols = 3
    if "ncols" in kwargs.keys():
        ncols = kwargs["ncols"]

    if nrows * ncols != 6:
        print("Need to have nrows * ncols == 6.")
        print("Got (nrows, ncols) = (%g, %g)" % (nrows, ncols))
        print("Aborting...")
        sys.exit(0)

    #===========================================================================
    # Read simulation data
    Bz = pu.readArrayFromHdf5(os.path.join(fielddir, "Bz_%s.h5" % it), "field")
    By = pu.readArrayFromHdf5(os.path.join(fielddir, "By_%s.h5" % it), "field")
    Bx = pu.readArrayFromHdf5(os.path.join(fielddir, "Bx_%s.h5" % it), "field")
    Ez = pu.readArrayFromHdf5(os.path.join(fielddir, "Ez_%s.h5" % it), "field")
    Ey = pu.readArrayFromHdf5(os.path.join(fielddir, "Ey_%s.h5" % it), "field")
    Ex = pu.readArrayFromHdf5(os.path.join(fielddir, "Ex_%s.h5" % it), "field")

    Bz = Bz[dslice].squeeze()
    By = By[dslice].squeeze()
    Bx = Bx[dslice].squeeze()
    Ez = Ez[dslice].squeeze()
    Ey = Ey[dslice].squeeze()
    Ex = Ex[dslice].squeeze()

    if need_transpose:
        Bz = np.transpose(Bz)
        By = np.transpose(By)
        Bx = np.transpose(Bx)
        Ez = np.transpose(Ez)
        Ey = np.transpose(Ey)
        Ex = np.transpose(Ex)

    #===========================================================================
    # Build figure
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

    cnorm = colors.Normalize(cmin, cmax)
    if clog and cmin*cmax > 0.0:
        cnorm = colors.LogNorm(cmin, cmax)
    if clog and cmin*cmax < 0.0:
        cnorm = colors.SymLogNorm(
            clinthresh, linscale=clinscale,
            vmin=cmin, vmax=cmax
        )


    #===========================================================================
    # Ex
    ax = grid[0]
    ax.grid(False)
    ax.set_title(titlestr, loc = "left")
    pm = ax.pcolormesh(
        hcoords/hnorm, vcoords/vnorm, Ex/B0,
        cmap = cmap,
        shading = "auto",
        norm = cnorm,
    )
    ax.set_xlabel(hlabel)
    ax.set_ylabel(vlabel)
    ax.set_xlim((hmin_plot, hmax_plot))
    ax.set_ylim((vmin_plot, vmax_plot))
    ax.text(
        texth, textv, r"$E_x/B_0$",
        ha = textha, va = textva,
        path_effects = path_effects,
    )
    ax.set_aspect("equal")

    #===========================================================================
    # Ey
    ax = grid[1]
    ax.grid(False)
    pm = ax.pcolormesh(
        hcoords/hnorm, vcoords/vnorm, Ey/B0,
        cmap = cmap,
        shading = "auto",
        norm = cnorm,
    )
    ax.set_xlabel(hlabel)
    ax.set_ylabel(vlabel)
    ax.set_xlim((hmin_plot, hmax_plot))
    ax.set_ylim((vmin_plot, vmax_plot))
    ax.text(
        texth, textv, r"$E_y/B_0$",
        ha = textha, va = textva,
        path_effects = path_effects,
    )
    ax.set_aspect("equal")

    #===========================================================================
    # Ez
    ax = grid[2]
    ax.grid(False)
    pm = ax.pcolormesh(
        hcoords/hnorm, vcoords/vnorm, Ez/B0,
        cmap = cmap,
        shading = "auto",
        norm = cnorm,
    )
    ax.set_xlabel(hlabel)
    ax.set_ylabel(vlabel)
    ax.set_xlim((hmin_plot, hmax_plot))
    ax.set_ylim((vmin_plot, vmax_plot))
    ax.text(
        texth, textv, r"$E_z/B_0$",
        ha = textha, va = textva,
        path_effects = path_effects,
    )
    ax.set_aspect("equal")
    ax.set_aspect("equal")

    #===========================================================================
    # Bx
    ax = grid[3]
    ax.grid(False)
    pm = ax.pcolormesh(
        hcoords/hnorm, vcoords/vnorm, Bx/B0,
        cmap = cmap,
        shading = "auto",
        norm = cnorm,
    )
    ax.set_xlabel(hlabel)
    ax.set_ylabel(vlabel)
    ax.set_xlim((hmin_plot, hmax_plot))
    ax.set_ylim((vmin_plot, vmax_plot))
    ax.text(
        texth, textv, r"$B_x/B_0$",
        ha = textha, va = textva,
        path_effects = path_effects,
    )
    ax.set_aspect("equal")

    #===========================================================================
    # By
    ax = grid[4]
    ax.grid(False)
    pm = ax.pcolormesh(
        hcoords/hnorm, vcoords/vnorm, By/B0,
        cmap = cmap,
        shading = "auto",
        norm = cnorm,
    )
    ax.set_xlabel(hlabel)
    ax.set_ylabel(vlabel)
    ax.set_xlim((hmin_plot, hmax_plot))
    ax.set_ylim((vmin_plot, vmax_plot))
    ax.text(
        texth, textv, r"$B_y/B_0$",
        ha = textha, va = textva,
        path_effects = path_effects,
    )
    ax.set_aspect("equal")

    #===========================================================================
    # Bz
    ax = grid[5]
    ax.grid(False)
    pm = ax.pcolormesh(
        hcoords/hnorm, vcoords/vnorm, Bz/B0,
        cmap = cmap,
        shading = "auto",
        norm = cnorm,
        path_effects = path_effects,
    )
    ax.set_xlabel(hlabel)
    ax.set_ylabel(vlabel)
    ax.set_xlim((hmin_plot, hmax_plot))
    ax.set_ylim((vmin_plot, vmax_plot))
    ax.text(
        texth, textv, r"$B_z/B_0$",
        ha = textha, va = textva,
        path_effects = path_effects,
    )
    ax.set_aspect("equal")

    #===========================================================================
    # Do the colorbar
    ax.cax.grid(False)
    ax.cax.colorbar(pm)

    #===========================================================================
    fname=".././data/plots/fields_%s_%s.png" % (slicestr, it)
    if "save" in kwargs.keys():
        fname = kwargs["save"]
    plt.savefig(fname, bbox_inches="tight")

args, kwargs = pu.getArgsAndKwargs(sys.argv[1:])
if len(args) == 0:
    print("usage: python plot_fields_slice.py timestep [keywords]")
    print("  Plot E&B fields on one 2D plane")
    print("  ix=;iy=;iz=")
    print("  Plane is specied by one cell index in x, y, or z")
    print("  Default is iy=0")
    print("  hmin/hmax = the min/max abscissa axis values")
    print("  vmin/vmax = the min/max ordinate axis values")
    print("  cmin/cmax = the min/max colorbar values")
    print("  clog = True/False(default) -- whether colorbar is logarithmic")
    print("  clin[thresh/scale] = linthresh/linscale for matplotlib symlognorm")
    print("  save = the file name (including extension) to save the figure to")
    print("  verbose = True/False -- print runtime diagnostic information")
    print("  nrows/ncols number of rows and columns in figure panels")
    print("    Defaults to nrows = 2; ncols = 3. Need to have nrows*ncols = 6.")
    # print("  fieldnorm = normalize E/B fields by this multiple of B0")
    # print("  reduceRes = 1 (no reduction), 2, 3, 4, ...")
    # print("  redResPow2 = 0 (no reduction), 1, 2, ...")
    # print("     reduce resolution in each dimension by specified factor of 2")
    # print("     'smooth' smooths is applied each time (including if redResPow2 = 0)")
    # print("  smooth = 0 (no smoothing), 1, 2, 3, ...")
    # print("  norm = [bunif (default), bcone, bparab]")
    # print("  smooth = 0 (no smoothing), 1, 2, 3, ...")
    sys.exit(0)
plot_fields_slice(*args,**kwargs)
