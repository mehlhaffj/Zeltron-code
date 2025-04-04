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
import math
import matplotlib.pyplot as plt
import sys

def plot_spectrum(*args, **kwargs):
    it = args[0]

    spec = args[1]
    
    spec = "both"
    if "species" in kwargs:
        spec = kwargs["species"]

    #===============================================================================
    # 4-velocity
    u=numpy.loadtxt(".././Zeltron3D/data/u.dat")
    u=u[0:len(u)-1]
    
    # Spectrum
    if spec == "both":
        dNdud =numpy.loadtxt(".././Zeltron3D/data/spectrum_electrons_drift"+it+".dat")
        dNdub =numpy.loadtxt(".././Zeltron3D/data/spectrum_electrons_bg"+it+".dat")
        dNdud+=numpy.loadtxt(".././Zeltron3D/data/spectrum_ions_drift"+it+".dat")
        dNdub+=numpy.loadtxt(".././Zeltron3D/data/spectrum_ions_bg"+it+".dat")
    else:
        dNdud=numpy.loadtxt(".././Zeltron3D/data/spectrum_"+spec+"_drift"+it+".dat")
        dNdub=numpy.loadtxt(".././Zeltron3D/data/spectrum_"+spec+"_bg"+it+".dat")
    
    # Total spectrum
    dNdu=dNdud+dNdub
    
    plt.plot(u,dNdu,color='blue',lw=2)
    plt.xscale('log')
    plt.yscale('log')

    plt.xlabel(r'$\gamma\beta$',fontsize=20)
    plt.ylabel(r'$\frac{dN}{d(\gamma\beta)}$',fontsize=20)
    plt.ylim([1e-3*numpy.max(dNdu),3.0*numpy.max(dNdu)])
    
    plt.title("Time step="+it+", Species="+spec, fontsize=18)
    
    #===============================================================================

    # plt.show()
    fname=".././data/plots/spectrum_%s.png" % (it)
    if "save" in kwargs.keys():
        fname = kwargs["save"]
    plt.savefig(fname, bbox_inches="tight")
    
    #===============================================================================

args, kwargs = pu.getArgsAndKwargs(sys.argv[1:])
if len(args) == 0:
    print("usage: python plot_densities_slice.py timestep [keywords]")
    print("  Plot E&B fields on one 2D plane")
    print("  ix=;iy=;iz=")
    print("  Plane is specied by one cell index in x, y, or z")
    print("  Default is iy=0")
    print("  species = [electrons|ions|both] default is 'both'")
    print("  hmin/hmax = the min/max abscissa axis values")
    print("  vmin/vmax = the min/max ordinate axis values")
    # print("  cmin/cmax = the min/max colorbar values")
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
plot_spectrum(*args, **kwargs)
