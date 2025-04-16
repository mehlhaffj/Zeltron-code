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

import sys
import numpy as np

import vtkutil
import plotutil as pu

args, kwargs = pu.getArgsAndKwargs(sys.argv[1:])
if len(args) == 0:
    print("usage: python h5tovtk.py timestep [keywords]")
    sys.exit(0)
vtkutil.h5tovtk(*args, **kwargs)
