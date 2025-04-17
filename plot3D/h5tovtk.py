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

import sys
import numpy as np

import vtkutil
import plotutil as pu

args, kwargs = pu.getArgsAndKwargs(sys.argv[1:])
if len(args) == 0:
    print("usage: python h5tovtk.py timestep [keywords]")
    print("  redResPow2 = 0 (no reduction), 1, 2, ...")
    print("     reduce resolution in each dimension by specified factor of 2")
    print("     nearest-neighbor averaged is performed for each power of 2")
    print("     so that reduced-res grids hold the average over neigboring")
    print("     cells in their parent grids")
    print("     redResPow2 will only work for (reduced) grid sizes of powers")
    print("     of 2.")
    sys.exit(0)
# vtkutil.h5tovtk_sgrid(*args, **kwargs)
# vtkutil.h5tovtk_imagedata(*args, **kwargs)
vtkutil.h5tovtk_rgrid(*args, **kwargs)
