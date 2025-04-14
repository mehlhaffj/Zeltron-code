import numpy
import os
import tables
import h5py
np=numpy

# # Class to import custom XML colormaps, of the kind you can download from
# # https://sciviscolor.org/colormaps/divergent/
# def cmap_from_xml(xml_cmap_file, cmap_name = "mycmap"):
#     try:
#         from bs4 import BeautifulSoup
#     except:
#         msg = "func cmap_from_xml requires missing module: "
#         msg += "bs4.BeautifulSoup"
#         raise ValueError(msg)
# 
#     try:
#         from matplotlib.colors import LinearSegmentedColormap
#     except:
#         raise ValueError(msg)
#         msg = "func cmap_from_xml requires missing module: "
#         msg += "matplotlib.colors.LinearSegmentedColormap"
# 
#     try:
#         from pathlib import Path
#     except:
#         msg = "func cmap_from_xml requires missing module: "
#         msg = "pathlib.Path; are you running on Python 3?"
#         raise ValueError(msg)
# 
#     path = Path(xml_cmap_file)
#     txt = path.read_text()
#     soup = BeautifulSoup(txt, "html.parser")
#     points = soup.find_all("point")
# 
#     levels = np.array( [ float(points[i].get("x")) for i in range(len(points))] )
#     colors = np.array( [[float(points[i].get("r")), float(points[i].get("g")),
#                         float(points[i].get("b"))] for i in range(len(points))] )
# 
#     return LinearSegmentedColormap.from_list(cmap_name, list(zip(levels, colors)))
 

# Stolen from Greg Werner
def getArgsAndKwargs(allArgs): #{
    """
    Take a list of strings (e.g., command line arguments)
    and split them into args and kwargs -- a kwarg has an equals sign in
    it, e.g., "name=val".  Returns the args as they are, in a list, and
    return a dictionary of kwargs, converting the value through evaluation.
    """
    args = []
    kwargs = dict()
    for a in allArgs: #{
        if '=' in a: #{
            parts = a.split("=")
            k = parts[0]
            vStr = '='.join(parts[1:])
            #k,vStr = a.split("=")
            try: #{
                # don't let local variables (like k and a and vStr) get evaluated
                #v = eval(vStr)
                v = eval(vStr, dict(), dict())
            except: #}{
                v = vStr
                #}
            kwargs[k] = v
        else: #}{
            args.append(a)
        #}
    #}
    return (args, kwargs)
#}


def get_sim_params_dict(top_sim_dir):
    directory = top_sim_dir
    inputp_file = os.path.join(directory, "data", "input_params.dat")
    physp_file  = os.path.join(directory, "data", "phys_params.dat")
    physp = open(physp_file, 'r')
    physl = physp.readlines()
    physp.close()
    inputp = open(inputp_file, 'r')
    inputl = inputp.readlines()
    inputp.close()
    param_names = "  ".join([physl[0], inputl[0]])
    param_vals  = "  ".join([physl[1], inputl[1]])
    param_names = param_names.split("  ")
    param_vals  = param_vals.split("  ")

    names_trim = []
    for name in param_names:
        if len(name) > 0 and not name.isspace():
            names_trim.append(name.strip())

    vals_trim = []
    for val in param_vals:
        if len(val) > 0 and not val.isspace():
            vals_trim.append(float(val.strip()))

    pdict = dict(list(zip(names_trim, vals_trim)))

    if "FDUMP_BASE" not in pdict and "FDUMP" in pdict:
        pdict["FDUMP_BASE"] = pdict["FDUMP"]

    return pdict


def readArrayFromHdf5(h5FileName, datasetName, slices = None,
  shapeOnly = False, verbose=False, alternateDatasetName = None): #{
  """
  Returns a (numpy float64) array of the field values
  from an hdf5 dataset

  If slices is set, this will retrieve the requested slice of a field;
  e.g., if the field is 4x5x3, then 
  slices = (slice(None), slice(2,5), slice(None)) will return a 
  4x3x3 subset.

  If shapeOnly, returns the field shape rather than the field, and 
    slices is ignored.
  """
  try:
    file = tables.open_file(h5FileName,'r')
    getNode = file.get_node
  except:
    try:
      file = tables.openFile(h5FileName,'r')
      getNode = file.getNode
    except:
      msg = "\n\nError: Failed to open hdf5 file: '" + h5FileName + "'"
      msg += " in " + os.getcwd() + " (is it an hdf5 file?)\n"
      raise ValueError(msg)

  groupName = '/' + datasetName

  if alternateDatasetName is not None and datasetName not in getNode("/"):
    if alternateDatasetName not in getNode("/"):
      msg = "Neither " + datasetName + " nor " + alternateDatasetName 
      msg += " names a datset in " + h5FileName
      raise ValueError(msg)
    groupName = '/' + alternateDatasetName

  try:
    if shapeOnly: #{
      res = getNode(groupName).shape
    else: #}{
      if (slices is None):
        res = numpy.array( getNode(groupName) )
      else:
        if verbose:
          print("Getting slices =", slices, "for", h5FileName)
        res = numpy.array( getNode(groupName)[tuple(slices)] )
        #print "Got slices =", slices, "->", res.shape
    #}
  except:
    file.close()
    raise
    msg = "\n\nError: Failed to open dataset '" + groupName + "'"
    msg += " in hdf5 file: '" + h5FileName + "'\n"
    raise ValueError(msg)
  file.close()

  return res
#}


def readSpectrumFromHdf5(ifile):
    h5file = h5py.File(ifile, 'r')
    output = h5file["field"][...]
    h5file.close()

    return output


class FieldNormalizer(object):
    allowed_norm_modes = ["bunif", "bcone", "bparab", "bnozzle"]

    def __init__(self, norm_mode, Bz0, z, r):
        if not (norm_mode in self.allowed_norm_modes):
            msg = "requsted normalization '%s' " % norm_mode
            msg += "not one of allowed avlues: %s" % str(self.allowed_norm_modes)
            raise ValueError(msg)
        
        self.mode = norm_mode

        if norm_mode == "bcone":
            # Extract the location of the monopole generating
            # the conical field from the magnetic field at time 0
            rat=np.sqrt(Bz0[1,0]/Bz0[0,0])
            zmp=(rat*z[1]-z[0])/(rat-1)

            self.z = zmp

        elif norm_mode == "bparab":
            # Extract the location of the equatorial toroidal current sheet
            # (which lives at some zeq<zmin).
            rat=np.sqrt(Bz0[1,0]/Bz0[0,0])
            z0=z[0]
            z1=z[1]
            r0=r[0]
            zeq=1./(1.-rat**2)*\
                (z0-rat*z1-
                 np.sqrt((z0-rat*z1)**2-
                         (1-rat**2)*(r0**2*(1.-rat**2)+z0**2-rat**2*z1**2)))

            self.z = zeq

        elif norm_mode == "bnozzle":
            # In this case, we are dealing with the superposition of a
            # monopolar and a uniform field.
            # The 3 underlying parameters (monopole location, monopole
            # strength, uniform field strength) can be deduced from the
            # initial Bz field after some terrible algebra.
            eval1 = Bz0[0,0]
            eval2 = Bz0[1,0]
            eval3 = Bz0[2,0]
            z1 = z[0]
            z2 = z[1]
            z3 = z[2]
            A = (eval1-eval2)/(eval1-eval3)

            zroots = np.roots(self._bnozzle_zmp_coeffs(z1,z2,z3,A))
            # Eliminate imaginary roots
            zroots = zroots[np.where(np.isreal(zroots))]
            # Eliminate roots > zmin
            zroots = zroots[np.where(zroots <= z[0])]

            # There should be just one feasible root
            # self.z is the monopole location
            self.z = zroots[0]
            zmp = self.z

            # self.b0 is the monopole strength
            self.b0 = (eval1-eval2)/( 1.0 - ((z1-zmp)/(z2-zmp))**2 )
            # self.b1 is the uniform field strength
            self.b1 = eval1 - self.b0

    @staticmethod
    def _bnozzle_zmp_coeffs(z1,z2,z3,A):
        # Used sympy to help derive these expressions
        cpow3 = -2*A*z1 + 2*A*z3 + 2*z1 - 2*z2
        cpow2 = A*z1**2 + 4*A*z1*z2 - 4*A*z2*z3 - A*z3**2 - z1**2 - 4*z1*z3 + z2**2 + 4*z2*z3
        cpow1 = -2*A*z1**2*z2 - 2*A*z1*z2**2 + 2*A*z2**2*z3 + 2*A*z2*z3**2 + 2*z1**2*z3 + 2*z1*z3**2 - 2*z2**2*z3 - 2*z2*z3**2
        cpow0 = A*z1**2*z2**2 - A*z2**2*z3**2 - z1**2*z3**2 + z2**2*z3**2
        return [cpow3,cpow2,cpow1,cpow0]

    @staticmethod
    def _bparab_normer(z, r, zeq):
        result = np.zeros(z.shape)
        idx1 = np.where(r < 0.01*(z-zeq))
        idx2 = np.where(r >=0.01*(z-zeq))
        r1 = r[idx1]
        z1 = z[idx1]
        r2 = r[idx2]
        z2 = z[idx2]
        result[idx1] = 1./(r1**2+(z1-zeq)**2) 
        result[idx2] = 2./(r2**2+(z2-zeq)**2)*\
            (1.+((z2-zeq)/r2)**2-(z2-zeq)/r2*np.sqrt(1.+((z2-zeq)/r2)**2))

        return result

    def field_comp(self, z, r, zmin = 0.0, rmin = 0.0):
        R, Z = _check_input_arrays(r, z)
        if self.mode == "bcone":
            zmp = self.z
            field_comp = (R**2+(Z-zmp)**2)/((zmin-zmp)**2)

        elif self.mode == "bparab":
            zeq = self.z
            numer = self._bparab_normer(zmin+0.0*R, rmin+0.0*R, zeq)
            denom = self._bparab_normer(Z, R, zeq)
            field_comp = np.sqrt(numer / denom)

        elif self.mode == "bnozzle":
            zmp = self.z
            b0 = self.b0
            b1 = self.b1
            fid_field_bz = b1 + b0*(zmin-zmp)**2/(R**2+(Z-zmp)**2) * (Z-zmp)/np.sqrt(R**2+(Z-zmp)**2)
            fid_field_br = b0*(zmin-zmp)**2/(R**2+(Z-zmp)**2) * R/np.sqrt(R**2+(Z-zmp)**2)
            fid_field_str = np.sqrt(fid_field_bz**2 + fid_field_br**2)
            field_comp = b0/fid_field_str

        else:
            field_comp = np.ones(R.shape)

        return field_comp

    def compstr(self):
        if self.mode == "bcone":
            compstr = r"(r/r_0)^2"
        elif self.mode == "bparab" or self.mode == "bnozzle":
            compstr = r"[norm.]"
        else:
            compstr = ""
        return compstr

    def compstr2(self):
        if self.mode == "bcone":
            compstr2 = r"(r/r_0)^4"
        elif self.mode == "bparab" or self.mode == "bnozzle":
            compstr2 = r"[norm.]"
        else:
            compstr2 = ""
        return compstr2


# Smooth a field in 2d according to the filter matrix.
#     |1 2 1|
# M = |2 4 2|
#     |1 2 1|
# pdim1, pdim2 -- Whether the first, second dimensions of 'field' are periodic,
#   respectively.
def smoothField2d(field, pdim1 = False, pdim2 = False):
    expanded = np.zeros( (field.shape[0]+2, field.shape[1]+2) )
    normalization = 16.*np.ones(field.shape)

    expanded[1:-1,1:-1] = field
    if pdim1:
        expanded[0,1:-1] = field[-1,:]
        expanded[-1,1:-1] = field[0,:]
    else:
        normalization[0,:] = 12.
        normalization[-1,:] = 12.

    if pdim2:
        expanded[1:-1,0] = field[:,-1]
        expanded[1:-1,-1] = field[:,0]
    else:
        normalization[:,0] = 12.
        normalization[:,-1] = 12.

    if pdim1 and pdim2:
        # Fill in the corners according to PBCs
        expanded[0,0] = field[-1,-1]
        expanded[0,-1] = field[-1,0]
        expanded[-1,0] = field[0,-1]
        expanded[-1,-1] = field[0,0]

    # If only one of the dimensions is periodic, the corners are already
    # correctly set to 16 or 12. Only if BOTH dimensions are not periodic do we
    # set the normalization at the corners to 9.
    if not pdim1 and not pdim2:
        normalization[0,0] = 9.
        normalization[0,-1] = 9.
        normalization[-1,0] = 9.
        normalization[-1,-1] = 9.

    return 1./normalization * (1.0*expanded[:-2,:-2] + 2.0*expanded[1:-1,:-2]  + 1.0*expanded[2:,:-2]
                             +  2.0*expanded[:-2,1:-1]+ 4.0*expanded[1:-1,1:-1] + 2.0*expanded[2:,1:-1]
                             +  1.0*expanded[:-2,2:]  + 2.0*expanded[1:-1,2:]   + 1.0*expanded[2:,2:])


# Take a cartesian derivative in 2d according to the sobel filter.
#     |-1 0 1|
# M = |-2 0 2|
#     |-1 0 1|
# x,y can be...
#    - 1d coordinate arrays such that data[i,j] is known at (x[i],y[j])
#    - 2d coordinate arrays such that data[i,j] is known at (x[i,j],y[i,j]).
# The derivative is known at the grid points specified by x,y
# The returned array has the same shape as the original data array, but
# contains zeros at the boundary where the derivative cannot be computed.
def sobelx(x,y,data):
    X, Y = _check_input_arrays(x,y)
    output = np.zeros(data.shape)
    output[1:-1,1:-1] = \
        0.25 * (data[2:  ,2:] - data[2:  ,:-2]) / (X[2:  ,2:] - X[2:  ,:-2]) + \
        0.50 * (data[1:-1,2:] - data[1:-1,:-2]) / (X[1:-1,2:] - X[1:-1,:-2]) + \
        0.25 * (data[:-2 ,2:] - data[:-2 ,:-2]) / (X[:-2 ,2:] - X[:-2 ,:-2]) 
    return output


# Take a cartesian derivative in 2d according to the sobel filter.
#     | 1  2  1|
# M = | 0  0  0|
#     |-1 -2 -1|
# x,y can be...
#    - 1d coordinate arrays such that data[i,j] is known at (x[i],y[j])
#    - 2d coordinate arrays such that data[i,j] is known at (x[i,j],y[i,j]).
# The derivative is known at the grid points specified by x,y
# The returned array has the same shape as the original data array, but
# contains zeros at the boundary where the derivative cannot be computed.
def sobely(x,y,data):
    X, Y = _check_input_arrays(x,y)
    output = np.zeros(data.shape)
    output[1:-1,1:-1] = \
        0.25 * (data[2:,2:  ] - data[:-2  ,2:  ]) / (Y[2:  ,2:  ] - Y[:-2  ,2:  ]) + \
        0.50 * (data[2:,1:-1] - data[:-2  ,1:-1]) / (Y[2:  ,1:-1] - Y[:-2  ,1:-1]) + \
        0.25 * (data[2:,:-2 ] - data[:-2  ,:-2 ]) / (Y[2:  ,:-2 ] - Y[:-2  ,:-2 ])
    return output


# Find local minima, maxima, and saddle points of the 2D array data
# x,y can be...
#    - 1d coordinate arrays such that data[i,j] is known at (x[i],y[j])
#    - 2d coordinate arrays such that data[i,j] is known at (x[i,j],y[i,j]).
#
# x,y are used here instead of r,z to emphasize that this is a cartesian
# algorithm, though I expect it will work fine in cylindrical coordinates
# r,z
#
# frac_tol - Fractional tolerance for the second derivative test.
#     The test is determined to have failed if the derivatives at the
#     location of the extremum don't give fractional differences in the
#     data above frac_tol over one grid cell
#     JM: I have found that playing with this doesn't improve results
#     much, and so suggest leaving it equal to zero.
# presmooths - How many times to smooth 'data' before searching for
#     extrema. Smoothing a little bit can have the advantageous effect
#     of reducing the number of spurious extrema found, while smoothing
#     too much can have the detrimental effect of washing out gentle
#     local extrema.
#
# Returns a list of elements of the form (xex,yex,type) where xex and yex are,
# respectively, the x and y coordinates of the extremum in the data. These
# can be continuous, since the exact extremum location is estimated using an
# interpolation procedure. Thus xex and yex should not be assumed to fall on
# the x,y grid provided as input. "type" indicates the type of extremum: max,
# min, saddle, or unknown. The latter case, unknown, indicates that the
# second derivative test has failed, a case in which the identity of the
# extremum is not determined by this algorithm.
def find_extrema(x,y,data,frac_tol=0.0,presmooths=10):
    X, Y = _check_input_arrays(x, y)
    dX = X[:,1:] - X[:,:-1]
    dY = Y[1:,:] - Y[:-1,:]

    for i in range(presmooths):
        data = smoothField2d(data)

    # dx_data = (data[:,1:] - data[:,:-1]) / dX
    # dy_data = (data[1:,:] - data[:-1,:]) / dY
    # JM: 2022-08-03: I have found that using Sobel filters makes the code
    # easier to read, generally gives fewer spurrious extrema, and tends
    # to give more trustworthy extrema classifications (e.g., mins are truly
    # mins, saddles are truly saddles, etc.). This may have to do with the
    # fact that the Sobel filters include an inherent smoothing step, which
    # prevents the derivatives from getting too noisy.
    # See: https://en.wikipedia.org/wiki/Image_derivative
    dx_data = sobelx(x,y,data)
    dy_data = sobely(x,y,data)
    
    dxy_data = sobely(x,y,dx_data)
    dxx_data = sobelx(x,y,dx_data)
    dyy_data = sobely(x,y,dy_data)
    
    # extremum bit flags
    # We flag extrema, locations where both components of the gradient
    # change sign, as 1's in this array
    # Include potentially exact extrema locations, where one of the
    # derivatives is exactly zero on a gridpoint.
    exbf = np.zeros(data.shape,dtype=int)
    # exbf[1:-1,1:-1] = ( (dx_data[1:-1,1:]*dx_data[1:-1,:-1] <= 0) &
    #                     (dy_data[1:,1:-1]*dy_data[:-1,1:-1] <= 0) )
    # The gradient goes to zero in this cell if both of its components
    # switch sign.
    # A component switches sign if its sign on one corner of the cell
    # does not match it's sign on all other corners of the cell.
    exbf[1:-2,1:-2] = ( 
                        ((dx_data[1:-2,1:-2]*dx_data[1:-2,2:-1] <= 0) |
                         (dx_data[1:-2,1:-2]*dx_data[2:-1,1:-2] <= 0) |
                         (dx_data[1:-2,1:-2]*dx_data[2:-1,2:-1] <= 0)) &
                        ((dy_data[1:-2,1:-2]*dy_data[1:-2,2:-1] <= 0) |
                         (dy_data[1:-2,1:-2]*dy_data[2:-1,1:-2] <= 0) |
                         (dy_data[1:-2,1:-2]*dy_data[2:-1,2:-1] <= 0))
                      )

    exidx = np.where(exbf)
    
    output = []
    for yidx, xidx in zip(*exidx):
        # # Compute half gridpoints adjacent to each extremum location
        # xmid1 = 0.5*(X[yidx,xidx-1]+X[yidx,xidx])
        # xmid2 = 0.5*(X[yidx,xidx]+X[yidx,xidx+1])
        # ymid1 = 0.5*(Y[yidx-1,xidx]+Y[yidx,xidx])
        # ymid2 = 0.5*(Y[yidx,xidx]+Y[yidx+1,xidx])

        # # Compute averaging weights from half-gridpoints to extremum point.
        # wx1 = (xmid2 - X[yidx,xidx]) / (xmid2 - xmid1)
        # wx2 = (X[yidx,xidx] - xmid1) / (xmid2 - xmid1)
        # wy1 = (ymid2 - Y[yidx,xidx]) / (ymid2 - ymid1)
        # wy2 = (Y[yidx,xidx] - ymid1) / (ymid2 - ymid1)

        # # The x/y derivatives of the data, averaged to the extremum point.
        # dx_data_ex = \
        #     wx1*(data[yidx,xidx]-data[yidx,xidx-1])/(X[yidx,xidx]-X[yidx,xidx-1]) + \
        #     wx2*(data[yidx,xidx+1]-data[yidx,xidx])/(X[yidx,xidx+1]-X[yidx,xidx]) 
        # dy_data_ex = \
        #     wy1*(data[yidx,xidx]-data[yidx-1,xidx])/(Y[yidx,xidx]-Y[yidx-1,xidx]) + \
        #     wy2*(data[yidx+1,xidx]-data[yidx,xidx])/(Y[yidx+1,xidx]-Y[yidx,xidx]) 

        # # The xx and yy double derivatives of the data. No averaging is
        # # required after the derivatives are evaluated: they are already known
        # # at the extremum location
        # dxx_data_ex = \
        #     ((data[yidx,xidx]-data[yidx,xidx-1])/(X[yidx,xidx]-X[yidx,xidx-1])
        #     -(data[yidx,xidx+1]-data[yidx,xidx])/(X[yidx,xidx+1]-X[yidx,xidx]))/\
        #      (xmid2 - xmid1)
        # dyy_data_ex = \
        #     ((data[yidx,xidx]-data[yidx-1,xidx])/(Y[yidx,xidx]-Y[yidx-1,xidx])
        #     -(data[yidx+1,xidx]-data[yidx,xidx])/(Y[yidx+1,xidx]-Y[yidx,xidx]))/\
        #      (ymid2 - ymid1)

        # # The four xy double derivatives that then need to be averaged to the
        # # extremum location
        # dxy_data_x1y1 = \
        #     ((data[yidx  ,xidx]-data[yidx  ,xidx-1])/(X[yidx  ,xidx]-X[yidx  ,xidx-1])
        #     -(data[yidx-1,xidx]-data[yidx-1,xidx-1])/(X[yidx-1,xidx]-X[yidx-1,xidx-1]))/\
        #      (Y[yidx,xidx]-Y[yidx-1,xidx])
        # dxy_data_x2y1 = \
        #     ((data[yidx  ,xidx+1]-data[yidx  ,xidx])/(X[yidx  ,xidx+1]-X[yidx  ,xidx])
        #     -(data[yidx-1,xidx+1]-data[yidx-1,xidx])/(X[yidx-1,xidx+1]-X[yidx-1,xidx]))/\
        #      (Y[yidx,xidx]-Y[yidx-1,xidx])
        # dxy_data_x1y2 = \
        #     ((data[yidx+1,xidx]-data[yidx+1,xidx-1])/(X[yidx+1,xidx]-X[yidx+1,xidx-1])
        #     -(data[yidx  ,xidx]-data[yidx  ,xidx-1])/(X[yidx  ,xidx]-X[yidx  ,xidx-1]))/\
        #      (Y[yidx+1,xidx]-Y[yidx,xidx])
        # dxy_data_x2y2 = \
        #     ((data[yidx+1,xidx+1]-data[yidx+1,xidx])/(X[yidx+1,xidx+1]-X[yidx+1,xidx])
        #     -(data[yidx  ,xidx+1]-data[yidx  ,xidx])/(X[yidx  ,xidx+1]-X[yidx  ,xidx]))/\
        #      (Y[yidx+1,xidx]-Y[yidx,xidx])

        # # The xy double derivatives averaged to the extremum location
        # dxy_data_ex = wx1*wy1*dxy_data_x1y1 + \
        #               wx2*wy1*dxy_data_x2y1 + \
        #               wx1*wy2*dxy_data_x1y2 + \
        #               wx2*wy2*dxy_data_x2y2

        dx_data_ex = dx_data[yidx,xidx]
        dy_data_ex = dy_data[yidx,xidx]
        dxy_data_ex = dxy_data[yidx,xidx]
        dxx_data_ex = dxx_data[yidx,xidx]
        dyy_data_ex = dyy_data[yidx,xidx]
        xmid2 = X[yidx,xidx+1]
        xmid1 = X[yidx,xidx-1]
        ymid2 = Y[yidx+1,xidx]
        ymid1 = Y[yidx-1,xidx]

        # The determinant of the Hessian matrix, which gives the identity of
        # this extremum.
        det = dxx_data_ex*dyy_data_ex - dxy_data_ex**2
        # Numerical noise will (probably?) always prevent the determinant from
        # being zero, even when the extremum in question should probably
        # fail the second derivative test.
        zero = frac_tol * (data[yidx,xidx] / ((ymid2-ymid1)*(xmid2-xmid1)))**2 
        extype = "unknown"
        if det > zero:
            if dxx_data_ex > 0:
                extype = "min"
            else:
                extype = "max"
        elif det < -zero:
            extype = "saddle"
        # print("determinant,zero,extype=%g,%g,%s" % (det,zero,extype))

        xex = X[yidx,xidx]
        yex = Y[yidx,xidx]
        # We can only refine the estimate of the extremum location if the
        # second derivative test is conclusive.
        if not (extype == "unknown"):
            xex = 0.5*(dy_data_ex*dxy_data_ex - dx_data_ex*dyy_data_ex)/det + X[yidx,xidx]
            yex = 0.5*(dx_data_ex*dxy_data_ex - dy_data_ex*dxx_data_ex)/det + Y[yidx,xidx]

        output.append((xex,yex,extype))

    return output


# Calculate the theta component of the vector potential (times r)
# Br,Bz are 2d arrays with values known at the corresponding gridpoints of r,z
# r,z can be...
#    - 1d arrays of cylindrical radius and z, respectively. In this case, it
#      is assumed that B[i,j] is known at the coordinate (r[i],z[j]).
#    - 2d arrays of cylindrical radius and z, respectively. In this case, it
#      is assumed that B[i,j] is known at the coorindate (r[i,j],z[i,j]).
#
# r*Ath is known on the same gridpoints r,z as supplied for Br/Bz.
# As a boundary condition, r*Ath[0,0] is set to zero.
def calc_rAth(r,z,Br,Bz):
    R, Z = _check_input_arrays(r, z)

    # Average to half-gridpoints for integration (so that Ath is known at the
    # same gridpoints as B).
    # First we need to know Bz at half-integer r-points to integrate in r.
    dR = R[:,1:] - R[:,:-1]
    R_halfr = 0.5 * (R[:,1:] + R[:,:-1])
    dR2_up = R[:,1:]**2 - R_halfr**2
    dR2_do = R_halfr**2 - R[:,:-1]**2
    dR2 = R[:,1:]**2 - R[:,:-1]**2
    vol_inv = 1. / (dR2)
    w1 = vol_inv * dR2_up
    w2 = vol_inv * dR2_do
    Bz_halfr = w1 * Bz[:,:-1] + w2 * Bz[:,1:]
    # w11 = vol_inv * dZ * dR2_up
    # w12 = vol_inv * dZ * dR2_do
    # w21 = vol_inv * dZ * dR2_up
    # w22 = vol_inv * dZ * dR2_do
    # Bz_mid = w11*Bz[:-1,:-1] + w12*Bz[:-1,1:] + w21*Bz[1:,:-1] + w22*Bz[1:,1:]
    # Br_mid = w11*Br[:-1,:-1] + w12*Br[:-1,1:] + w21*Br[1:,:-1] + w22*Br[1:,1:]

    # Next, we need to know Br at half-integer z-points to integrate in z.
    dZ = Z[1:,:] - Z[:-1,:]
    w1 = 0.5
    w2 = 0.5
    Br_halfz = w1*Br[:-1,:] + w2*Br[1:,:]
    R_halfz = 0.5 * (R[1:,:] + R[:-1,:])  # == R[1:,:]

    # rAth is set to zero at these points
    # r0 = R[0,0]
    # z0 = Z[0,0]

    # 2pi*(B flux through a disk from r0 to R at a height of Z)
    flux_r0toR_Z = np.cumsum(dR * R_halfr * Bz_halfr, axis = 1)
    # cumsum doesn't insert a leading zero, which we need as a boundary
    # condition (there's no flux from r0 to r0...)
    zeros = np.zeros( (flux_r0toR_Z.shape[0], 1) )
    flux_r0toR_Z = np.concatenate( (zeros, flux_r0toR_Z), axis = 1)
    # 2pi*(B flux through a cylindrical wall from z0 to Z at a radius R)
    flux_R_z0toZ = np.cumsum(R_halfz * Br_halfz * dZ, axis = 0)
    # Same insertion of 0 as above (there's no flux from z0 to z0)
    zeros = np.zeros( (1, flux_R_z0toZ.shape[1]) )
    flux_R_z0toZ = np.concatenate( (zeros, flux_R_z0toZ), axis = 0)

    # We only actually need one slice of this array for the boundary condition
    flux_R_z0toZ = flux_R_z0toZ[:,0][:,np.newaxis]

    rAth = flux_r0toR_Z - flux_R_z0toZ

    return rAth


# Interpolate the 2d field "field" known at the gridpoints r,z to the value at
# rinterp, zinterp
#
# r,z can be...
#    - 1d arrays of cylindrical radius and z, respectively. In this case, it
#      is assumed that field[i,j] is known at the coordinate (r[i],z[j]).
#    - 2d arrays of cylindrical radius and z, respectively. In this case, it
#      is assumed that field[i,j] is known at the coorindate (r[i,j],z[i,j]).
def interp_field_2d(rinterp, zinterp, r, z, field):
    R, Z = _check_input_arrays(r, z)
    r, z = R[0,:], Z[:,0]

    ridx = np.argmin(np.abs(r - rinterp))
    zidx = np.argmin(np.abs(z - zinterp))
    if r[ridx] > rinterp:
        ridx -= 1
    if z[zidx] > zinterp:
        zidx -= 1

    ru = r[ridx+1]
    rl = r[ridx]
    zu = z[zidx+1]
    zl = z[zidx]

    atot = (zu - zl) * (ru**2 - rl**2)
    aruzu = (zu - zinterp) * (ru**2 - rinterp**2) / atot
    aruzl = (zinterp - zl) * (ru**2 - rinterp**2) / atot
    arlzu = (zu - zinterp) * (rinterp**2 - rl**2) / atot
    arlzl = (zinterp - zl) * (rinterp**2 - rl**2) / atot
    # atot = (zu - zl) * (ru - rl)
    # aruzu = (zu - zinterp) * (ru - rinterp) / atot
    # aruzl = (zinterp - zl) * (ru - rinterp) / atot
    # arlzu = (zu - zinterp) * (rinterp - rl) / atot
    # arlzl = (zinterp - zl) * (rinterp - rl) / atot

    out = aruzu * field[zidx  ,ridx  ] + \
          aruzl * field[zidx+1,ridx  ] + \
          arlzu * field[zidx  ,ridx+1] + \
          arlzl * field[zidx+1,ridx+1]

    return out

# Create an evenly spaced array from a to b with n entries.
# If c == None, this defaults to np.linspace(a, b, num = nlevels)
# Optionally specify some c between a and b that is demanded to be in the
# output array. In that case, the array will have nlevels entries, contain a,
# contain c, and approximately terminate at b.
def choose_nice_contour_levels(a, b, nlevels, c = None):
    # Do something that doesn't give NaNs when a == b
    if np.abs(b - a) <= 1e-14 * np.abs(a):
        return np.array([a])

    if c is None:
        return np.linspace(a, b, num = nlevels)

    # c is not None
    nsteps = nlevels - 1
    step = (c - a) / np.round( (c - a) / (b - a) * nsteps)

    return step * np.arange(nlevels) + a
