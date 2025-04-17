import numpy as np
import h5py
from vtkmodules.numpy_interface import algorithms as algs
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtk import *

import os

import plotutil as pu

c = 299792458e2
e = 4.8032068e-10
me=9.1093897e-28


def avg_neighbors(arr3d):
    return 1./8. * (
        arr3d[1:  ,1:  ,1:  ] + arr3d[1:  ,1:  , :-1] +
        arr3d[1:  , :-1,1:  ] + arr3d[1:  , :-1, :-1] +
        arr3d[ :-1,1:  ,1:  ] + arr3d[ :-1,1:  , :-1] +
        arr3d[ :-1, :-1,1:  ] + arr3d[ :-1, :-1, :-1]
    )

# Cut the resolution of a 3d array in half in each dimension
# The shape of arr3d in each dimension needs to be even
# A smoothing pass is done so that each cell in the new array
# contains the average of the corresponding 8 neighboring cells in the
# original array.
def half_res(arr3d):
    output = np.empty(np.array(arr3d.shape)//2)
    smoothed = avg_neighbors(arr3d)
    output[:,:,:] = smoothed[::2,::2,::2]
    return output

def h5tovtk_sgrid(*args,**kwargs):
    print("Converting to vtk Structured Grid")

    simdir     = "./../"
    datadir    = os.path.join(simdir, "data")
    fielddir   = os.path.join(datadir, "fields")
    densitydir = os.path.join(datadir, "densities")
    params     = pu.get_sim_params_dict(simdir)

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
    # 'Nx' in this function is what Zeltron calls 'NCX'
    # I.e., the number of cells (not grid points) in x
    Nx   = int(params["NX"] - 1)
    Ny   = int(params["NY"] - 1)
    Nz   = int(params["NZ"] - 1)

    B0      = params["B0 [G]"]
    rlc     = params["R_LC [cm]"]
    rnozzle = params["R_nozzle [cm]"]
    zmono   = params["z_monopole [cm]"]

    # Derived parameters
    omega = c / rlc
    ratio = rnozzle / zmono
    ngj = B0 * omega / (2 * np.pi * c * e)
    Lbz = c * B0**2 * zmono**4 / (3. * rlc**2)
    # Lbz in the above line is the Poynting power for a semisphere
    # However, the Poynting flux is only injected through a cone with opening
    # angle given by tan(theta) = rnozzle / zmono = ratio
    # We need to calculate the fraction ('correction' below) of the hemispherical
    # injected Poynting power that passes through this cone.
    correction = (1. - 1. / np.sqrt(1. + ratio**2) * (3. * ratio**2 / 2. + 1.) / (1. + ratio**2))
    Lbz = Lbz * correction
    # Sznorm = B0**2 * zmono**4 * c / (4 * np.pi)
    Sznorm = B0**2 * c / (4 * np.pi)

    #===========================================================================
    # Parse args
    it = args[0]

    rrp2 = 0
    if "redResPow2" in kwargs.keys():
        rrp2 = kwargs["redResPow2"]

    #===========================================================================
    # Grid reconstruction

    # Nph, Nth, Nr = np.array(
    #     h5py.File(dir + f"/densities/densities_{it}.h5")["electrons"]
    # ).shape
    # h5file = h5py.File(datadir + f"/densities/mapxyz_electrons_bg_{it}.h5")
    # Nz, Ny, Nx = np.array(h5file["field"]).shape
    # h5file.close()

    # input_params = np.loadtxt(dir + "input_params.dat", skiprows=1)
    # phys_params = np.loadtxt(dir + "phys_params.dat", skiprows=1)
    # B0 = phys_params[0]
    # n0 = phys_params[2]
    # a = phys_params[3]
    # rh = 1.0 + np.sqrt(1.0 - a * a)
    # rmin = input_params[12]
    # rmax = input_params[13]
    # thmin = input_params[15]
    # thmax = input_params[16]
    # r = np.geomspace(rmin, rmax, Nr)
    # th = np.linspace(thmin, thmax, Nth)
    # ph = np.linspace(0, 2 * np.pi, Nph)
    x = np.linspace(xmin, xmax, Nx + 1)
    y = np.linspace(ymin, ymax, Ny + 1)
    z = np.linspace(zmin, zmax, Nz + 1)
    # Densities have one fewer cell in each dimension with respect to fields
    # and the grid. Fix this by shaving off the last cell in the quantities
    # with extra data.
    x = x[:-1]
    y = y[:-1]
    z = z[:-1]
    # R, Th, Ph = np.meshgrid(r, th, ph, indexing="ij")
    # costh, sinth = np.cos(Th), np.sin(Th)
    # cosph, sinph = np.cos(Ph), np.sin(Ph)

    # filename = dir + f"fields/fields_{it}.h5"
    # h5file = h5py.File(filename)

    # dens_elec = np.array(h5py.File(dir + f"/densities/densities_{it}.h5")["electrons"])
    # dens_posi = np.array(h5py.File(dir + f"/densities/densities_{it}.h5")["positrons"])

    # dens_elec = dens_elec.swapaxes(0, 2)
    # dens_posi = dens_posi.swapaxes(0, 2)

    #===========================================================================
    # Read simulation data

    print("Loading densities...")
    nbge = pu.readArrayFromHdf5(os.path.join(densitydir, "mapxyz_electrons_bg_%s.h5" % it), "field")
    nbgi = pu.readArrayFromHdf5(os.path.join(densitydir, "mapxyz_ions_bg_%s.h5" % it), "field")
    nje  = pu.readArrayFromHdf5(os.path.join(densitydir, "mapxyz_electrons_drift_%s.h5" % it), "field")
    nji  = pu.readArrayFromHdf5(os.path.join(densitydir, "mapxyz_ions_drift_%s.h5" % it), "field")
    print("Loading B field...")
    Bz   = pu.readArrayFromHdf5(os.path.join(fielddir, "Bz_%s.h5" % it), "field")
    By   = pu.readArrayFromHdf5(os.path.join(fielddir, "By_%s.h5" % it), "field")
    Bx   = pu.readArrayFromHdf5(os.path.join(fielddir, "Bx_%s.h5" % it), "field")
    print("Loading E field...")
    Ez   = pu.readArrayFromHdf5(os.path.join(fielddir, "Ez_%s.h5" % it), "field")
    Ey   = pu.readArrayFromHdf5(os.path.join(fielddir, "Ey_%s.h5" % it), "field")
    Ex   = pu.readArrayFromHdf5(os.path.join(fielddir, "Ex_%s.h5" % it), "field")

    print("Data loaded.")

    # zB = pu.readArrayFromHdf5(os.path.join(fielddir, "By_%s.h5" % it), "axis0coords")
    # yB = pu.readArrayFromHdf5(os.path.join(fielddir, "By_%s.h5" % it), "axis1coords")
    # xB = pu.readArrayFromHdf5(os.path.join(fielddir, "Bx_%s.h5" % it), "axis2coords")
    # zE = pu.readArrayFromHdf5(os.path.join(fielddir, "Ey_%s.h5" % it), "axis0coords")
    # yE = pu.readArrayFromHdf5(os.path.join(fielddir, "Ey_%s.h5" % it), "axis1coords")
    # xE = pu.readArrayFromHdf5(os.path.join(fielddir, "Ex_%s.h5" % it), "axis2coords")
    # zn = pu.readArrayFromHdf5(os.path.join(densitydir, "mapxyz_electrons_bg_%s.h5" % it), "axis0coords")
    # yn = pu.readArrayFromHdf5(os.path.join(densitydir, "mapxyz_electrons_bg_%s.h5" % it), "axis1coords")
    # xn = pu.readArrayFromHdf5(os.path.join(densitydir, "mapxyz_electrons_bg_%s.h5" % it), "axis2coords")

    nbg = nbge + nbgi
    nj  = nje + nji
    ntot = nbg + nj

    ntot = ntot.swapaxes(0, 2)
    nbg  = nbg.swapaxes(0, 2)
    nj   = nj.swapaxes(0, 2)
    Bx   = Bx.swapaxes(0, 2)
    By   = By.swapaxes(0, 2)
    Bz   = Bz.swapaxes(0, 2)
    Ex   = Ex.swapaxes(0, 2)
    Ey   = Ey.swapaxes(0, 2)
    Ez   = Ez.swapaxes(0, 2)

    # Densities have one fewer cell in each dimension with respect to fields
    # and the grid. Fix this by shaving off the last cell in the quantities
    # with extra data.
    Bx = Bx[:-1,:-1,:-1]
    By = By[:-1,:-1,:-1]
    Bz = Bz[:-1,:-1,:-1]
    Ex = Ex[:-1,:-1,:-1]
    Ey = Ey[:-1,:-1,:-1]
    Ez = Ez[:-1,:-1,:-1]

    # import pdb; pdb.set_trace()

    #===========================================================================
    # Reduce resolution
    # To avoid averaging noisy quantities, it is probably best to reduce the
    # resolution before calculating derived quantities below (rather than
    # first calculating derived quantities and then reducing the grid)
    if rrp2 > 0:
        print("Cutting resolution by 2^%d" % rrp2)

    for reduction_step in range(rrp2):
        Bx   = half_res(Bx)
        By   = half_res(By)
        Bz   = half_res(Bz)
        Ex   = half_res(Ex)
        Ey   = half_res(Ey)
        Ez   = half_res(Ez)
        ntot = half_res(ntot)
        nbg  = half_res(nbg)
        nj   = half_res(nj)

        x = x[::2]
        y = y[::2]
        z = z[::2]
        dx = 2 * dx
        dy = 2 * dy
        dz = 2 * dz

    X, Y, Z = np.meshgrid(x, y, z, indexing = "ij")

    #===========================================================================
    # Derived fields

    print("Calculating derived fields...")
    Sz = c * (Ex * By - Ey * Bx) / (4 * np.pi)
    Bsq = Bx * Bx + By * By + Bz * Bz
    sigma = np.zeros(Bsq.shape)
    ixnz = np.where(ntot > 0.0)
    sigma[ixnz] = Bsq[ixnz] / (4 * np.pi * me * c * c * ntot[ixnz])

    # # Jru = np.array(h5py.File(dir+f'currents/currents_{it}.h5','r')['Jr']).swapaxes(0,2)
    # # Jthu = np.array(h5py.File(dir+f'currents/currents_{it}.h5','r')['Jth']).swapaxes(0,2)
    # # Jphu = np.array(h5py.File(dir+f'currents/currents_{it}.h5','r')['Jph']).swapaxes(0,2)

    # Bru = np.array(h5file["Bru"]).swapaxes(0, 2)
    # Bthu = np.array(h5file["Bthu"]).swapaxes(0, 2)
    # Bphu = np.array(h5file["Bphu"]).swapaxes(0, 2)

    # # Erd = np.array(h5file["Erd"]).swapaxes(0, 2)
    # Ethd = np.array(h5file["Ethd"]).swapaxes(0, 2)
    # Ephd = np.array(h5file["Ephd"]).swapaxes(0, 2)

    # # Hrd = np.array(h5file["Hrd"]).swapaxes(0, 2)
    # Hthd = np.array(h5file["Hthd"]).swapaxes(0, 2)
    # Hphd = np.array(h5file["Hphd"]).swapaxes(0, 2)

    # Sru = (Ethd * Hphd - Ephd * Hthd) / (4 * np.pi * sgam)

    # Bsq = gr * Bru * Bru + gth * Bthu * Bthu + gph * Bphu * Bphu + 2 * grph * Bru * Bphu
    # sigma = Bsq / (4 * np.pi * (dens_elec + dens_posi))

    # LBZ = B0**2 * a**2 / 96

    # # Setting field in a Cartesian grid

    # Bxu = sinth * cosph * Bru + R * costh * cosph * Bthu - R * sinth * sinph * Bphu
    # Byu = sinth * sinph * Bru + R * costh * sinph * Bthu + R * sinth * cosph * Bphu
    # Bzu = costh * Bru - R * sinth * Bthu

    # Jxu = sinth*cosph*Jru + R*costh*cosph*Jthu - R*sinth*sinph*Jphu
    # Jyu = sinth*sinph*Jru + R*costh*sinph*Jthu + R*sinth*cosph*Jphu
    # Jzu = costh*Jru - R*sinth*Jthu

    #===========================================================================
    # Initialize vtk arrays + point grid
    print("Creating the grid...")

    dims = Bx.shape
    # dims = Bru.shape
    sgrid = vtkStructuredGrid()
    # sgrid = vtkImageData()
    sgrid.SetDimensions(dims)    
    # sgrid.SetOrigin(xmin,ymin,zmin)
    # sgrid.SetSpacing(dx,dy,dz)

    data_points = algs.make_vector(
        X.flatten('F'), Y.flatten('F'), Z.flatten('F')
    )
    points = vtkPoints()
    points.Allocate(dims[0] * dims[1] * dims[2])
    points.SetData(dsa.numpyTovtkDataArray(data_points,'Points'))

    np_Bfield = algs.make_vector(
        Bx.flatten('F') / B0, By.flatten('F') / B0, Bz.flatten('F') / B0
    )
    Bfield = dsa.numpyTovtkDataArray(np_Bfield, "B field")

    density = dsa.numpyTovtkDataArray(ntot.flatten('F') / ngj, "Density")

    magnetization = dsa.numpyTovtkDataArray(sigma.flatten('F'), "Cold sigma")

    poynting = dsa.numpyTovtkDataArray(Sz.flatten('F') / Sznorm, "Poynting z")

    sgrid.SetPoints(points)
    sgrid.GetPointData().SetVectors(Bfield)
    sgrid.GetPointData().AddArray(density)
    sgrid.GetPointData().AddArray(poynting)
    sgrid.GetPointData().AddArray(magnetization)

    print("Grid created")

    #===========================================================================
    # Save to vtk file
    print("Writing file...")

    # Writing data
    if not os.path.isdir("./vtk"):
        os.mkdir("./vtk")
    writer = vtkStructuredGridWriter()
    # writer = vtkImageWriter()
    writer.SetFileTypeToBinary()
    it = int(it)
    writer.SetFileName(f"vtk/data_{it:06d}.vtk")
    writer.SetInputData(sgrid)
    writer.Write()
    # h5file.close()

def h5tovtk_rgrid(*args,**kwargs):
    print("Converting to vtk Rectilinear Grid")

    simdir     = "./../"
    datadir    = os.path.join(simdir, "data")
    fielddir   = os.path.join(datadir, "fields")
    densitydir = os.path.join(datadir, "densities")
    params     = pu.get_sim_params_dict(simdir)

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
    # 'Nx' in this function is what Zeltron calls 'NCX'
    # I.e., the number of cells (not grid points) in x
    Nx   = int(params["NX"] - 1)
    Ny   = int(params["NY"] - 1)
    Nz   = int(params["NZ"] - 1)

    B0      = params["B0 [G]"]
    rlc     = params["R_LC [cm]"]
    rnozzle = params["R_nozzle [cm]"]
    zmono   = params["z_monopole [cm]"]

    # Derived parameters
    omega = c / rlc
    ratio = rnozzle / zmono
    ngj = B0 * omega / (2 * np.pi * c * e)
    # Sznorm = B0**2 * zmono**4 * c / (4 * np.pi)
    Sznorm = B0**2 * c / (4 * np.pi)

    #===========================================================================
    # Parse args
    it = args[0]

    rrp2 = 0
    if "redResPow2" in kwargs.keys():
        rrp2 = kwargs["redResPow2"]

    #===========================================================================
    # Grid reconstruction

    x = np.linspace(xmin, xmax, Nx + 1)
    y = np.linspace(ymin, ymax, Ny + 1)
    z = np.linspace(zmin, zmax, Nz + 1)
    # Densities have one fewer cell in each dimension with respect to fields
    # and the grid. Fix this by shaving off the last cell in the quantities
    # with extra data.
    # x = x[:-1]
    # y = y[:-1]
    # z = z[:-1]

    #===========================================================================
    # Read simulation data

    print("Loading densities...")
    nbge = pu.readArrayFromHdf5(os.path.join(densitydir, "mapxyz_electrons_bg_%s.h5" % it), "field")
    nbgi = pu.readArrayFromHdf5(os.path.join(densitydir, "mapxyz_ions_bg_%s.h5" % it), "field")
    nje  = pu.readArrayFromHdf5(os.path.join(densitydir, "mapxyz_electrons_drift_%s.h5" % it), "field")
    nji  = pu.readArrayFromHdf5(os.path.join(densitydir, "mapxyz_ions_drift_%s.h5" % it), "field")
    print("Loading B field...")
    Bz   = pu.readArrayFromHdf5(os.path.join(fielddir, "Bz_%s.h5" % it), "field")
    By   = pu.readArrayFromHdf5(os.path.join(fielddir, "By_%s.h5" % it), "field")
    Bx   = pu.readArrayFromHdf5(os.path.join(fielddir, "Bx_%s.h5" % it), "field")
    print("Loading E field...")
    Ez   = pu.readArrayFromHdf5(os.path.join(fielddir, "Ez_%s.h5" % it), "field")
    Ey   = pu.readArrayFromHdf5(os.path.join(fielddir, "Ey_%s.h5" % it), "field")
    Ex   = pu.readArrayFromHdf5(os.path.join(fielddir, "Ex_%s.h5" % it), "field")

    print("Data loaded.")

    nbg = nbge + nbgi
    nj  = nje + nji
    ntot = nbg + nj

    ntot = ntot.swapaxes(0, 2)
    nbg  = nbg.swapaxes(0, 2)
    nj   = nj.swapaxes(0, 2)
    Bx   = Bx.swapaxes(0, 2)
    By   = By.swapaxes(0, 2)
    Bz   = Bz.swapaxes(0, 2)
    Ex   = Ex.swapaxes(0, 2)
    Ey   = Ey.swapaxes(0, 2)
    Ez   = Ez.swapaxes(0, 2)

    # Densities have one fewer cell in each dimension with respect to fields
    # and the grid. Fix this by shaving off the last cell in the quantities
    # with extra data.
    Bx = Bx[:-1,:-1,:-1]
    By = By[:-1,:-1,:-1]
    Bz = Bz[:-1,:-1,:-1]
    Ex = Ex[:-1,:-1,:-1]
    Ey = Ey[:-1,:-1,:-1]
    Ez = Ez[:-1,:-1,:-1]

    #===========================================================================
    # Reduce resolution
    # To avoid averaging noisy quantities, it is probably best to reduce the
    # resolution before calculating derived quantities below (rather than
    # first calculating derived quantities and then reducing the grid)
    if rrp2 > 0:
        print("Cutting resolution by 2^%d" % rrp2)

    for reduction_step in range(rrp2):
        Bx   = half_res(Bx)
        By   = half_res(By)
        Bz   = half_res(Bz)
        Ex   = half_res(Ex)
        Ey   = half_res(Ey)
        Ez   = half_res(Ez)
        ntot = half_res(ntot)
        nbg  = half_res(nbg)
        nj   = half_res(nj)

        x = x[::2]
        y = y[::2]
        z = z[::2]
        dx = 2 * dx
        dy = 2 * dy
        dz = 2 * dz

    X, Y, Z = np.meshgrid(x, y, z, indexing = "ij")

    #===========================================================================
    # Derived fields

    print("Calculating derived fields...")
    Sz = c * (Ex * By - Ey * Bx) / (4 * np.pi)
    Bsq = Bx * Bx + By * By + Bz * Bz
    # sigma = np.zeros(Bsq.shape)
    # ixnz = np.where(ntot > 0.0)
    # sigma[ixnz] = Bsq[ixnz] / (4 * np.pi * me * c * c * ntot[ixnz])
    sigma = Bsq / (4 * np.pi * me * c * c * ntot)

    #===========================================================================
    # Initialize vtk arrays + point grid
    print("Creating the grid...")

    dims = Bx.shape
    rgrid = vtkRectilinearGrid()
    rgrid.SetDimensions(dims)
    # xlin = dsa.numpyTovtkDataArray(X.flatten('F'))
    # ylin = dsa.numpyTovtkDataArray(Y.flatten('F'))
    # zlin = dsa.numpyTovtkDataArray(Z.flatten('F'))
    xlin = dsa.numpyTovtkDataArray(x)
    ylin = dsa.numpyTovtkDataArray(y)
    zlin = dsa.numpyTovtkDataArray(z)
    rgrid.SetXCoordinates(xlin)
    rgrid.SetYCoordinates(ylin)
    rgrid.SetZCoordinates(zlin)

    np_Bfield = algs.make_vector(
        Bx.flatten('F') / B0, By.flatten('F') / B0, Bz.flatten('F') / B0
    )
    Bfield = dsa.numpyTovtkDataArray(np_Bfield, "B field")

    density = dsa.numpyTovtkDataArray(ntot.flatten('F') / ngj, "Density")

    magnetization = dsa.numpyTovtkDataArray(sigma.flatten('F'), "Cold sigma")

    poynting = dsa.numpyTovtkDataArray(Sz.flatten('F') / Sznorm, "Poynting z")

    rgrid.GetPointData().SetVectors(Bfield)
    rgrid.GetPointData().AddArray(density)
    rgrid.GetPointData().AddArray(poynting)
    rgrid.GetPointData().AddArray(magnetization)

    print("Grid created")

    #===========================================================================
    # Save to vtk file
    print("Writing file...")

    if not os.path.isdir("./vtk"):
        os.mkdir("./vtk")

    it = int(it)
    fname = f"vtk/data_{it:06d}.vtk"
    if "save" in kwargs.keys():
        fname = kwargs["save"]

    # Writing data
    writer = vtkRectilinearGridWriter()
    writer.SetFileTypeToBinary()
    writer.SetFileName(fname)
    writer.SetInputData(rgrid)
    writer.Write()


def h5tovtk_imagedata(*args,**kwargs):
    print("Converting to vtk Image Data.")
    print("This is currently broken.")

    simdir     = "./../"
    datadir    = os.path.join(simdir, "data")
    fielddir   = os.path.join(datadir, "fields")
    densitydir = os.path.join(datadir, "densities")
    params     = pu.get_sim_params_dict(simdir)

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
    # 'Nx' in this function is what Zeltron calls 'NCX'
    # I.e., the number of cells (not grid points) in x
    Nx   = int(params["NX"] - 1)
    Ny   = int(params["NY"] - 1)
    Nz   = int(params["NZ"] - 1)

    B0      = params["B0 [G]"]
    rlc     = params["R_LC [cm]"]
    rnozzle = params["R_nozzle [cm]"]
    zmono   = params["z_monopole [cm]"]

    # Derived parameters
    omega = c / rlc
    ratio = rnozzle / zmono
    ngj = B0 * omega / (2 * np.pi * c * e)
    # Sznorm = B0**2 * zmono**4 * c / (4 * np.pi)
    Sznorm = B0**2 * c / (4 * np.pi)

    #===========================================================================
    # Parse args
    it = args[0]

    rrp2 = 0
    if "redResPow2" in kwargs.keys():
        rrp2 = kwargs["redResPow2"]

    #===========================================================================
    # Grid reconstruction

    x = np.linspace(xmin, xmax, Nx + 1)
    y = np.linspace(ymin, ymax, Ny + 1)
    z = np.linspace(zmin, zmax, Nz + 1)
    # Densities have one fewer cell in each dimension with respect to fields
    # and the grid. Fix this by shaving off the last cell in the quantities
    # with extra data.
    x = x[:-1]
    y = y[:-1]
    z = z[:-1]

    #===========================================================================
    # Read simulation data

    print("Loading densities...")
    nbge = pu.readArrayFromHdf5(os.path.join(densitydir, "mapxyz_electrons_bg_%s.h5" % it), "field")
    nbgi = pu.readArrayFromHdf5(os.path.join(densitydir, "mapxyz_ions_bg_%s.h5" % it), "field")
    nje  = pu.readArrayFromHdf5(os.path.join(densitydir, "mapxyz_electrons_drift_%s.h5" % it), "field")
    nji  = pu.readArrayFromHdf5(os.path.join(densitydir, "mapxyz_ions_drift_%s.h5" % it), "field")
    print("Loading B field...")
    Bz   = pu.readArrayFromHdf5(os.path.join(fielddir, "Bz_%s.h5" % it), "field")
    By   = pu.readArrayFromHdf5(os.path.join(fielddir, "By_%s.h5" % it), "field")
    Bx   = pu.readArrayFromHdf5(os.path.join(fielddir, "Bx_%s.h5" % it), "field")
    print("Loading E field...")
    Ez   = pu.readArrayFromHdf5(os.path.join(fielddir, "Ez_%s.h5" % it), "field")
    Ey   = pu.readArrayFromHdf5(os.path.join(fielddir, "Ey_%s.h5" % it), "field")
    Ex   = pu.readArrayFromHdf5(os.path.join(fielddir, "Ex_%s.h5" % it), "field")

    print("Data loaded.")

    nbg = nbge + nbgi
    nj  = nje + nji
    ntot = nbg + nj

    ntot = ntot.swapaxes(0, 2)
    nbg  = nbg.swapaxes(0, 2)
    nj   = nj.swapaxes(0, 2)
    Bx   = Bx.swapaxes(0, 2)
    By   = By.swapaxes(0, 2)
    Bz   = Bz.swapaxes(0, 2)
    Ex   = Ex.swapaxes(0, 2)
    Ey   = Ey.swapaxes(0, 2)
    Ez   = Ez.swapaxes(0, 2)

    # Densities have one fewer cell in each dimension with respect to fields
    # and the grid. Fix this by shaving off the last cell in the quantities
    # with extra data.
    Bx = Bx[:-1,:-1,:-1]
    By = By[:-1,:-1,:-1]
    Bz = Bz[:-1,:-1,:-1]
    Ex = Ex[:-1,:-1,:-1]
    Ey = Ey[:-1,:-1,:-1]
    Ez = Ez[:-1,:-1,:-1]

    #===========================================================================
    # Reduce resolution
    # To avoid averaging noisy quantities, it is probably best to reduce the
    # resolution before calculating derived quantities below (rather than
    # first calculating derived quantities and then reducing the grid)
    if rrp2 > 0:
        print("Cutting resolution by 2^%d" % rrp2)

    for reduction_step in range(rrp2):
        Bx   = half_res(Bx)
        By   = half_res(By)
        Bz   = half_res(Bz)
        Ex   = half_res(Ex)
        Ey   = half_res(Ey)
        Ez   = half_res(Ez)
        ntot = half_res(ntot)
        nbg  = half_res(nbg)
        nj   = half_res(nj)

        x = x[::2]
        y = y[::2]
        z = z[::2]
        dx = 2 * dx
        dy = 2 * dy
        dz = 2 * dz

    X, Y, Z = np.meshgrid(x, y, z, indexing = "ij")

    #===========================================================================
    # Derived fields

    print("Calculating derived fields...")
    Sz = c * (Ex * By - Ey * Bx) / (4 * np.pi)
    Bsq = Bx * Bx + By * By + Bz * Bz
    sigma = np.zeros(Bsq.shape)
    ixnz = np.where(ntot > 0.0)
    sigma[ixnz] = Bsq[ixnz] / (4 * np.pi * me * c * c * ntot[ixnz])

    #===========================================================================
    # Initialize vtk arrays + point grid
    print("Creating the grid...")

    dims = Bx.shape
    idata = vtkImageData()
    idata.SetDimensions(dims)    
    idata.SetOrigin(xmin,ymin,zmin)
    idata.SetSpacing(dx,dy,dz)

    # data_points = algs.make_vector(
    #     X.flatten('F'), Y.flatten('F'), Z.flatten('F')
    # )
    # points = vtkPoints()
    # points.Allocate(dims[0] * dims[1] * dims[2])
    # points.SetData(dsa.numpyTovtkDataArray(data_points,'Points'))

    np_Bfield = algs.make_vector(
        Bx.flatten('F') / B0, By.flatten('F') / B0, Bz.flatten('F') / B0
    )
    Bfield = dsa.numpyTovtkDataArray(np_Bfield, "B field")

    density = dsa.numpyTovtkDataArray(ntot.flatten('F') / ngj, "Density")

    magnetization = dsa.numpyTovtkDataArray(sigma.flatten('F'), "Cold sigma")

    poynting = dsa.numpyTovtkDataArray(Sz.flatten('F') / Sznorm, "Poynting z")

    idata.GetPointData().SetVectors(Bfield)
    idata.GetPointData().AddArray(density)
    idata.GetPointData().AddArray(poynting)
    idata.GetPointData().AddArray(magnetization)

    print("Grid created")

    #===========================================================================
    # Save to vtk file
    print("Writing file...")

    # Writing data
    if not os.path.isdir("./vtk"):
        os.mkdir("./vtk")
    # writer = vtkStructuredGridWriter()
    # writer = vtkImageWriter()
    writer = vtkXMLImageDataWriter()
    # writer.SetFileTypeToBinary()
    it = int(it)
    writer.SetFileName(f"vtk/data_{it:06d}.xml")
    writer.SetInputData(idata)
    writer.Write()
    # h5file.close()

