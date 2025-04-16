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

def h5tovtk(it, if_einf=False):

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
    Sznorm = B0**2 * zmono**4 * c / (4 * np.pi)

    #===========================================================================
    # Grid reconstruction

    # Nph, Nth, Nr = np.array(
    #     h5py.File(dir + f"/densities/densities_{it}.h5")["electrons"]
    # ).shape
    h5file = h5py.File(datadir + f"/densities/mapxyz_electrons_bg_{it}.h5")
    Nz, Ny, Nx = np.array(h5file["field"]).shape
    h5file.close()

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
    x = np.linspace(xmin, xmax, Nx)
    y = np.linspace(ymin, ymax, Ny)
    z = np.linspace(zmin, zmax, Nz)
    # R, Th, Ph = np.meshgrid(r, th, ph, indexing="ij")
    # costh, sinth = np.cos(Th), np.sin(Th)
    # cosph, sinph = np.cos(Ph), np.sin(Ph)
    X, Y, Z = np.meshgrid(x, y, z, indexing = "ij")

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

    print("Data loaded.")

    #===========================================================================
    # Initialize vtk arrays + point grid
    print("Creating the grid...")

    dims = Bx.shape
    # dims = Bru.shape
    sgrid = vtkStructuredGrid()
    sgrid.SetDimensions(dims)    

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
    writer.SetFileTypeToBinary()
    it = int(it)
    writer.SetFileName(f"vtk/data_{it:06d}.vtk")
    writer.SetInputData(sgrid)
    writer.Write()
    h5file.close()
