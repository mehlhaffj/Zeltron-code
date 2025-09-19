!***********************************************************************!
!                       The Zeltron code project.                       !
!***********************************************************************!
! Copyright (C) 2012-2015. Authors: Benoît Cerutti & Greg Werner        !
!                                                                       !
! This program is free software: you can redistribute it and/or modify  !
! it under the terms of the GNU General Public License as published by  !
! the Free Software Foundation, either version 3 of the License, or     !
! (at your option) any later version.                                   !
!                                                                       !
! This program is distributed in the hope that it will be useful,       !
! but WITHOUT ANY WARRANTY; without even the implied warranty of        !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         !
! GNU General Public License for more details.                          !
!                                                                       !
! You should have received a copy of the GNU General Public License     !
! along with this program. If not, see <http://www.gnu.org/licenses/>.  !
!***********************************************************************!

MODULE MOD_INPUT

IMPLICIT NONE

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++ CONSTANTS +++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Speed of light [cm/s]
DOUBLE PRECISION, PARAMETER, PUBLIC           :: c=299792458d2
! Fundamental charge [esu]
DOUBLE PRECISION, PARAMETER, PUBLIC           :: e=4.8032068d-10
! Boltzmann constant [erg/K]
DOUBLE PRECISION, PARAMETER, PUBLIC           :: k=1.380658d-16
! Planck constant [erg.s]
DOUBLE PRECISION, PARAMETER, PUBLIC           :: h=6.6260755d-27
! Mass of the electron [g]
DOUBLE PRECISION, PARAMETER, PUBLIC           :: me=9.1093897d-28
! Mass of the proton [g]
DOUBLE PRECISION, PARAMETER, PUBLIC           :: mp=1.6726231d-24
! 1 eV in erg [erg]
DOUBLE PRECISION, PARAMETER, PUBLIC           :: evtoerg=1.602177d-12
! pi
DOUBLE PRECISION, PARAMETER, PUBLIC           :: pi=dacos(-1d0)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++ DIMENSION ++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! spatial dimension
INTEGER, PARAMETER, PUBLIC :: NDIM=3

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++ SAVE/RESTORE ++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Save particle and field data for future restoration of the simulation
LOGICAL, PARAMETER, PUBLIC :: CHECKPOINT=.FALSE.

! Restore a simulation where it stopped
LOGICAL, PARAMETER, PUBLIC :: RESTORE=.FALSE.

! Saving frequency, in seconds (the simulation will try to save its 
! current state *before* elapsed time n*FSAVE, for integer n).
! N.B. It is highly recommended to have more than 1 checkpoint during
! the simulation, unless it is certain that the simulation will end
! naturally (as opposed to being killed before it is finished).
DOUBLE PRECISION, PARAMETER, PUBLIC :: FSAVE=0.

! Give the time step from which the simulation should restart
INTEGER, PARAMETER, PUBLIC :: time_ref=0

! Whether to test restore capability: this restores data from step 
!   time_ref, and then immediately write data (as if dumping
!   for a checkpoint) under the next time step; the simulation then halts.
LOGICAL, PARAMETER, PUBLIC :: TEST_RESTORE = .FALSE.

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++ BOUNDARY CONDITIONS ++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Specify the boundary conditions for the fields:
! 1. "PERIODIC": Periodic boundary conditions
! 2. "METAL": Perfect metal with infinite conductivity
! 3. "NOZZLE": Impose B-line rotation across a conducting nozzle
!   NOZZLE only supported for [X|Y|Z]MIN boundaries !!
CHARACTER(LEN=10), PARAMETER, PUBLIC :: BOUND_FIELD_XMIN="OPEN"
CHARACTER(LEN=10), PARAMETER, PUBLIC :: BOUND_FIELD_XMAX="OPEN"
CHARACTER(LEN=10), PARAMETER, PUBLIC :: BOUND_FIELD_YMIN="OPEN"
CHARACTER(LEN=10), PARAMETER, PUBLIC :: BOUND_FIELD_YMAX="OPEN"
CHARACTER(LEN=10), PARAMETER, PUBLIC :: BOUND_FIELD_ZMIN="NOZZLE"
CHARACTER(LEN=10), PARAMETER, PUBLIC :: BOUND_FIELD_ZMAX="OPEN"

! Specify the boundary conditions for the particles:
! 1. "PERIODIC": Periodic boundary conditions
! 2. "REFLECT": Particles are elastically reflected at the wall
! 3. "ABSORB": Particles are absorbed at the wall
! 4. "INJECT": Inject particles across the NOZZLE where field-
!     -line rotation is imposed
CHARACTER(LEN=10), PARAMETER, PUBLIC :: BOUND_PART_XMIN="PERIODIC"
CHARACTER(LEN=10), PARAMETER, PUBLIC :: BOUND_PART_XMAX="PERIODIC"
CHARACTER(LEN=10), PARAMETER, PUBLIC :: BOUND_PART_YMIN="PERIODIC"
CHARACTER(LEN=10), PARAMETER, PUBLIC :: BOUND_PART_YMAX="PERIODIC"
CHARACTER(LEN=10), PARAMETER, PUBLIC :: BOUND_PART_ZMIN="PERIODIC"
CHARACTER(LEN=10), PARAMETER, PUBLIC :: BOUND_PART_ZMAX="PERIODIC"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++++ INITIAL CONDITIONS +++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 1. "RECONN": Initial reconnection fields
! 2. "UNIFORM": Initial uniform \vec{B} = B_0 \hat{z}
! 3. "MONOPOLE": Place magnetic monopole on z-axis below zmin
!     Strength of B-field at (xmiddle, ymiddle, zmin) is B_0
CHARACTER(LEN=10), PARAMETER, PUBLIC :: INIT="MONOPOLE"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++ INPUT PARAMETERS ++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Number of cells in X
INTEGER*8, PARAMETER, PUBLIC :: NCX=64

! Number of cells in Y
INTEGER*8, PARAMETER, PUBLIC :: NCY=NCX

! Number of cells in Z
INTEGER*8, PARAMETER, PUBLIC :: NCZ=64

! Number of particles per cell per species
INTEGER*8, PARAMETER, PUBLIC :: PPC=1

! Number of process (domain decomposition in the X- Y- and Z-directions)
INTEGER, PARAMETER, PUBLIC :: NPX=1
INTEGER, PARAMETER, PUBLIC :: NPY=1
INTEGER, PARAMETER, PUBLIC :: NPZ=1

! Mass ratio IONS/ELECTRONS
DOUBLE PRECISION, PARAMETER, PUBLIC :: mass_ratio=1d0

! Spatial boundaries in the X-direction
! _stat quantities can be used to set other PARAMETERS in this input file
! xpml_stat is the fraction of xmax-xmin that will be used for the x-pml
DOUBLE PRECISION, PARAMETER, PUBLIC :: xmin_stat=-16d0
DOUBLE PRECISION, PARAMETER, PUBLIC :: xmax_stat=16d0
DOUBLE PRECISION, PARAMETER, PUBLIC :: xpml_stat=0.1
DOUBLE PRECISION, PUBLIC            :: xmin=xmin_stat
DOUBLE PRECISION, PUBLIC            :: xmax=xmax_stat
DOUBLE PRECISION, PUBLIC            :: xpml1=xpml_stat
DOUBLE PRECISION, PUBLIC            :: xpml2=xpml_stat
! DOUBLE PRECISION, PARAMETER, PUBLIC :: xmin=-32d0,xmax=32d0

! Spatial boundaries in the Y-direction
! _stat quantities can be used to set other PARAMETERS in this input file
! ypml_stat is the fraction of xmax-xmin that will be used for the y-pml
DOUBLE PRECISION, PARAMETER, PUBLIC :: ymin_stat=xmin_stat
DOUBLE PRECISION, PARAMETER, PUBLIC :: ymax_stat=xmax_stat
DOUBLE PRECISION, PARAMETER, PUBLIC :: ypml_stat=xpml_stat
DOUBLE PRECISION, PUBLIC            :: ymin=ymin_stat
DOUBLE PRECISION, PUBLIC            :: ymax=ymax_stat
DOUBLE PRECISION, PUBLIC            :: ypml1=ypml_stat
DOUBLE PRECISION, PUBLIC            :: ypml2=ypml_stat
! DOUBLE PRECISION, PARAMETER, PUBLIC :: ymin=xmin,ymax=xmax

! Spatial boundaries in the Z-direction
! Spatial boundaries in the Z-direction
! _stat quantities can be used to set other PARAMETERS in this input file
! zpml_stat is the fraction of xmax-xmin that will be used for the z-pml
DOUBLE PRECISION, PARAMETER, PUBLIC :: zmin_stat=0d0
DOUBLE PRECISION, PARAMETER, PUBLIC :: zmax_stat=(xmax_stat-xmin_stat)*NCZ/NCX
DOUBLE PRECISION, PARAMETER, PUBLIC :: zpml_stat=xpml_stat
DOUBLE PRECISION, PUBLIC            :: zmin=zmin_stat
DOUBLE PRECISION, PUBLIC            :: zmax=zmax_stat
DOUBLE PRECISION, PUBLIC            :: zpml1=zpml_stat
DOUBLE PRECISION, PUBLIC            :: zpml2=zpml_stat
! DOUBLE PRECISION, PARAMETER, PUBLIC :: zmin=0d0,zmax=(xmax-xmin)*NCZ/NCX

! Dump data frequency in terms of timesteps
INTEGER, PARAMETER, PUBLIC :: FDUMP=10

! Number of data dumps
INTEGER, PARAMETER, PUBLIC :: NDUMP=20

! Number of time steps
INTEGER, PARAMETER, PUBLIC :: NT=FDUMP*NDUMP

! Poisson solver calling frequency in terms of timesteps
INTEGER, PARAMETER, PUBLIC :: FREQ_POISSON=25

! Number of iterations to solve Poisson's equation
INTEGER, PARAMETER, PUBLIC :: NIT=500

! Switch ON/OFF the radiation reaction force (synchrotron radiation and
! inverse Compton scattering in the Thomson regime only)
LOGICAL, PARAMETER, PUBLIC :: RAD_FORCE=.FALSE.

! Switch ON/OFF spatial filtering the electric field
LOGICAL, PARAMETER, PUBLIC :: FILTER=.FALSE.

! Smoothing parameter
DOUBLE PRECISION, PARAMETER, PUBLIC :: alpha=1d0/64d0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Number of grid cells in X
INTEGER, PARAMETER, PUBLIC :: NX=NCX+1

! Number of grid cells in Y
INTEGER, PARAMETER, PUBLIC :: NY=NCY+1

! Number of grid cells in Z
INTEGER, PARAMETER, PUBLIC :: NZ=NCZ+1

! Number of cell per process in the X-direction
INTEGER, PARAMETER, PUBLIC :: NCXP=NCX/NPX

! Number of cell per process in the Y-direction
INTEGER, PARAMETER, PUBLIC :: NCYP=NCY/NPY

! Number of cell per process in the Z-direction
INTEGER, PARAMETER, PUBLIC :: NCZP=NCZ/NPZ

! Number of grid cells per process in X
INTEGER, PARAMETER, PUBLIC :: NXP=NCXP+1

! Number of grid cells per process in Y
INTEGER, PARAMETER, PUBLIC :: NYP=NCYP+1

! Number of grid cells per process in Z
INTEGER, PARAMETER, PUBLIC :: NZP=NCZP+1

! Total number of particles per species NP=NCX*NCY*NCZ*PPC
INTEGER*8, PARAMETER, PUBLIC :: NP=NCX*NCY*NCZ*PPC

! Spatial step
DOUBLE PRECISION, PARAMETER, PUBLIC :: dx=(xmax_stat-xmin_stat)/NCX
DOUBLE PRECISION, PARAMETER, PUBLIC :: dy=(ymax_stat-ymin_stat)/NCY
DOUBLE PRECISION, PARAMETER, PUBLIC :: dz=(zmax_stat-zmin_stat)/NCZ

! Time step (Courant-Friedrichs-Lewy timestep)
DOUBLE PRECISION, PARAMETER, PUBLIC :: dt=0.99*1d0/sqrt(1d0/dx**2d0+&
                                          1d0/dy**2d0+1d0/dz**2d0)/c

! Nodal grid in each domain
DOUBLE PRECISION, DIMENSION(1:NXP) :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP) :: ygp
DOUBLE PRECISION, DIMENSION(1:NZP) :: zgp

! Yee grid in each domain
DOUBLE PRECISION, DIMENSION(1:NXP) :: xyeep
DOUBLE PRECISION, DIMENSION(1:NYP) :: yyeep
DOUBLE PRECISION, DIMENSION(1:NZP) :: zyeep

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Initial (tearing-mode) perturbation amplitude
DOUBLE PRECISION, PARAMETER, PUBLIC :: perturb_amp=0.0

! Ratio of Bz/B0 (Guide field / reconnecting field)
DOUBLE PRECISION, PARAMETER, PUBLIC :: guide_field=0.0
                                          
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Temperature in unit of me*c^2/k of the DRIFTING ELECTRONS in the co-moving frame
DOUBLE PRECISION, PARAMETER, PUBLIC :: thde=1d0

! Temperature in unit of mi*c^2/k of the DRIFTING IONS in the co-moving frame
DOUBLE PRECISION, PARAMETER, PUBLIC :: thdi=1d0

! Temperature in unit of me*c^2/k of the BACKGROUND ELECTRONS in the lab frame
DOUBLE PRECISION, PARAMETER, PUBLIC :: thbe=1d0

! Temperature in unit of mi*c^2/k of the BACKGROUND IONS in the lab frame
DOUBLE PRECISION, PARAMETER, PUBLIC :: thbi=1d0

! Minimum Larmor radius of the electrons
! For INIT="RECONN", this is thde*me*c^2/(e*B0)
! For INIT="UNIFORM" or "MONOPOLE", this is me*c^2/(e*B0)
DOUBLE PRECISION, PARAMETER, PUBLIC :: rhoc=1d0

! Drift velocity betad=vdrift/c
! For "INJECT" particle BC, this is the drift speed of the injected particles
DOUBLE PRECISION, PARAMETER, PUBLIC :: betad=0.5d0

! Ratio of initial BACKGROUND to DRIFTING particle number densities
DOUBLE PRECISION, PARAMETER, PUBLIC :: ratio_nb2nd=1.0/SQRT(2d0**4)

! Energy density ratio between external radiation field and the magnetic field
! udens_ratio=Uph/Ub, where Ub=B0^2/8*pi
DOUBLE PRECISION, PARAMETER, PUBLIC :: udens_ratio=0d0

! 4-velocity boundaries for the particles
DOUBLE PRECISION, PARAMETER, PUBLIC :: umin=1d-2,umax=1d2

! 4-velocity boundaries for the drifting electrons
DOUBLE PRECISION, PARAMETER, PUBLIC :: udemin=-1d2,udemax=1d2
! 4-velocity boundaries for the drifting ions
DOUBLE PRECISION, PARAMETER, PUBLIC :: udpmin=-1d2,udpmax=1d2

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Extra parameters for INIT="MONOPOLE" or "UNIFORM"
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! The magnetization injected across NOZZLE with INJECT particle BC
! This is the total magnetization (including electrons and positrons)
DOUBLE PRECISION, PARAMETER, PUBLIC :: sigma_inj = 1.0/(ratio_nb2nd**2)
! The multiplicity injected across NOZZLE with INJECT particle BC
! This 
DOUBLE PRECISION, PARAMETER, PUBLIC :: kappa  = 12.0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Extra parameters for BOUND_FIELD_[X|Y|Z]MIN = "NOZZLE"
! This might only be implemented for BOUND_FIELD_ZMIN="NOZZLE"
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Convenience parameter: magnetization of background plasma
DOUBLE PRECISION, PARAMETER, PRIVATE :: sigma_bg = sigma_inj/ratio_nb2nd
! Light cylinder radius defined by nozzle solid-body rotation
DOUBLE PRECISION, PARAMETER, PUBLIC :: rlc=12*rhoc  ! 2*kappa*sigma_bg*rhoc
! Nozzle radius
! (in plane perp. to the coord. for which BOUND_FIELD_MIN = "NOZZLE")
DOUBLE PRECISION, PARAMETER, PUBLIC :: rnozzle=0.9*rlc
! Nozzle is shut off at rnozzle across half-thickness delta_nozzle
! Nozzle ang. freq. is omega = 0.5*c/rlc*(1.0-TANH((R-rnozzle)/delta_nozzle))
DOUBLE PRECISION, PARAMETER, PUBLIC :: delta_nozzle=0.1*rnozzle

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Extra parameters for BOUND_PART_[X|Y|Z]MIN="INJECT"
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Plasma injection rate
! Changes the physical number of particles injected per unit time without
! modifying the number of macroparticles injected.
! Hence, this changes physics but not simulation speed, memory cost, or
! particle statistics.
DOUBLE PRECISION, PARAMETER, PUBLIC :: rate=1.0
! A factor by which to amplify injection of macroparticles without changing
! the physical particle density.
! Hence, this does NOT change physics, but may change simulation speed, memory
! cost, and particle statistics.
DOUBLE PRECISION, PARAMETER, PUBLIC :: over_inject=1.0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Extra parameters for INIT="MONOPOLE"
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! The location, (0,0,zmp), where the monopole is placed
DOUBLE PRECISION, PARAMETER, PUBLIC :: zmp=zmin_stat-rlc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!++++++++++++++++++++++ DEFINE PARTICLE ARRAYS +++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! pcl_species(1,:)=x      : x position of the particle
! pcl_species(2,:)=y      : y position of the particle
! pcl_species(3,:)=z      : z position of the particle
! pcl_species(4,:)=ux     : x-component of the reduced particle 4-velocity
! pcl_species(5,:)=uy     : y-component of the reduced particle 4-velocity
! pcl_species(6,:)=uz     : z-component of the reduced particle 4-velocity
! pcl_species(7,:)=weight : weight of the particle
!
! Defition of the particle distribution data function: it gives extra
! information regarding the particles, e.g., electric/magnetic field/force
! at the location of the particle.
!
! pcl_data_species(1,:)=Ell : Electric field parallel to the particle velocity
! pcl_data_species(2,:)=Bpp : Magnetic field perp. to the particle velocity
! pcl_data_species(3,:)=Fr  : Synchrotron radiation reaction force
! pcl_data_species(4,:)=Fe  : Electric force
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! DRIFTING ELECTRONS distribution function components
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: pcl_ed(:,:)
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: pcl_data_ed(:,:)
INTEGER*8, ALLOCATABLE, PUBLIC        :: taged(:)

! DRIFTING IONS distribution function components
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: pcl_pd(:,:)
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: pcl_data_pd(:,:)
INTEGER*8, ALLOCATABLE, PUBLIC        :: tagpd(:)

! BACKGROUND ELECTRONS distribution function components
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: pcl_eb(:,:)
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: pcl_data_eb(:,:)
INTEGER*8, ALLOCATABLE, PUBLIC        :: tageb(:)

! BACKGROUND IONS distribution function components
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: pcl_pb(:,:)
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: pcl_data_pb(:,:)
INTEGER*8, ALLOCATABLE, PUBLIC        :: tagpb(:)

! Temporary particle distribution function
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: pcl_f(:,:)
DOUBLE PRECISION, ALLOCATABLE, PUBLIC :: pcl_data_f(:,:)
INTEGER*8, ALLOCATABLE, PUBLIC        :: tagf(:)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++ INITIAL FIELD ARRAYS ++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! INITIAL Magnetic and Electric fields components Yee lattice
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP) :: Bx0,By0,Bz0
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP) :: Ex0,Ey0,Ez0

! INITIAL Magnetic and Electric fields components at nodes
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP) :: Bxg0,Byg0,Bzg0
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP) :: Exg0,Eyg0,Ezg0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++ ANALYSIS ++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Data formatting (for data in text files) for integers and doubles
CHARACTER(LEN=6), PARAMETER, PUBLIC  :: FMT_INT="I26"
CHARACTER(LEN=11), PARAMETER, PUBLIC :: FMT_DBL="E26.16E3"

! For text output of fields/density/pressure/current: (irrelevant for hdf5 output)
! If .TRUE., each MPI process will write to a separate file.
! This can be faster (for large simulations), but may create too many files.
LOGICAL, PARAMETER, PUBLIC :: DUMP_FIELDS_PARALLEL=.FALSE.

! Switch ON/OFF write fields to disk
LOGICAL, PARAMETER, PUBLIC :: WRITE_FIELDS=.TRUE.

! Switch ON/OFF write the charge and current densities to disk
LOGICAL, PARAMETER, PUBLIC :: WRITE_RHOJ=.TRUE.

! For text output of particle data: (irrelevant for hdf5 output)
! If .TRUE., each MPI process will write to its own file.
! This can be faster (for large simulations), but may create too many files.
LOGICAL, PARAMETER, PUBLIC :: DUMP_PARTICLES_PARALLEL=.FALSE.

! Switch ON/OFF write raw particle data to disk
LOGICAL, PARAMETER, PUBLIC :: WRITE_PARTICLES=.FALSE.

! Switch ON/OFF particle tracker
LOGICAL, PARAMETER, PUBLIC :: TRACK_PARTICLES=.FALSE.

! Switch ON/OFF analysis of particle energy and angular distributions
LOGICAL, PARAMETER, PUBLIC :: ANALYZE_DISTRIBUTIONS=.TRUE.

! Switch ON/OFF analysis the plasma density
LOGICAL, PARAMETER, PUBLIC :: ANALYZE_DENSITIES=.TRUE.

! Switch ON/OFF analysis the relativistic macroscopic fluid quantities
LOGICAL, PARAMETER, PUBLIC :: ANALYZE_FLUID=.FALSE.
! Finer granularity of fluid quantity analysis
! ANALYZE_FLUID takes precedence over these quantities and must be TRUE for
! any of them to have an effect.
! Choose from among the possible fluid quantities to dump by toggling
! the following.
LOGICAL, PARAMETER, PUBLIC :: ANALYZE_EN_DENS=.TRUE.
LOGICAL, PARAMETER, PUBLIC :: ANALYZE_PTCL_FLUX=.FALSE.
LOGICAL, PARAMETER, PUBLIC :: ANALYZE_MOM_DENS=.FALSE.
LOGICAL, PARAMETER, PUBLIC :: ANALYZE_PRESSURE=.FALSE.

! Switch ON/OFF analysis the radiation spectra and anisotropies
LOGICAL, PARAMETER, PUBLIC :: ANALYZE_RADIATION=.FALSE.

! Total number of tracked particles
INTEGER, PARAMETER, PUBLIC :: NSAMPLE=100

! Give the time step from which the tracker should restart
INTEGER, PARAMETER, PUBLIC :: time_track=1

! The tracking frequency (timesteps between recording of tracking data)
! E.g., 1 = record every step; 2 = record every other step; ....
INTEGER, PARAMETER, PUBLIC :: FTRACK=1

! Particle's Lorentz factor grid for the spectrum
INTEGER, PARAMETER, PUBLIC :: NU=100

! Spatial grid size for dumping EM-related fields (E,B,J)
! Spatial grid size for dumping EM-related fields
! If NFIELDX = NCX, etc., Zeltron dumps exactly the (nodal) field it uses.
INTEGER, PARAMETER, PUBLIC :: NFIELDX=NCX,NFIELDY=NCY,NFIELDZ=NCZ

! Spatial grid size for dumping particle-related fields (densities, etc.)
INTEGER, PARAMETER, PUBLIC :: NDX=NCX,NDY=NCY,NDZ=NCZ

! Angular distribution grid
INTEGER, PARAMETER, PUBLIC :: NPHI=50,NLBA=90

! Synchrotron radiation frequency
INTEGER, PARAMETER, PUBLIC :: NNU=100

! Frequency boundaries for the synchrotron radiation
DOUBLE PRECISION, PARAMETER, PUBLIC :: nmin=1d9,nmax=1d14

! Latitude boundaries for the angular distribution (phi)
DOUBLE PRECISION, PARAMETER, PUBLIC :: pmin=-89.99d0,pmax=89.99d0

! Longitude boundaries for the angular distribution (lambda)
DOUBLE PRECISION, PARAMETER, PUBLIC :: lmin=-179.99d0,lmax=179.99d0

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++ TUNING ++++++++++++++++++++++++++++++++++++
! You might want to adjust these parameters, but they should mostly  
! be left as they are. 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! (Inverse) frequency, in steps, with which all ranks communicate
! general info, such as whether to checkpoint, or the contents of the
! command file.
INTEGER, PARAMETER, PUBLIC :: FGENCOMM=20

! Fraction of FSAVE by which a checkpoint time can be shifted (earlier)
! so that it occurs as the same time as a dump.
DOUBLE PRECISION, PARAMETER, PUBLIC :: FSAVE_SHIFT_FRAC=0.25

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++ DEBUGGING +++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! For normal runs, RANDOMIZE = true ensures a different seed each run
! For testing new code, can set to false for reproducibility
LOGICAL, PARAMETER, PUBLIC :: RANDOMIZE=.TRUE.

! whether to write one diagnostic file per rank
LOGICAL, PARAMETER, PUBLIC :: writePerRankFiles=.FALSE.

! the name of the per-rank file
CHARACTER(len=18), PUBLIC  :: perRankFile

! Do a simulation without particles
LOGICAL, PARAMETER, PUBLIC :: VACUUM_SIM=.FALSE.

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE MOD_INPUT
