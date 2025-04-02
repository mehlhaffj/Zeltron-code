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

MODULE MOD_INITIAL

USE MOD_INPUT

IMPLICIT NONE

PRIVATE

PUBLIC :: COM_TOPOLOGY ! Initialize the virtual topology of the processes
PUBLIC :: SET_FIELDS ! Generation of the initial fields E and B
PUBLIC :: WEIGHT ! This subroutine weight each super-particles
PUBLIC :: SET_TAG ! This subroutine assigns a tag to each super-particles
PUBLIC :: SET_MAXWELLIAN ! Generation of a relativistic Maxwellian distribution function
PUBLIC :: DIST_MAXWELL ! Integration of the Maxwell distribution function
PUBLIC :: GEN_U ! Generates the particle 4-velocity
PUBLIC :: SET_DRIFT_MAXWELLIAN ! Generation of a drifting Maxwellian distribution function
PUBLIC :: INIT_DRIFT_MAXWELLIAN ! Computes the cumulative distribution functions
PUBLIC :: GEN_UP ! Generates particle parallel momentum
PUBLIC :: GEN_PS ! Generates particle perpendicular momentum
PUBLIC :: INIT_RANDOM_SEED ! Avoid repeating series of random numbers
PUBLIC :: INJ_DRIFT_MAXWELLIAN ! Inject drifting Maxwellian from simulation boundary

 CONTAINS

!***********************************************************************
! This subroutine initializes a virtual cartesian topology of the processes.
!
! INPUT:
!
! - ngh: neighbor array
! - coords: coordinates of the process id in the cartesian topology
! - NPROC: Total number of processes
! - id: process rank
! - COMM: communicator
! - ierr: error code
!
! OUTPUT: ngh,coords,COMM
!
!         -BACK-             -CENTER-             -FRONT-
!
!    BNW----BN----BNE      NW----N-----NE     FNW----FN----FNE
!     |     |      |       |     |     |       |     |      |
!    BW-----B------BE      W-----id----E      FW-----F------FE
!     |     |      |       |     |     |       |     |      |
!    BSW----BS----BSE      SW----S-----SE     FSW----FS----FSE
!
!***********************************************************************

SUBROUTINE COM_TOPOLOGY(ngh,coords,NPROC,id,COMM,ierr)

IMPLICIT NONE

INCLUDE 'mpif.h'

INTEGER, DIMENSION(MPI_STATUS_SIZE) :: stat
INTEGER                             :: NPROC,COMM,direction,step,ierr,id
LOGICAL                             :: reorder

INTEGER, PARAMETER                  :: BNorth=1,BEast=2,BSouth=3,BWest=4
INTEGER, PARAMETER                  :: BNEast=5,BSEast=6,BSWest=7,BNWest=8
INTEGER, PARAMETER                  :: BCenter=9

INTEGER, PARAMETER                  :: North=10,East=11,South=12,West=13
INTEGER, PARAMETER                  :: NEast=14,SEast=15,SWest=16,NWest=17

INTEGER, PARAMETER                  :: FNorth=18,FEast=19,FSouth=20,FWest=21
INTEGER, PARAMETER                  :: FNEast=22,FSEast=23,FSWest=24,FNWest=25
INTEGER, PARAMETER                  :: FCenter=26

INTEGER, DIMENSION(3)               :: dims,coords
LOGICAL, DIMENSION(3)               :: periods
INTEGER, DIMENSION(26)              :: ngh

!***********************************************************************

! Initialization of the cartesian topology
periods(1)=.TRUE.
periods(2)=.TRUE.
periods(3)=.TRUE.
reorder=.FALSE.
dims(1)=NPX
dims(2)=NPY
dims(3)=NPZ

! Creation of the dimension in each direction
CALL MPI_DIMS_CREATE(NPROC,3,dims,ierr)

! Creation of the topology
CALL MPI_CART_CREATE(MPI_COMM_WORLD,3,dims,periods,reorder,COMM,ierr)

! To obtain the ID number of each process
CALL MPI_COMM_RANK(COMM,id,ierr)

! To obtain the coordinates of the process
CALL MPI_CART_COORDS(COMM,id,3,coords,ierr)

!PRINT*,'*******************************************************'
!PRINT*,'Processus: ',id,', Coordinates: ',coords(1),coords(2),coords(3)
!PRINT*,'*******************************************************'

! Search neighbors

! West-East
direction=0
step=1

CALL MPI_CART_SHIFT(COMM,direction,step,ngh(West),ngh(East),ierr)

! North-South
direction=1
step=1

CALL MPI_CART_SHIFT(COMM,direction,step,ngh(South),ngh(North),ierr)

! Back-Front
direction=2
step=1

CALL MPI_CART_SHIFT(COMM,direction,step,ngh(BCenter),ngh(FCenter),ierr)

! Diagonal neighbors

! South-West
CALL MPI_SENDRECV(ngh(South),1,MPI_INTEGER,ngh(East),1,ngh(SWest),1,&
MPI_INTEGER,ngh(West),1,COMM,stat,ierr)

! North-West
CALL MPI_SENDRECV(ngh(North),1,MPI_INTEGER,ngh(East),2,ngh(NWest),1,&
MPI_INTEGER,ngh(West),2,COMM,stat,ierr)

! South-East
CALL MPI_SENDRECV(ngh(South),1,MPI_INTEGER,ngh(West),3,ngh(SEast),1,&
MPI_INTEGER,ngh(East),3,COMM,stat,ierr)

! North-East
CALL MPI_SENDRECV(ngh(North),1,MPI_INTEGER,ngh(West),4,ngh(NEast),1,&
MPI_INTEGER,ngh(East),4,COMM,stat,ierr)

! Back-North
CALL MPI_SENDRECV(ngh(North),1,MPI_INTEGER,ngh(FCenter),1,ngh(BNorth),1,&
MPI_INTEGER,ngh(BCenter),1,COMM,stat,ierr)

! Back-South
CALL MPI_SENDRECV(ngh(South),1,MPI_INTEGER,ngh(FCenter),1,ngh(BSouth),1,&
MPI_INTEGER,ngh(BCenter),1,COMM,stat,ierr)

! Back-East
CALL MPI_SENDRECV(ngh(East),1,MPI_INTEGER,ngh(FCenter),1,ngh(BEast),1,&
MPI_INTEGER,ngh(BCenter),1,COMM,stat,ierr)

! Back-West
CALL MPI_SENDRECV(ngh(West),1,MPI_INTEGER,ngh(FCenter),1,ngh(BWest),1,&
MPI_INTEGER,ngh(BCenter),1,COMM,stat,ierr)

! Front-North
CALL MPI_SENDRECV(ngh(North),1,MPI_INTEGER,ngh(BCenter),1,ngh(FNorth),1,&
MPI_INTEGER,ngh(FCenter),1,COMM,stat,ierr)

! Front-South
CALL MPI_SENDRECV(ngh(South),1,MPI_INTEGER,ngh(BCenter),1,ngh(FSouth),1,&
MPI_INTEGER,ngh(FCenter),1,COMM,stat,ierr)

! Front-East
CALL MPI_SENDRECV(ngh(East),1,MPI_INTEGER,ngh(BCenter),1,ngh(FEast),1,&
MPI_INTEGER,ngh(FCenter),1,COMM,stat,ierr)

! Front-West
CALL MPI_SENDRECV(ngh(West),1,MPI_INTEGER,ngh(BCenter),1,ngh(FWest),1,&
MPI_INTEGER,ngh(FCenter),1,COMM,stat,ierr)

! Back-NEast
CALL MPI_SENDRECV(ngh(NEast),1,MPI_INTEGER,ngh(FCenter),1,ngh(BNEast),1,&
MPI_INTEGER,ngh(BCenter),1,COMM,stat,ierr)

! Back-SEast
CALL MPI_SENDRECV(ngh(SEast),1,MPI_INTEGER,ngh(FCenter),1,ngh(BSEast),1,&
MPI_INTEGER,ngh(BCenter),1,COMM,stat,ierr)

! Back-SWest
CALL MPI_SENDRECV(ngh(SWest),1,MPI_INTEGER,ngh(FCenter),1,ngh(BSWest),1,&
MPI_INTEGER,ngh(BCenter),1,COMM,stat,ierr)

! Back-NWest
CALL MPI_SENDRECV(ngh(NWest),1,MPI_INTEGER,ngh(FCenter),1,ngh(BNWest),1,&
MPI_INTEGER,ngh(BCenter),1,COMM,stat,ierr)

! Front-NEast
CALL MPI_SENDRECV(ngh(NEast),1,MPI_INTEGER,ngh(BCenter),1,ngh(FNEast),1,&
MPI_INTEGER,ngh(FCenter),1,COMM,stat,ierr)

! Front-SEast
CALL MPI_SENDRECV(ngh(SEast),1,MPI_INTEGER,ngh(BCenter),1,ngh(FSEast),1,&
MPI_INTEGER,ngh(FCenter),1,COMM,stat,ierr)

! Front-SWest
CALL MPI_SENDRECV(ngh(SWest),1,MPI_INTEGER,ngh(BCenter),1,ngh(FSWest),1,&
MPI_INTEGER,ngh(FCenter),1,COMM,stat,ierr)

! Front-NWest
CALL MPI_SENDRECV(ngh(NWest),1,MPI_INTEGER,ngh(BCenter),1,ngh(FNWest),1,&
MPI_INTEGER,ngh(FCenter),1,COMM,stat,ierr)

!PRINT*,'*******************************************************'
!PRINT*,'Processus: ',id,', my neighbor W is   : ',ngh(West)
!PRINT*,'Processus: ',id,', my neighbor E is   : ',ngh(East)
!PRINT*,'Processus: ',id,', my neighbor S is   : ',ngh(South)
!PRINT*,'Processus: ',id,', my neighbor N is   : ',ngh(North)
!PRINT*,'Processus: ',id,', my neighbor NE is  : ',ngh(NEast)
!PRINT*,'Processus: ',id,', my neighbor SE is  : ',ngh(SEast)
!PRINT*,'Processus: ',id,', my neighbor SW is  : ',ngh(SWest)
!PRINT*,'Processus: ',id,', my neighbor NW is  : ',ngh(NWest)
!PRINT*,'Processus: ',id,', my neighbor BC is  : ',ngh(BCenter)
!PRINT*,'Processus: ',id,', my neighbor FC is  : ',ngh(FCenter)
!PRINT*,'Processus: ',id,', my neighbor BN is  : ',ngh(BNorth)
!PRINT*,'Processus: ',id,', my neighbor BS is  : ',ngh(BSouth)
!PRINT*,'Processus: ',id,', my neighbor BE is  : ',ngh(BEast)
!PRINT*,'Processus: ',id,', my neighbor BW is  : ',ngh(BWest)
!PRINT*,'Processus: ',id,', my neighbor FN is  : ',ngh(FNorth)
!PRINT*,'Processus: ',id,', my neighbor FS is  : ',ngh(FSouth)
!PRINT*,'Processus: ',id,', my neighbor FE is  : ',ngh(FEast)
!PRINT*,'Processus: ',id,', my neighbor FW is  : ',ngh(FWest)
!PRINT*,'Processus: ',id,', my neighbor BNE is : ',ngh(BNEast)
!PRINT*,'Processus: ',id,', my neighbor BSE is : ',ngh(BSEast)
!PRINT*,'Processus: ',id,', my neighbor BSW is : ',ngh(BSWest)
!PRINT*,'Processus: ',id,', my neighbor BNW is : ',ngh(BNWest)
!PRINT*,'Processus: ',id,', my neighbor FNE is : ',ngh(FNEast)
!PRINT*,'Processus: ',id,', my neighbor FSE is : ',ngh(FSEast)
!PRINT*,'Processus: ',id,', my neighbor FSW is : ',ngh(FSWest)
!PRINT*,'Processus: ',id,', my neighbor FNW is : ',ngh(FNWest)
!PRINT*,'*******************************************************'

END SUBROUTINE COM_TOPOLOGY

!***********************************************************************
! Subroutine SET_FIELDS
! This subroutine generates the initial electromagnetic fields E and B and 
! initial current density J in the Yee lattice.
!
! INPUT:
! 
! - Bx: x-component of B
! - By: y-component of B
! - Bz: z-component of B
! - Ex: x-component of E
! - Ey: y-component of E
! - Ez: z-component of E
! - Jx: x-component of J
! - Jy: y-component of J
! - Jz: z-component of J
!
! OUTPUT: Bx,By,Bz,Ex,Ey,Ez,Jx,Jy,Jz
!***********************************************************************

SUBROUTINE SET_FIELDS(nd0,delta,de,B0,grad,sigma,Bx,By,Bz,Ex,Ey,Ez,Jx,Jy,Jz,&
                      xgp,ygp,zgp,xyeep,yyeep,zyeep)

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP) :: Bx,By,Bz
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP) :: Ex,Ey,Ez
DOUBLE PRECISION, DIMENSION(1:NXP,1:NYP,1:NZP) :: Jx,Jy,Jz

! Nodal grid in each domain
DOUBLE PRECISION, DIMENSION(1:NXP) :: xgp
DOUBLE PRECISION, DIMENSION(1:NYP) :: ygp
DOUBLE PRECISION, DIMENSION(1:NZP) :: zgp

! Yee grid in each domain
DOUBLE PRECISION, DIMENSION(1:NXP) :: xyeep
DOUBLE PRECISION, DIMENSION(1:NYP) :: yyeep
DOUBLE PRECISION, DIMENSION(1:NZP) :: zyeep

! Other plasma paramters
DOUBLE PRECISION :: delta,de,B0,J0,nd0,gamd,Az,sigma,grad
DOUBLE PRECISION :: y12,y14,y34,x14,x34,xm,ym

DOUBLE PRECISION :: bmag,r0,rmp,omega,rcyl,rot_shutoff_factor
DOUBLE PRECISION :: bxtemp,bytemp,bztemp

! Loop indexes
INTEGER :: ix,iy,iz,id

IF (INIT.EQ."RECONN") THEN
!***********************************************************************
! Initial plasma parameters according to the Harris equilibrium

! Bulk Lorentz factor drifting particles
gamd=1.0/sqrt(1.0-betad*betad)

! Upstream initial reconnecting magnetic field
B0=thde*me*c*c/(e*rhoc)

! Initial layer thickness
delta=2.0*thde*me*c*c/(betad*gamd*e*B0)

! Electron density of the drifting particles in the layer.
nd0=(thde*me*c*c)/(4.0*pi*e*e*betad*betad*gamd*delta*delta)

! Magnetization parameter
sigma=B0*B0/(4.0*pi*me*c*c*2.0*nd0*ratio_nb2nd)

! Electron skin depth.
de=sqrt(thde*me*c*c/(4.0*pi*nd0*e*e))

! Initial current in the +/-z-direction
J0=(c*thde*me*c*c)/(2.0*pi*e*betad*gamd*delta*delta)

! Electron maximum radiation reaction limit energy
grad=sqrt(3.0*e/(2.0*e**4.0/(me**2.0*c**4.0)*B0))

DO ix=1,NXP
DO iy=1,NYP
DO iz=1,NZP

! Initial magnetic field components

y12=(ymax-ymin)/2.0+ymin
y14=(ymax-ymin)/4.0+ymin
y34=3.0*(ymax-ymin)/4.0+ymin
x14=(xmax-xmin)/4.0+xmin
x34=3.0*(xmax-xmin)/4.0+xmin
xm=xmax-xmin
ym=ymax-ymin

IF (ygp(iy).LT.y12) THEN

Bx(ix,iy,iz)=-B0*tanh((yyeep(iy)-y14)/delta)*(1.0+perturb_amp*cos(2.0*pi*(xgp(ix)-&
x14)/xm)*(cos(2.0*pi*(yyeep(iy)-y14)/ym))**2.0)-2.0*(B0*delta*&
log(cosh(y14/delta))-delta*B0*log(cosh((yyeep(iy)-y14)/delta)))*&
perturb_amp*cos(2.0*pi*(xgp(ix)-x14)/xm)*2.0*pi/ym*cos(2.0*pi*(yyeep(iy)-&
y14)/ym)*sin(2.0*pi*(yyeep(iy)-y14)/ym)

By(ix,iy,iz)=(B0*delta*log(cosh(y14/delta))-delta*B0*log(cosh((ygp(iy)-&
y14)/delta)))*perturb_amp*2.0*pi/xm*sin(2.0*pi*(xyeep(ix)-&
x14)/xm)*(cos(2.0*pi*(ygp(iy)-y14)/ym))**2.0

ELSE

Bx(ix,iy,iz)=B0*tanh((yyeep(iy)-y34)/delta)*(1.0+perturb_amp*cos(2.0*pi*(xgp(ix)-&
x34)/xm)*(cos(2.0*pi*(yyeep(iy)-y34)/ym))**2.0)+2.0*(B0*delta*&
log(cosh(y14/delta))-delta*B0*log(cosh((yyeep(iy)-y34)/delta)))*&
perturb_amp*cos(2.0*pi*(xgp(ix)-x34)/xm)*2.0*pi/ym*cos(2.0*pi*(yyeep(iy)-&
y34)/ym)*sin(2.0*pi*(yyeep(iy)-y34)/ym)

By(ix,iy,iz)=-(B0*delta*log(cosh(y14/delta))-delta*B0*log(cosh((ygp(iy)-&
y34)/delta)))*perturb_amp*2.0*pi/xm*sin(2.0*pi*(xyeep(ix)-&
x34)/xm)*(cos(2.0*pi*(ygp(iy)-y34)/ym))**2.0

END IF

Bz(ix,iy,iz)=B0*guide_field

! Initial electric field components
Ex(ix,iy,iz)=0.0
Ey(ix,iy,iz)=0.0
Ez(ix,iy,iz)=0.0

! Initial current density components
Jx(ix,iy,iz)=0.0
Jy(ix,iy,iz)=0.0

IF (ygp(iy).LT.y12) THEN
! z-component potential vecteur
Az=(thde*me*c*c)/(e*betad*gamd)*log((cosh((ygp(iy)-y14)/delta))**(-2.0))
Jz(ix,iy,iz)=J0*exp(e*betad*gamd*Az/(thde*me*c*c))
ELSE
! z-component potential vecteur
Az=(thde*me*c*c)/(e*betad*gamd)*log((cosh((ygp(iy)-y34)/delta))**(-2.0))
Jz(ix,iy,iz)=-1.0*J0*exp(e*betad*gamd*Az/(thde*me*c*c))
END IF

ENDDO
ENDDO
ENDDO

ELSE IF (.FALSE.) THEN
! ELSE IF (INIT.EQ."UNIFORM") THEN
!***********************************************************************
! Initial plasma parameters according to a uniform field 
! The initial plasma particles are uniformly distributed and do not drift
! Optionally, drifting particles can be injected through a NOZZLE at velocity
! betad.

! Bulk Lorentz factor drifting particles
gamd=1.0/sqrt(1.0-betad*betad)

! Upstream initial reconnecting magnetic field
B0=me*c*c/(e*rhoc)

! Initial layer thickness
! Unused for INIT.NEQ."RECONN"
delta=2.0*thde*me*c*c/(betad*gamd*e*B0)

! Electron density of the drifting (i.e., injected) particles.
! nd0=(thde*me*c*c)/(4.0*pi*e*e*betad*betad*gamd*delta*delta)
nd0=B0*B0/(4.0*pi*me*c*c*2.0*sigma_inj)

! Magnetization parameter of initial (background) particles
sigma=B0*B0/(4.0*pi*me*c*c*2.0*nd0*ratio_nb2nd)

! Electron skin depth.
! Report the non-relativistic skin-depth for INIT.NEQ."RECONN"
! de=sqrt(thde*me*c*c/(4.0*pi*nd0*e*e))
de=sqrt(me*c*c/(4.0*pi*nd0*e*e))

! All currents are initially 0
J0=0d0 

! Electron maximum radiation reaction limit energy
grad=sqrt(3.0*e/(2.0*e**4.0/(me**2.0*c**4.0)*B0))

! Put nd0 to zero in the end if we're doing a vacuum run
IF (VACUUM_SIM) THEN
    nd0=0d0
END IF

DO ix=1,NXP
DO iy=1,NYP
DO iz=1,NZP

Bx(ix,iy,iz)=0.0
By(ix,iy,iz)=0.0
Bz(ix,iy,iz)=B0

Bx0(ix,iy,iz)=0.0
By0(ix,iy,iz)=0.0
Bz0(ix,iy,iz)=B0

! Set Ex
bmag=B0
bztemp=bmag
Ex(ix,iy,iz)=-xyeep(ix)*bztemp/rlc

rcyl=SQRT(xyeep(ix)*xyeep(ix) + ygp(iy)*ygp(iy))
rot_shutoff_factor=0.5*(1.0-TANH((rcyl-rnozzle)/delta_nozzle))
Ex(ix,iy,iz)=Ex(ix,iy,iz)*rot_shutoff_factor

! Set Ey
bmag=B0
bztemp=bmag
Ey(ix,iy,iz)=-yyeep(iy)*bztemp/rlc

rcyl=SQRT(xgp(ix)*xgp(ix) + yyeep(iy)*yyeep(iy))
rot_shutoff_factor=0.5*(1.0-TANH((rcyl-rnozzle)/delta_nozzle))
Ey(ix,iy,iz)=Ey(ix,iy,iz)*rot_shutoff_factor

! Set Ez
bxtemp=0.0
bytemp=0.0
Ez(ix,iy,iz)=(xgp(ix)*bxtemp + ygp(iy)*bytemp)/rlc

rcyl=SQRT(xgp(ix)*xgp(ix) + ygp(iy)*ygp(iy))
rot_shutoff_factor=0.5*(1.0-TANH((rcyl-rnozzle)/delta_nozzle))
Ez(ix,iy,iz)=Ez(ix,iy,iz)*rot_shutoff_factor

Ex0(ix,iy,iz)=Ex(ix,iy,iz)
Ey0(ix,iy,iz)=Ey(ix,iy,iz)
Ez0(ix,iy,iz)=Ez(ix,iy,iz)

! Zeroing out the initial field prevents a transient in the case that
! BOUND_FIELD_ZMIN is not set equal to NOZZLE.
! This is because the nontrivial, rotation-inducing E-field components
! are contained in the Ei0, which are only injected on the boundary if
! the NOZZLE BC is active.
! Comment this to be able to plot the initial E-field, though!
Ex(ix,iy,iz)=0.0
Ey(ix,iy,iz)=0.0
Ez(ix,iy,iz)=0.0

! Same as above but now just on grid nodes
bmag=B0
Bxg0(ix,iy,iz)=0.0
Byg0(ix,iy,iz)=0.0
Bzg0(ix,iy,iz)=bmag

Exg0(ix,iy,iz)=-xgp(ix)*Bzg0(ix,iy,iz)/rlc
Eyg0(ix,iy,iz)=-ygp(iy)*Bzg0(ix,iy,iz)/rlc
Ezg0(ix,iy,iz)=(xgp(ix)*Bxg0(ix,iy,iz) + ygp(iy)*Byg0(ix,iy,iz))/rlc

rcyl=SQRT(xgp(ix)*xgp(ix) + ygp(iy)*ygp(iy))
rot_shutoff_factor=0.5*(1.0-TANH((rcyl-rnozzle)/delta_nozzle))
Exg0(ix,iy,iz)=Exg0(ix,iy,iz)*rot_shutoff_factor
Eyg0(ix,iy,iz)=Eyg0(ix,iy,iz)*rot_shutoff_factor
Ezg0(ix,iy,iz)=Ezg0(ix,iy,iz)*rot_shutoff_factor

Jx(ix,iy,iz)=0.0
Jy(ix,iy,iz)=0.0
Jz(ix,iy,iz)=0.0

ENDDO
ENDDO
ENDDO

ELSE IF (.FALSE.) THEN
!ELSE IF (INIT.EQ."MONOPOLE") THEN
!***********************************************************************
! Initial plasma parameters according to a monopolar field 
! The initial plasma particles are uniformly distributed and do not drift
! Optionally, drifting particles can be injected through a NOZZLE at velocity
! betad.

! Bulk Lorentz factor drifting particles
gamd=1.0/sqrt(1.0-betad*betad)

! Upstream initial reconnecting magnetic field
B0=me*c*c/(e*rhoc)

! Initial layer thickness
! Unused for INIT.NEQ."RECONN"
delta=2.0*thde*me*c*c/(betad*gamd*e*B0)

! Electron density of the drifting (i.e., injected) particles.
! nd0=(thde*me*c*c)/(4.0*pi*e*e*betad*betad*gamd*delta*delta)
nd0=B0*B0/(4.0*pi*me*c*c*2.0*sigma_inj)

! Magnetization parameter of initial (background) particles
sigma=B0*B0/(4.0*pi*me*c*c*2.0*nd0*ratio_nb2nd)

! Electron skin depth.
! Report the non-relativistic skin-depth for INIT.NEQ."RECONN"
! de=sqrt(thde*me*c*c/(4.0*pi*nd0*e*e))
de=sqrt(me*c*c/(4.0*pi*nd0*e*e))

! All currents are initially 0
J0=0d0 

! Electron maximum radiation reaction limit energy
grad=sqrt(3.0*e/(2.0*e**4.0/(me**2.0*c**4.0)*B0))

! Put nd0 to zero in the end if we're doing a vacuum run
IF (VACUUM_SIM) THEN
    nd0=0d0
END IF

r0=(zmin-zmp)

DO ix=1,NXP
DO iy=1,NYP
DO iz=1,NZP

! Set Bx
rmp=SQRT(xgp(ix)*xgp(ix) + yyeep(iy)*yyeep(iy) + (zyeep(iz)-zmp)*(zyeep(iz)-zmp))
bmag=B0*r0*r0/(rmp*rmp)
Bx(ix,iy,iz)=bmag*(xgp(ix) - 0.0)/rmp

! Set By
rmp=SQRT(xyeep(ix)*xyeep(ix) + ygp(iy)*ygp(iy) + (zyeep(iz)-zmp)*(zyeep(iz)-zmp))
bmag=B0*r0*r0/(rmp*rmp)
By(ix,iy,iz)=bmag*(ygp(iy) - 0.0)/rmp

! Set Bz
rmp=SQRT(xyeep(ix)*xyeep(ix) + yyeep(iy)*yyeep(iy) + (zgp(iz)-zmp)*(zgp(iz)-zmp))
bmag=B0*r0*r0/(rmp*rmp)
Bz(ix,iy,iz)=bmag*(zgp(iz) - zmp)/rmp

Bx0(ix,iy,iz)=Bx(ix,iy,iz)
By0(ix,iy,iz)=By(ix,iy,iz)
Bz0(ix,iy,iz)=Bz(ix,iy,iz)

! Set Ex
rmp=SQRT(xyeep(ix)*xyeep(ix) + ygp(iy)*ygp(iy) + (zgp(iz)-zmp)*(zgp(iz)-zmp))
bmag=B0*r0*r0/(rmp*rmp)
bztemp=bmag*(zgp(iz) - zmp)/rmp
Ex(ix,iy,iz)=-xyeep(ix)*bztemp/rlc

rcyl=SQRT(xyeep(ix)*xyeep(ix) + ygp(iy)*ygp(iy))
rot_shutoff_factor=0.5*(1.0-TANH((rcyl-rnozzle)/delta_nozzle))
Ex(ix,iy,iz)=Ex(ix,iy,iz)*rot_shutoff_factor

! Set Ey
rmp=SQRT(xgp(ix)*xgp(ix) + yyeep(iy)*yyeep(iy) + (zgp(iz)-zmp)*(zgp(iz)-zmp))
bmag=B0*r0*r0/(rmp*rmp)
bztemp=bmag*(zgp(iz) - zmp)/rmp
Ey(ix,iy,iz)=-yyeep(iy)*bztemp/rlc

rcyl=SQRT(xgp(ix)*xgp(ix) + yyeep(iy)*yyeep(iy))
rot_shutoff_factor=0.5*(1.0-TANH((rcyl-rnozzle)/delta_nozzle))
Ey(ix,iy,iz)=Ey(ix,iy,iz)*rot_shutoff_factor

! Set Ez
rmp=SQRT(xgp(ix)*xgp(ix) + ygp(iy)*ygp(iy) + (zyeep(iz)-zmp)*(zyeep(iz)-zmp))
bmag=B0*r0*r0/(rmp*rmp)
bxtemp=bmag*(xgp(ix) - 0.0)/rmp
bytemp=bmag*(ygp(iy) - 0.0)/rmp
Ez(ix,iy,iz)=(xgp(ix)*bxtemp + ygp(iy)*bytemp)/rlc

rcyl=SQRT(xgp(ix)*xgp(ix) + ygp(iy)*ygp(iy))
rot_shutoff_factor=0.5*(1.0-TANH((rcyl-rnozzle)/delta_nozzle))
Ez(ix,iy,iz)=Ez(ix,iy,iz)*rot_shutoff_factor

Ex0(ix,iy,iz)=Ex(ix,iy,iz)
Ey0(ix,iy,iz)=Ey(ix,iy,iz)
Ez0(ix,iy,iz)=Ez(ix,iy,iz)

! Zeroing out the initial field prevents a transient in the case that
! BOUND_FIELD_ZMIN is not set equal to NOZZLE.
! This is because the nontrivial, rotation-inducing E-field components
! are contained in the Ei0, which are only injected on the boundary if
! the NOZZLE BC is active.
Ex(ix,iy,iz)=0.0
Ey(ix,iy,iz)=0.0
Ez(ix,iy,iz)=0.0

! Same as above but now just on grid nodes
rmp=SQRT(xgp(ix)*xgp(ix) + ygp(iy)*ygp(iy) + (zgp(iz)-zmp)*(zgp(iz)-zmp))
bmag=B0*r0*r0/(rmp*rmp)
Bxg0(ix,iy,iz)=bmag*(xgp(ix) - 0.0)/rmp
Byg0(ix,iy,iz)=bmag*(ygp(iy) - 0.0)/rmp
Bzg0(ix,iy,iz)=bmag*(zgp(iz) - zmp)/rmp

Exg0(ix,iy,iz)=-xgp(ix)*Bzg0(ix,iy,iz)/rlc
Eyg0(ix,iy,iz)=-ygp(iy)*Bzg0(ix,iy,iz)/rlc
Ezg0(ix,iy,iz)=(xgp(ix)*Bxg0(ix,iy,iz) + ygp(iy)*Byg0(ix,iy,iz))/rlc

rcyl=SQRT(xgp(ix)*xgp(ix) + ygp(iy)*ygp(iy))
rot_shutoff_factor=0.5*(1.0-TANH((rcyl-rnozzle)/delta_nozzle))
Exg0(ix,iy,iz)=Exg0(ix,iy,iz)*rot_shutoff_factor
Eyg0(ix,iy,iz)=Eyg0(ix,iy,iz)*rot_shutoff_factor
Ezg0(ix,iy,iz)=Ezg0(ix,iy,iz)*rot_shutoff_factor

Jx(ix,iy,iz)=0.0
Jy(ix,iy,iz)=0.0
Jz(ix,iy,iz)=0.0

ENDDO
ENDDO
ENDDO

ELSE IF ((INIT.EQ."UNIFORM").OR.(INIT.EQ."MONOPOLE")) THEN
!***********************************************************************
! Initial plasma parameters according to a monopolar or uniform initial field 
! The initial plasma particles are uniformly distributed and do not drift
! Optionally, drifting particles can be injected through a NOZZLE at velocity
! betad.

! Bulk Lorentz factor drifting particles
gamd=1.0/sqrt(1.0-betad*betad)

! Upstream initial reconnecting magnetic field
B0=me*c*c/(e*rhoc)

! Initial layer thickness
! Unused for INIT.NEQ."RECONN"
delta=2.0*thde*me*c*c/(betad*gamd*e*B0)

! Electron density of the drifting (i.e., injected) particles.
! nd0=(thde*me*c*c)/(4.0*pi*e*e*betad*betad*gamd*delta*delta)
nd0=B0*B0/(4.0*pi*me*c*c*2.0*sigma_inj)

! Magnetization parameter of initial (background) particles
sigma=B0*B0/(4.0*pi*me*c*c*2.0*nd0*ratio_nb2nd)

! Electron skin depth.
! Report the non-relativistic skin-depth for INIT.NEQ."RECONN"
! de=sqrt(thde*me*c*c/(4.0*pi*nd0*e*e))
de=sqrt(me*c*c/(4.0*pi*nd0*e*e))

! All currents are initially 0
J0=0d0 

! Electron maximum radiation reaction limit energy
grad=sqrt(3.0*e/(2.0*e**4.0/(me**2.0*c**4.0)*B0))

! Put nd0 to zero in the end if we're doing a vacuum run
IF (VACUUM_SIM) THEN
    nd0=0d0
END IF

DO ix=1,NXP
DO iy=1,NYP
DO iz=1,NZP

! Set Bx
CALL INITIAL_BCOMPONENTS(xgp(ix),yyeep(iy),zyeep(iz),bxtemp,bytemp,bztemp)
Bx(ix,iy,iz) = bxtemp
! rmp=SQRT(xgp(ix)*xgp(ix) + yyeep(iy)*yyeep(iy) + (zyeep(iz)-zmp)*(zyeep(iz)-zmp))
! bmag=B0*r0*r0/(rmp*rmp)
! Bx(ix,iy,iz)=bmag*(xgp(ix) - 0.0)/rmp

! Set By
CALL INITIAL_BCOMPONENTS(xyeep(ix),ygp(iy),zyeep(iz),bxtemp,bytemp,bztemp)
By(ix,iy,iz) = bytemp
! rmp=SQRT(xyeep(ix)*xyeep(ix) + ygp(iy)*ygp(iy) + (zyeep(iz)-zmp)*(zyeep(iz)-zmp))
! bmag=B0*r0*r0/(rmp*rmp)
! By(ix,iy,iz)=bmag*(ygp(iy) - 0.0)/rmp

! Set Bz
CALL INITIAL_BCOMPONENTS(xyeep(ix),yyeep(iy),zgp(iz),bxtemp,bytemp,bztemp)
Bz(ix,iy,iz) = bztemp
! rmp=SQRT(xyeep(ix)*xyeep(ix) + yyeep(iy)*yyeep(iy) + (zgp(iz)-zmp)*(zgp(iz)-zmp))
! bmag=B0*r0*r0/(rmp*rmp)
! Bz(ix,iy,iz)=bmag*(zgp(iz) - zmp)/rmp

Bx0(ix,iy,iz) = Bx(ix,iy,iz)
By0(ix,iy,iz) = By(ix,iy,iz)
Bz0(ix,iy,iz) = Bz(ix,iy,iz)

! Set Ex
CALL INITIAL_BCOMPONENTS(xyeep(ix),ygp(iy),zgp(iz),bxtemp,bytemp,bztemp)
! rmp=SQRT(xyeep(ix)*xyeep(ix) + ygp(iy)*ygp(iy) + (zgp(iz)-zmp)*(zgp(iz)-zmp))
! bmag=B0*r0*r0/(rmp*rmp)
! bztemp=bmag*(zgp(iz) - zmp)/rmp
Ex(ix,iy,iz) = -xyeep(ix)*bztemp/rlc

! rcyl=SQRT(xyeep(ix)*xyeep(ix) + ygp(iy)*ygp(iy))
! rot_shutoff_factor=0.5*(1.0-TANH((rcyl-rnozzle)/delta_nozzle))
rot_shutoff_factor = ROTATION_PROFILE(xyeep(ix),ygp(iy))
Ex(ix,iy,iz) = Ex(ix,iy,iz)*rot_shutoff_factor

! Set Ey
CALL INITIAL_BCOMPONENTS(xgp(ix),yyeep(iy),zgp(iz),bxtemp,bytemp,bztemp)
! rmp=SQRT(xgp(ix)*xgp(ix) + yyeep(iy)*yyeep(iy) + (zgp(iz)-zmp)*(zgp(iz)-zmp))
! bmag=B0*r0*r0/(rmp*rmp)
! bztemp=bmag*(zgp(iz) - zmp)/rmp
Ey(ix,iy,iz) = -yyeep(iy)*bztemp/rlc

! rcyl=SQRT(xgp(ix)*xgp(ix) + yyeep(iy)*yyeep(iy))
! rot_shutoff_factor=0.5*(1.0-TANH((rcyl-rnozzle)/delta_nozzle))
rot_shutoff_factor = ROTATION_PROFILE(xgp(ix),yyeep(iy))
Ey(ix,iy,iz) = Ey(ix,iy,iz)*rot_shutoff_factor

! Set Ez
CALL INITIAL_BCOMPONENTS(xgp(ix),ygp(iy),zyeep(iz),bxtemp,bytemp,bztemp)
! rmp=SQRT(xgp(ix)*xgp(ix) + ygp(iy)*ygp(iy) + (zyeep(iz)-zmp)*(zyeep(iz)-zmp))
! bmag=B0*r0*r0/(rmp*rmp)
! bxtemp=bmag*(xgp(ix) - 0.0)/rmp
! bytemp=bmag*(ygp(iy) - 0.0)/rmp
Ez(ix,iy,iz) = (xgp(ix)*bxtemp + ygp(iy)*bytemp)/rlc

! rcyl=SQRT(xgp(ix)*xgp(ix) + ygp(iy)*ygp(iy))
! rot_shutoff_factor=0.5*(1.0-TANH((rcyl-rnozzle)/delta_nozzle))
rot_shutoff_factor = ROTATION_PROFILE(xgp(ix),ygp(iy))
Ez(ix,iy,iz) = Ez(ix,iy,iz)*rot_shutoff_factor

Ex0(ix,iy,iz) = Ex(ix,iy,iz)
Ey0(ix,iy,iz) = Ey(ix,iy,iz)
Ez0(ix,iy,iz) = Ez(ix,iy,iz)

! Zeroing out the initial field prevents a transient in the case that
! BOUND_FIELD_ZMIN is not set equal to NOZZLE.
! This is because the nontrivial, rotation-inducing E-field components
! are contained in the Ei0, which are only injected on the boundary if
! the NOZZLE BC is active.
Ex(ix,iy,iz)=0.0
Ey(ix,iy,iz)=0.0
Ez(ix,iy,iz)=0.0

! Same as above but now just on grid nodes
CALL INITIAL_BCOMPONENTS(xgp(ix),ygp(iy),zgp(iz),bxtemp,bytemp,bztemp)
! rmp=SQRT(xgp(ix)*xgp(ix) + ygp(iy)*ygp(iy) + (zgp(iz)-zmp)*(zgp(iz)-zmp))
! bmag=B0*r0*r0/(rmp*rmp)
! Bxg0(ix,iy,iz)=bmag*(xgp(ix) - 0.0)/rmp
! Byg0(ix,iy,iz)=bmag*(ygp(iy) - 0.0)/rmp
! Bzg0(ix,iy,iz)=bmag*(zgp(iz) - zmp)/rmp
Bxg0(ix,iy,iz) = bxtemp
Byg0(ix,iy,iz) = bytemp
Bzg0(ix,iy,iz) = bztemp

Exg0(ix,iy,iz) = -xgp(ix)*Bzg0(ix,iy,iz)/rlc
Eyg0(ix,iy,iz) = -ygp(iy)*Bzg0(ix,iy,iz)/rlc
Ezg0(ix,iy,iz) = (xgp(ix)*Bxg0(ix,iy,iz) + ygp(iy)*Byg0(ix,iy,iz))/rlc

! rcyl=SQRT(xgp(ix)*xgp(ix) + ygp(iy)*ygp(iy))
! rot_shutoff_factor=0.5*(1.0-TANH((rcyl-rnozzle)/delta_nozzle))
rot_shutoff_factor = ROTATION_PROFILE(xgp(ix),ygp(iy))
Exg0(ix,iy,iz) = Exg0(ix,iy,iz)*rot_shutoff_factor
Eyg0(ix,iy,iz) = Eyg0(ix,iy,iz)*rot_shutoff_factor
Ezg0(ix,iy,iz) = Ezg0(ix,iy,iz)*rot_shutoff_factor

Jx(ix,iy,iz) = 0.0
Jy(ix,iy,iz) = 0.0
Jz(ix,iy,iz) = 0.0

ENDDO
ENDDO
ENDDO

ELSE

PRINT *, "Unrecognized initial condition: ", INIT

END IF

END SUBROUTINE SET_FIELDS

!***********************************************************************
! Subroutine INITIAL_BCOMPONENTS
! Note that the output of this function depends on whether the LAST ARGUMENT
! is supplied! Read carefully.
!  
! INPUT:
! - x,y,z: spatial position
! 
! OUTPUT:
! - bhatx,bhaty,bhatz,bmag: components of the magnetic field initial condition
!     at position (x,y,z)
!     The magnetic field vector is (bhatx,bhaty,bhatz)*bmag
!     The magnetic field magnitude is bmag
!     The magnetic field unit vector is (bhatx,bhaty,bhatz)
! 
!   OR, if bmag IS NOT SUPPLIED,
!
! - bx,by,bz: components of the magnetic field initial condition at (x,y,z)
!     The magnetic field vector is (bhatx,bhaty,bhatz)
!***********************************************************************
SUBROUTINE INITIAL_BCOMPONENTS(x,y,z,bhatx,bhaty,bhatz,bmag)

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN)            :: x,y,z
DOUBLE PRECISION, INTENT(OUT)           :: bhatx,bhaty,bhatz
DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: bmag

DOUBLE PRECISION :: r0,rmp,bmagtemp,B0

B0=me*c*c/(e*rhoc)

IF (INIT.EQ."UNIFORM") THEN

bmagtemp = B0
bhatx    = 0d0
bhaty    = 0d0
bhatz    = 1d0

ELSE IF (INIT.EQ."MONOPOLE") THEN

! Distance to the monopole at (x=0,y=0,zmin)
r0 = (zmin-zmp)
! Distance to the monopole at (x,y,z)
rmp = SQRT((x - 0d0)*(x - 0d0) + (y - 0d0)*(y - 0d0) + (z - zmp)*(z - zmp))

bmagtemp = B0*r0*r0/(rmp*rmp)
bhatx    = (x - 0d0)/rmp
bhaty    = (y - 0d0)/rmp
bhatz    = (z - zmp)/rmp

ELSE

PRINT *, "Unrecognized initial condition: ", INIT

ENDIF

IF (PRESENT(bmag)) THEN
    bmag = bmagtemp
ELSE
    bhatx = bhatx*bmagtemp
    bhaty = bhaty*bmagtemp
    bhatz = bhatz*bmagtemp
END IF

END SUBROUTINE INITIAL_BCOMPONENTS

!***********************************************************************
! Subroutine ROTATION_PROFILE
! For the NOZZLE boundary condition applied at z=zmin, the rotation cuts
! off smoothly, such that the angular frequency is just
! 
!     (c/rcyl) * ROTATION_PROFILE(x,y,z=zmin)
!  
! INPUT:
! - x,y: spatial position
! 
! OUTPUT:
! - nozzle rotation profile applied at (x,y,z=zmin)
!***********************************************************************
FUNCTION ROTATION_PROFILE(x,y)

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN)  :: x,y
DOUBLE PRECISION              :: ROTATION_PROFILE

DOUBLE PRECISION :: rcyl
DOUBLE PRECISION :: clamp_threshold = 4.0

rcyl = SQRT( (x - 0d0)*(x - 0d0) + (y - 0d0)*(y - 0d0) )

! Clamp to 0 or 1 in case delta_nozzle is 0
! (avoid div-by-zero)
IF (rcyl - rnozzle < -clamp_threshold * delta_nozzle) THEN
    ROTATION_PROFILE=1d0
ELSE IF (rcyl - rnozzle > clamp_threshold * delta_nozzle) THEN
    ! This has the nice side-effect of shutting off macroparticle injection
    ! once rcyl > clamp_threshold * delta_nozzle + rnozzle
    ROTATION_PROFILE=0d0
ELSE
    ROTATION_PROFILE=0.5*(1d0 - TANH( (rcyl - rnozzle)/delta_nozzle ))
END IF

END FUNCTION ROTATION_PROFILE


!***********************************************************************
! Subroutine WEIGHT
! This subroutine weight each super-particles
!
! INPUT:
! 
! - pcl0: Initial particle distribution function
! - delta: initial layer thickness
! - NPP: Number of initial particles per process
!
! OUTPUT: wght=pcl0(7,:)
!***********************************************************************

SUBROUTINE WEIGHT(pcl0,delta,NPP)

IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(1:7,1:NPP) :: pcl0
DOUBLE PRECISION                       :: delta,y,y12,y14,y34
INTEGER*8                              :: ip,NPP

!***********************************************************************

y12=(ymax-ymin)/2.0+ymin
y14=(ymax-ymin)/4.0+ymin
y34=3.0*(ymax-ymin)/4.0+ymin

! From the Harris equilibirum

DO ip=1,NPP

y=pcl0(2,ip)

IF (y.LT.y12) THEN
pcl0(7,ip)=(cosh((y-y14)/delta))**(-2.0)
ELSE
pcl0(7,ip)=(cosh((y-y34)/delta))**(-2.0)
END IF

ENDDO

END SUBROUTINE WEIGHT

!***********************************************************************
! Subroutine SET_TAG
! This subroutine assigns a tag to each super-particles
!
! INPUT:
! 
! - tag: tag of the super-particle
! - id: process rank
! - NPP: Number of initial particles per process
!
! OUTPUT: tag
!***********************************************************************

SUBROUTINE SET_TAG(tag,id,NPP)

IMPLICIT NONE

INTEGER*8, DIMENSION(1:NPP) :: tag
INTEGER*8 :: ip,NPP
INTEGER   :: id

!***********************************************************************

tag=0

DO ip=1,NPP
tag(ip)=id*NPP+ip ! tag must be > 0 so -tag indicates a tracked particle
ENDDO

END SUBROUTINE SET_TAG

!***********************************************************************
! Subroutine SET_MAXWELLIAN
! This subroutine generates the initial distribution function of the
! particles with a Maxwellian distribution in energy.
!
! INPUT: 
! - pcl0: Initial distribution function of the particles
! - theta: temperature of plasma in mc^2/k units
! - xminp,yminp,zminp: lower spatial boundary for each domain
! - NPP: Number of initial particles per process
!
! OUTPUT: Particle distribution function
!***********************************************************************

SUBROUTINE SET_MAXWELLIAN(pcl0,theta,xminp,yminp,zminp,NPP)

IMPLICIT NONE

! Temperature in mc^2/k units
DOUBLE PRECISION :: theta,umint,umaxt

! Cosine of the polar angle between the z-axis and the velocity
DOUBLE PRECISION, PARAMETER           :: cmin=-1.0,cmax=1.0

! Azimuth between x-axis and the projected velocity in the xy-plane
DOUBLE PRECISION, PARAMETER           :: pmin=0.0,pmax=2.0*pi

! Number of elements for integration distribution function
INTEGER, PARAMETER :: ND=500

! Number of initial particles per process
INTEGER*8 :: NPP

! Lower spatial boundary for each domain
DOUBLE PRECISION :: xminp,yminp,zminp

DOUBLE PRECISION, DIMENSION(1:ND)           :: ud,gFu
DOUBLE PRECISION, DIMENSION(1:PPC)          :: x0c,y0c,z0c,cth0,phi0
DOUBLE PRECISION, DIMENSION(1:PPC)          :: u0,ux0c,uy0c,uz0c
DOUBLE PRECISION, DIMENSION(1:7,1:NPP)      :: pcl0

! Loop indexes
INTEGER   :: ic,id,ix,iy,iz
INTEGER*8 :: ip

!***********************************************************************

CALL INIT_RANDOM_SEED()

ip=1

pcl0=0.0
ud=0.0

! Cumulative distribution function
gFu=0.0

umint=1d-2*sqrt(theta)
umaxt=1d2*sqrt(theta)

DO id=1,ND
ud(id)=1d1**((id-1)*1d0/((ND-1)*1d0)*(log10(umaxt)-log10(umint))+log10(umint))
ENDDO

! Numerical calculation of the cumulative distribution function
CALL DIST_MAXWELL(ud,gFu,theta,ND)

DO ix=1,NCXP
DO iy=1,NCYP
DO iz=1,NCZP

!***********************************************************************
! Definition of the initial position x0 as a uniform random value

x0c=0.0

CALL RANDOM_NUMBER(x0c)

x0c=xminp+(ix-1)*dx+x0c*dx

! Definition of the initial position y0 as a uniform random value

y0c=0.0

CALL RANDOM_NUMBER(y0c)

y0c=yminp+(iy-1)*dy+y0c*dy

! Definition of the initial position z0 as a uniform random value

z0c=0.0

CALL RANDOM_NUMBER(z0c)

z0c=zminp+(iz-1)*dz+z0c*dz

!***********************************************************************
! Definition of the initial particle 4-velocity

u0=0.0

CALL GEN_U(u0,ud,gFu,ND)

!***********************************************************************
! Definition of the initial phi0 as a uniform random value

phi0=0.0

CALL RANDOM_NUMBER(phi0)

 phi0=phi0*(pmax-pmin)+pmin

! Definition of the initial cth0=cos(theta0) as a uniform random value

 cth0=0.0

CALL RANDOM_NUMBER(cth0)

 cth0=cth0*(cmax-cmin)+cmin

! Initial 4-velocity components
ux0c=u0*sqrt(1.0-cth0*cth0)*cos(phi0)
uy0c=u0*sqrt(1.0-cth0*cth0)*sin(phi0)
uz0c=u0*cth0

DO ic=1,PPC
pcl0(1,ip)=x0c(ic)
pcl0(2,ip)=y0c(ic)
pcl0(3,ip)=z0c(ic)
pcl0(4,ip)=ux0c(ic)
pcl0(5,ip)=uy0c(ic)
pcl0(6,ip)=uz0c(ic)
ip=ip+1
ENDDO

ENDDO
ENDDO
ENDDO

END SUBROUTINE SET_MAXWELLIAN

!***********************************************************************
! Subroutine DIST_MAXWELL
! This subroutine generates a random distribution of particle Lorentz 
! factor gamma following a relativistic Maxwellian distribution.
!
! INPUT: 
! - ud: reference 4-velocity array
! - gFu: pre-calculated integrated distribution function
! - theta: temperature of plasma in mc^2/k units
! - ND: number of elements in the array ud and gFu
!
! OUTPUT: Integrated distribution function gFu
! :TODO: parallelize the outer loops?
!***********************************************************************

SUBROUTINE DIST_MAXWELL(ud,gFu,theta,ND)

IMPLICIT NONE

! Input parameter
DOUBLE PRECISION :: theta
INTEGER          :: ND
DOUBLE PRECISION, DIMENSION(1:ND) :: ud,gFu,fu,udp,fup
DOUBLE PRECISION :: norm,hu,upmax,hup

! Loop indexes
INTEGER :: id1,id2

!***********************************************************************

! Relativistic Maxwellian distribution
fu=ud*ud*exp(-(sqrt(1.0+ud*ud)-1.0)/theta)*ud

norm=0.0

DO id1=1,ND-1
hu=(log(ud(id1+1))-log(ud(id1)))
norm=norm+hu/2.0*(fu(id1+1)+fu(id1))
ENDDO

DO id1=1,ND

upmax=ud(id1)

  DO id2=1,ND
  udp(id2)=1d1**((id2-1)*1d0/((ND-1)*1d0)*(log10(upmax)-&
           log10(minval(ud)))+log10(minval(ud)))
  ENDDO

  ! Relativistic Maxwellian distribution
  fup=udp*udp*exp(-(sqrt(1.0+udp*udp)-1.0)/theta)*udp

  DO id2=1,ND-1
  hup=(log(udp(id2+1))-log(udp(id2)))
  gFu(id1)=gFu(id1)+hup/2.0*(fup(id2+1)+fup(id2))
  ENDDO

ENDDO

gFu=gFu/norm

END SUBROUTINE DIST_MAXWELL

!***********************************************************************
! Subroutine GEN_U
! This subroutine generates a random distribution of particle 4-velocity
!
! INPUT: 
! - u0: 4-velocity array
! - ud: reference 4-velocity array
! - uFg: pre-calculated integrated distribution function
! - theta: temperature of plasma in mc^2/k units
! - ND: number of elements in the array gd and gFg
!
! OUTPUT: distribution of u0
!***********************************************************************

SUBROUTINE GEN_U(u0,ud,gFu,ND)

IMPLICIT NONE

INTEGER :: ND
DOUBLE PRECISION, DIMENSION(1:ND)  :: ud,gFu
DOUBLE PRECISION, DIMENSION(1:PPC) :: u0,Ru
DOUBLE PRECISION                   :: gFu1,gFu2,u1,u2
INTEGER, DIMENSION(1)              :: minu

! Loop indexes
INTEGER :: ic,iu

!***********************************************************************

Ru=0.0

CALL RANDOM_NUMBER(Ru)
Ru=Ru*0.9999

DO ic=1,PPC

minu=minloc(abs(gFu-Ru(ic)))
iu=minu(1)

IF (iu.EQ.ND) THEN
iu=iu-1
END IF

gFu1=gFu(iu)
gFu2=gFu(iu+1)

u1=ud(iu)
u2=ud(iu+1)

u0(ic)=exp((Ru(ic)*(log(u2)-log(u1))-log(u2)*gFu1+log(u1)*gFu2)/(gFu2-gFu1))

ENDDO

END SUBROUTINE GEN_U

!***********************************************************************
! Subroutine SET_DRIFT_MAXWELLIAN
! This subroutine generates the initial distribution function of the
! particles with a drifting Maxwellian distribution in energy.
!
! INPUT: 
! - q: electric charge sign
! - pcl0: Initial distribution function of the particles
! - theta: Temperature of plasma in mc^2/k units
! - up,gFp,ps,gFs,ND: Particle distribution generator quantities
! - xminp,yminp: Lower spatial boundary for each domain
! - NPP: Number of initial particles per process
!
! OUTPUT: Particle distribution function
!***********************************************************************

SUBROUTINE SET_DRIFT_MAXWELLIAN(q,pcl0,theta,up,gFp,ps,gFs,xminp,yminp,zminp,ND,NPP)

IMPLICIT NONE

! Temperature in mc^2/k units
DOUBLE PRECISION :: theta,q

! Angle defined in the plane perpendicular to the direction of the bulk motion
DOUBLE PRECISION, PARAMETER           :: pmin=0.0,pmax=2.0*pi

! Number of elements in the pre-calculated cumulative distribution functions
INTEGER :: ND

! Number of initial particles per process
INTEGER*8 :: NPP

! Lower spatial boundary for each domain
DOUBLE PRECISION :: xminp,yminp,zminp,gamd,y12

DOUBLE PRECISION, DIMENSION(1:ND)           :: up,gFp,ps
DOUBLE PRECISION, DIMENSION(1:ND,1:ND)      :: gFs
DOUBLE PRECISION, DIMENSION(1:PPC)          :: x0c,y0c,z0c,phi0
DOUBLE PRECISION, DIMENSION(1:PPC)          :: up0,ps0,uperp0,ux0c,uy0c,uz0c
DOUBLE PRECISION, DIMENSION(1:7,1:NPP)      :: pcl0

! Loop indexes
INTEGER   :: ic,ix,iy,iz
INTEGER*8 :: ip

!***********************************************************************

! Drift Lorentz factor
gamd=1.0/sqrt(1.0-betad*betad)

! Mid-box
y12=(ymax-ymin)/2.0+ymin

! Initialization of the seed for the generation of random variables
CALL INIT_RANDOM_SEED()

ip=1

pcl0=0.0

DO ix=1,NCXP
DO iy=1,NCYP
DO iz=1,NCZP

!***********************************************************************
! Definition of the initial position x0 as a uniform random value

x0c=0.0

CALL RANDOM_NUMBER(x0c)

x0c=xminp+(ix-1)*dx+x0c*dx

! Definition of the initial position y0 as a uniform random value

y0c=0.0

CALL RANDOM_NUMBER(y0c)

y0c=yminp+(iy-1)*dy+y0c*dy

! Definition of the initial position z0 as a uniform random value

z0c=0.0

CALL RANDOM_NUMBER(z0c)

z0c=zminp+(iz-1)*dz+z0c*dz

!***********************************************************************
! Definition of the initial momentum parallel to the drift velocity

up0=0.0

CALL GEN_UP(up0,up,gFp,ND,PPC)

!***********************************************************************
! Definition of the initial perpendicular momentum

ps0=0.0

CALL GEN_PS(ps0,ps,up0,up,gFs,ND,PPC)

uperp0=ps0*sqrt(1.0+up0*up0)
  
!***********************************************************************
! Definition random angle in the xy-plane

phi0=0.0

CALL RANDOM_NUMBER(phi0)

 phi0=phi0*(pmax-pmin)+pmin

! Initial 4-velocity
ux0c=uperp0*cos(phi0)
uy0c=uperp0*sin(phi0)
uz0c=up0

DO ic=1,PPC
pcl0(1,ip)=x0c(ic)
pcl0(2,ip)=y0c(ic)
pcl0(3,ip)=z0c(ic)
pcl0(4,ip)=ux0c(ic)
pcl0(5,ip)=uy0c(ic)

  IF (y0c(ic).LT.y12) THEN
  pcl0(6,ip)=q*uz0c(ic)
  ELSE
  pcl0(6,ip)=-q*uz0c(ic)
  END IF

ip=ip+1
ENDDO

ENDDO
ENDDO
ENDDO

END SUBROUTINE SET_DRIFT_MAXWELLIAN

!***********************************************************************
! Subroutine INIT_DRIFT_MAXWELLIAN
! This subroutine computes the cumulative distribution functions for
! generating a drifting relativistic maxwellian distribution of particles.
! The distributions are taken from Swisdak (2013), PoP.
!
! INPUT:
! - up: parallel momentum array
! - gFp: cumulative distribution function for the parallel momentum
! - ps=uperp/sqrt(1+up*up)
! - gFs: conditional distribution function for ps knowing up
! - theta: temperature of plasma in mc^2/k units
! - ND: number of elements for all arrays
! - COMM: mpi comm
!
! OUTPUT: Integrated distribution function gFp and gFs
!***********************************************************************

SUBROUTINE INIT_DRIFT_MAXWELLIAN(up,gFp,ps,gFs,theta,ND, COMM)

IMPLICIT NONE

INCLUDE 'mpif.h'

! Input parameter
DOUBLE PRECISION :: theta
INTEGER          :: ND, COMM
DOUBLE PRECISION, DIMENSION(1:ND)      :: up,gFp,fg1,up2,fg2,g1,g2
DOUBLE PRECISION, DIMENSION(1:ND)      :: ps,gs1,gs2,ps2,fs1,fs2
DOUBLE PRECISION, DIMENSION(1:ND,1:ND) :: gFs, gFs2
DOUBLE PRECISION :: norm1,norm2,hp,hp2,gamd,pu,umint,hs,hs2

INTEGER :: rank, numRanks, mpiErr
! Loop indexes
INTEGER :: i1,i2,i3

!***********************************************************************

CALL MPI_COMM_RANK(COMM, rank, mpiErr)
CALL MPI_COMM_SIZE(COMM, numRanks, mpiErr)

! Drift Lorentz factor
gamd=1.0/sqrt(1.0-betad*betad)

! Drift 4-velocity
pu=gamd*betad

g1=sqrt(1.0+up*up)

! f(p||) see Eq.(12) Swisdak (2013)
fg1=(1.0+gamd*g1/theta)*exp(-(up-pu)**2.0/(g1*gamd+up*pu+1d0)/theta)

norm1=0.0

DO i1=1,ND-1
hp=up(i1+1)-up(i1)
norm1=norm1+hp/2.0*(fg1(i1+1)+fg1(i1))
ENDDO

umint=minval(up)

!***********************************************************************
! CALCULATION OF THE PARALLEL MOMENTUM CUMULATIVE DISTRIBUTION gFp
!***********************************************************************

gFp=0.0

DO i1=1,ND

  DO i2=1,ND
  up2(i2)=(i2-1)*1d0/((ND-1)*1d0)*(up(i1)-umint)+umint
  ENDDO

  g2=sqrt(1.0+up2*up2)
  
  ! f(p||) see Eq.(12) Swisdak (2013)
  fg2=(1.0+gamd*g2/theta)*exp(-(up2-pu)**2.0/(g2*gamd+up2*pu+1d0)/theta)
  
  DO i2=1,ND-1
  hp2=up2(i2+1)-up2(i2)
  gFp(i1)=gFp(i1)+hp2/2.0*(fg2(i2+1)+fg2(i2))
  ENDDO

ENDDO

gFp=gFp/norm1

!***********************************************************************
! CALCULATION OF THE CONDITIONAL PERPENDICULAR MOMENTUM CUMULATIVE
! DISTRIBUTION gFs
!***********************************************************************

gs1=sqrt(1.0+ps*ps)
gFs2=0.0

! GRW: this triple loop is time-consuming; it decreases parallel 
! efficiency when run on many processors, but I believe each outer-loop
! is independent, so farm it out.
! N.B. For ND=1000, on verus, this subroutine took about 97 s, in serial.
!   After parallelizing, it took about 3 s on 48 procs.

DO i1=1,ND
  
  IF (MOD(i1-1, numRanks) == rank) THEN !{
  
  ! f(ps|up) see Eq.(15) Swisdak (2013)
  fs1=ps*exp(-((up(i1)-pu)**2.0+g1(i1)*g1(i1)*gamd*gamd*ps*ps)/&
              (g1(i1)*gamd*gs1+up(i1)*pu+1.0)/theta)
  
  norm2=0.0
  
  DO i2=1,ND-1
    hs=ps(i2+1)-ps(i2)
    norm2=norm2+hs/2.0*(fs1(i2+1)+fs1(i2))
  ENDDO
  
  !***********************************************************************

  DO i2=1,ND

     DO i3=1,ND
     ps2(i3)=(i3-1)*1d0/((ND-1)*1d0)*(ps(i2)-0.0)+0.0
     ENDDO
     
     gs2=sqrt(1.0+ps2*ps2)
     
     ! f(ps|p||) see Eq.(15) Swisdak (2013)
     fs2=ps2*exp(-((up(i1)-pu)**2.0+g1(i1)*g1(i1)*gamd*gamd*ps2*ps2)/&
                  (g1(i1)*gamd*gs2+up(i1)*pu+1.0)/theta)
     
     DO i3=1,ND-1
     hs2=ps2(i3+1)-ps2(i3)
     gFs2(i2,i1)=gFs2(i2,i1)+hs2/2.0*(fs2(i3+1)+fs2(i3))
     ENDDO

     gFs2(i2,i1)=gFs2(i2,i1)/norm2     

  ENDDO

  ENDIF !}
     
ENDDO

IF (numRanks > 1) THEN
  CALL MPI_ALLREDUCE(gFs2, gFs, SIZE(gFs), MPI_DOUBLE_PRECISION, MPI_SUM,&
    COMM, mpiErr)
ELSE
  gFs = gFs2
ENDIF

END SUBROUTINE INIT_DRIFT_MAXWELLIAN

!***********************************************************************
! Subroutine GEN_UP
! This subroutine generates a random distribution of particle momenta 
! parallel to the drift velocity.
!
! INPUT: 
! - up0: Particle momentum parallel to the drift velocity
! - up: reference momentum array
! - gFp: pre-calculated integrated distribution function
! - ND: number of elements in the array ug and gFg
!
! OUTPUT: distribution of up0
!***********************************************************************

SUBROUTINE GEN_UP(up0,up,gFp,ND,NINJ)

IMPLICIT NONE

INTEGER   :: ND
INTEGER*8 :: NINJ
DOUBLE PRECISION, DIMENSION(1:ND)   :: up,gFp
DOUBLE PRECISION, DIMENSION(1:NINJ) :: up0,Rp
DOUBLE PRECISION                    :: gFp1,gFp2,u1,u2
INTEGER, DIMENSION(1)               :: minu

! Loop indexes
INTEGER :: ic,iu

!***********************************************************************

Rp=0.0

CALL RANDOM_NUMBER(Rp)
Rp=Rp*0.999

DO ic=1,NINJ  ! PPC

minu=minloc(abs(gFp-Rp(ic)))
iu=minu(1)

IF (iu.EQ.ND) THEN
iu=iu-1
END IF

gFp1=gFp(iu)
gFp2=gFp(iu+1)

u1=up(iu)
u2=up(iu+1)

up0(ic)=(Rp(ic)*(u2-u1)-u2*gFp1+u1*gFp2)/(gFp2-gFp1)

ENDDO

END SUBROUTINE GEN_UP

!***********************************************************************
! Subroutine GEN_PS
! This subroutine generates a random distribution of particle momenta ps.
!
! INPUT: 
! - ps0: Particle momentum ps=uperp/sqrt(1+up*up)
! - ps: reference momentum array
! - up0: Particle momentum parallel to the drift velocity
! - up: reference momentum array
! - gFs: pre-calculated integrated distribution function
! - ND: number of elements in the array ug and gFg
!
! OUTPUT: distribution of ps0
!***********************************************************************

SUBROUTINE GEN_PS(ps0,ps,up0,up,gFs,ND,NINJ)

IMPLICIT NONE

INTEGER   :: ND
INTEGER*8 :: NINJ
DOUBLE PRECISION, DIMENSION(1:ND)      :: ps,up
DOUBLE PRECISION, DIMENSION(1:ND,1:ND) :: gFs
DOUBLE PRECISION, DIMENSION(1:NINJ)    :: ps0,up0,Rp
DOUBLE PRECISION                       :: F11,F12,F21,F22,fp,fq
INTEGER, DIMENSION(1)                  :: minu,minf

! Loop indexes
INTEGER :: ic,i1,i2

!***********************************************************************

Rp=0.0

CALL RANDOM_NUMBER(Rp)
Rp=Rp*0.999

DO ic=1,NINJ  ! PPC

minu=minloc(abs(up-up0(ic)))
i1=minu(1)

minf=minloc(abs(gFs(:,i1)-Rp(ic)))
i2=minf(1)

IF (i1.EQ.ND) THEN
i1=i1-1
END IF

IF (i2.EQ.ND) THEN
i2=i2-1
END IF

F11=gFs(i2,i1)
F12=gFs(i2,i1+1)
F21=gFs(i2+1,i1)
F22=gFs(i2+1,i1+1)

fq=(up0(ic)-up(i1))/(up(i1+1)-up(i1))
fp=(Rp(ic)-(1.0-fq)*F11-fq*F12)/((1.0-fq)*(F21-F11)+fq*(F22-F12))

ps0(ic)=fp*(ps(i2+1)-ps(i2))+ps(i2)

ENDDO

END SUBROUTINE GEN_PS

!***********************************************************************
! Subroutine INJ_DRIFT_MAXWELLIAN
! This subroutine injects a wind of particles along field lines where the
! NOZZLE field boundary condition is specified.
!
! INPUT: 
! - speed: Field-parallel bulk drift velocity of the wind in units of c 
! - theta: Temperature of plasma, in units of mc^2/k, in wind rest frame.
!     theta can be equal to 0, but numbers bigger than 0 and smaller than
!     something like 0.01 may cause strange behavior
! - n0: Target number density of the wind
!     The actual density will be close to, but possibly not exactly equal to,
!     this value
!     Specifying n0 <= 0 means no injection
! - [x|y|z]minp: Lower spatial boundary for each domain
! - up,gFp,ps,gFs,ND: Particle distribution generator quantities
!
! OUTPUT: Particle distribution function (weights set, too)
! - NINJ: The number of injected particles
! - pcl_inj: Injected distribution function of the particles
!   pcl_inj is never deallocated, but will be resized as necessary to
!   contain all of the injected particles.
!***********************************************************************
SUBROUTINE INJ_DRIFT_MAXWELLIAN(pcl_inj,NINJ,n0,speed,theta,xminp,yminp,zminp,up,gFp,ps,gFs,ND)

IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN) :: speed,theta,n0
DOUBLE PRECISION, INTENT(IN) :: xminp,yminp,zminp

INTEGER*8, INTENT(OUT) :: NINJ

DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: pcl_inj

! Pre-calculated CDF arrays
DOUBLE PRECISION, DIMENSION(:),   INTENT(IN) :: up,gFp,ps
DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: gFs
! Number of elements in the pre-calculated cumulative distribution functions
INTEGER,                          INTENT(IN) :: ND

!******************************
! Variables local to subroutine
!******************************

! Unit vectors defined by local field direction
DOUBLE PRECISION :: bhatx,bhaty,bhatz,bmag
DOUBLE PRECISION :: e1hatx,e1haty,e1hatz
DOUBLE PRECISION :: e2hatx,e2haty,e2hatz

DOUBLE PRECISION, DIMENSION(1:NCXP,1:NCYP)  :: ninj_per_cell_dbl
DOUBLE PRECISION, DIMENSION(1:NCXP,1:NCYP)  :: rand_draws, fractional_ptcls
INTEGER,          DIMENSION(1:NCXP,1:NCYP)  :: ninj_per_cell

DOUBLE PRECISION :: shutoff_factor, targ_dens

DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: up0,ps0,uperp0,phi0

DOUBLE PRECISION :: x0c,y0c,z0c

! Loop indices
INTEGER   :: ic,ix,iy,iz
INTEGER*8 :: ip

!******************************
! Begin code
!******************************

CALL INIT_RANDOM_SEED()

iz = 1

NINJ = 0
IF (ALLOCATED(pcl_inj)) THEN
    pcl_inj = 0.0
END IF

! Don't inject new particles if n0 <= 0d0
IF (n0 > 0d0) THEN

!*****************************************************
! Figure out how many particles to inject in each cell
!*****************************************************
DO iy=1,NCYP
DO ix=1,NCXP

    CALL INITIAL_BCOMPONENTS(xyeep(ix),yyeep(iy),zyeep(iz),bhatx,bhaty,bhatz,bmag)
    ninj_per_cell_dbl(ix,iy) = over_inject*PPC*(speed*c*dt/dz)*ABS(bhatz)

ENDDO
ENDDO

! The number of injected macro-particles per cell calculated above is not
! necessarily an integer. Here, we convert it into an integer by flipping
! a coin to see if the non-integer part of the number also becomes injected
! as a macroparticle. In this way, the number of particles injected per
! timestep in cell ir will be equal, on average, to ninj_per_cell_dbl(ir).
fractional_ptcls=ninj_per_cell_dbl-FLOOR(ninj_per_cell_dbl)
CALL RANDOM_NUMBER(rand_draws)

! The 0.0*x+1.0 is a hack to get an array of ones of the same dimensions as x.
ninj_per_cell=FLOOR(ninj_per_cell_dbl)+&
    0.5*(1.0+SIGN(0.0*fractional_ptcls+1.0,fractional_ptcls-rand_draws))

! Total number of injected particles
NINJ=SUM(ninj_per_cell)

!***********************************************************************
! Populate pcl_inj
!***********************************************************************
IF (ALLOCATED(pcl_inj)) THEN
    IF (SIZE(pcl_inj,DIM=2) < NINJ) THEN
        DEALLOCATE(pcl_inj)
        ALLOCATE(pcl_inj(1:7,1:NINJ))
    ENDIF
ELSE
    ALLOCATE(pcl_inj(1:7,1:NINJ))
ENDIF

pcl_inj=0.0

!***********************************************************************
! Sample particle velocities
!***********************************************************************
ALLOCATE(up0(NINJ))
ALLOCATE(ps0(NINJ))
ALLOCATE(uperp0(NINJ))
ALLOCATE(phi0(NINJ))

IF (theta.GT.0) THEN

    up0=0.0

    CALL GEN_UP(up0,up,gFp,ND,NINJ)

    ! Definition of the initial perpendicular momentum
    ps0=0.0

    CALL GEN_PS(ps0,ps,up0,up,gFs,ND,NINJ)

    uperp0=ps0*sqrt(1.0+up0*up0)

    ! Definition of the initial perpendicular momentum azimuthal angle phi0 as a
    ! uniform random value
    phi0=0.0

    CALL RANDOM_NUMBER(phi0)
    phi0=phi0*2.0*pi

ELSE
    ! Do theta.EQ.0 separately
    ! This avoids bogus from the routines used to sample the drifting particle
    ! distribution at low temperatures
    up0=speed
    uperp0=0.0
    phi0=0.0
    CALL RANDOM_NUMBER(phi0)
    phi0=phi0*2.0*pi

ENDIF

ip = 1
DO iy=1,NCYP
DO ix=1,NCXP
DO ic=1,ninj_per_cell(ix,iy)

!***********************************************************************
! Sample particle positions
!***********************************************************************
! Definition of the positions as uniform random values
x0c = 0.0
CALL RANDOM_NUMBER(x0c)
x0c = xgp(ix) + x0c*dx

y0c = 0.0
CALL RANDOM_NUMBER(y0c)
y0c = ygp(iy) + y0c*dy

z0c = 0.0
CALL RANDOM_NUMBER(z0c)
z0c = zgp(iz) + z0c*dz

pcl_inj(1,ip)=x0c
pcl_inj(2,ip)=y0c
pcl_inj(3,ip)=z0c

!***********************************************************************
! Rotate particle velocities along local magnetic field
!***********************************************************************
IF (ISNAN(uperp0(ip))) THEN
    PRINT *, "Found uperp0 NaN in injector"
ENDIF
IF (ISNAN(up0(ip))) THEN
    PRINT *, "Found up0 NaN in injector"
ENDIF

! Retrieve a unit vector, bhat = (bhatx,bhaty,bhatz), along the local magnetic field
CALL INITIAL_BCOMPONENTS(x0c,y0c,z0c,bhatx,bhaty,bhatz,bmag)
! Construct an orthonormal basis using bhat. Calculate two unit vectors e1 and
! e2 such that e1 \times e2 = bhat. Then interpret phi0(ip) as the angle of the
! perpendicular velocity in the e1-e2 plane from the e1 axis.
CALL ORTHONORMAL_TRIAD(bhatx,bhaty,bhatz,e1hatx,e1haty,e1hatz,e2hatx,e2haty,e2hatz)

pcl_inj(4,ip) = uperp0(ip)*( COS(phi0(ip))*e1hatx + SIN(phi0(ip))*e2hatx ) + up0(ip)*bhatx
pcl_inj(5,ip) = uperp0(ip)*( COS(phi0(ip))*e1haty + SIN(phi0(ip))*e2haty ) + up0(ip)*bhaty
pcl_inj(6,ip) = uperp0(ip)*( COS(phi0(ip))*e1hatz + SIN(phi0(ip))*e2hatz ) + up0(ip)*bhatz

! Ensure that particles are injected into the simulation
IF (pcl_inj(6,ip)<0d0) THEN
    pcl_inj(4,ip) = -pcl_inj(4,ip)
    pcl_inj(5,ip) = -pcl_inj(5,ip)
    pcl_inj(6,ip) = -pcl_inj(6,ip)
ENDIF

shutoff_factor = ROTATION_PROFILE(x0c,y0c)
! Assign weights so that, downstream of the injection region where the average
! number of macro-particles per cell is (over_inject*PPC), the equivalent physical
! number density is (shutoff_factor*rate*)n0.
pcl_inj(7,ip) = shutoff_factor*rate*n0*dx*dy*dz/(over_inject*PPC)

ip = ip + 1

ENDDO
ENDDO
ENDDO

DEALLOCATE(up0)
DEALLOCATE(ps0)
DEALLOCATE(uperp0)
DEALLOCATE(phi0)

END IF  ! n0 > 0d0

END SUBROUTINE INJ_DRIFT_MAXWELLIAN

!***********************************************************************
! Subroutine ORTHONORMAL_TRIAD
! Given one (unit) vector, construct two transverse unit vectors to complete
! an orthonormal triad.
!
! INPUT:
! - nx, ny, nz: Components of a vector defining one axis of the triad.
!     Requirement: nx^2 + ny^2 + nz^2 > 0.0,
!     But can have nx^2 + ny^2 + nz^2 /= 1.0
!
! OUTPUT:
! - b1x, b1y, b1z, b2x, b2y, b2z: Components of the other two transverse
!     unit vectors in the triad. Properties:
!     norm(b1) == norm(b2) == norm(n / norm(n)) == 1
!     cross(b1, b2) == n / norm(n)
!***********************************************************************
SUBROUTINE ORTHONORMAL_TRIAD(nx, ny, nz, b1x, b1y, b1z, b2x, b2y, b2z)

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN)  :: nx, ny, nz
  DOUBLE PRECISION, INTENT(OUT) :: b1x, b1y, b1z, b2x, b2y, b2z

  DOUBLE PRECISION :: singularity_cutoff = -0.999999
  DOUBLE PRECISION :: denom_fac, nhatx, nhaty, nhatz, nmag, nmaginv

  ! Begin code
  nmag = SQRT(nx*nx + ny*ny + nz*nz)
  IF (nmag == 0.0) THEN
     PRINT *, "ORTHONORMAL_TRIAD received zero-magnitude vector."
     PRINT *, "You're in trouble..."
  ENDIF
  nmaginv = 1.0 / nmag
  nhatx = nx * nmaginv
  nhaty = ny * nmaginv
  nhatz = nz * nmaginv
  IF (nhatz < singularity_cutoff) THEN
     b1x = 0.0
     b1y = -1.0
     b1z = 0.0
     b2x = -1.0
     b2y = 0.0
     b2z = 0.0
  ELSE
     denom_fac = 1.0 / (1.0 + nhatz)
     b1x = 1.0 - nhatx * nhatx * denom_fac
     b1y = -nhatx * nhaty * denom_fac
     b1z = -nhatx
     b2x = -nhatx * nhaty * denom_fac
     b2y = 1.0 - nhaty * nhaty * denom_fac
     b2z = -nhaty
  ENDIF
  
END SUBROUTINE ORTHONORMAL_TRIAD

!***********************************************************************
! init_random_seed() subroutine enables to avoid the repeating series of
! number given by random_number.
! Reference: http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
!***********************************************************************

SUBROUTINE INIT_RANDOM_SEED()

INTEGER :: i, n, clock
INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
IF (RANDOMIZE.EQV..TRUE.) THEN
CALL RANDOM_SEED(size = n)
ALLOCATE(seed(n))
        
CALL SYSTEM_CLOCK(COUNT=clock)
         
seed = clock + 37 * (/ (i - 1, i = 1, n) /)
CALL RANDOM_SEED(PUT = seed)
          
DEALLOCATE(seed)
ENDIF

END SUBROUTINE INIT_RANDOM_SEED

!***********************************************************************

END MODULE MOD_INITIAL
