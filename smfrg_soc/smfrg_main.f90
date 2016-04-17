module standard_derived_types

  type generation
    integer np                                    !number of k-vectors at a particular generation
    real*8 w                                      !energy scale at a particular generation
    real*8 dp                                     !half of the linear size of the box containing a k-vector
    integer, allocatable, dimension (:) :: info   !info=>starting index at the next generation /0 if a k-point has/hasn't off-generations
    real*8, allocatable, dimension (:,:) :: p     !the collection of k-vectors at this generation
  end type generation

  type meshconfig 
    integer ng                                         !number of mesh generations
    integer npoints                                    !number of end generation points
    type (generation), allocatable, dimension (:) :: g !generation collection in a log mesh
  end type meshconfig

  type Lformconfig
    character*3 form
    integer nr,ndim
	  real*8, allocatable, dimension (:,:) :: r
	  real*8, allocatable, dimension (:) :: f
  end type Lformconfig

  type Mformconfig
    integer atom,nbatom
    integer spin,nbspin
    integer orb,nborb
    integer Lform
    character*3 form
    character*2 Sform
    type (Lformconfig) :: Lformfunc
  end type Mformconfig

  type modelblock
    integer norb,natom                                !number of (total) orbits, atoms
	  real*8, allocatable, dimension (:,:) :: ra        !atom positions within a unit cell
  	real*8, dimension (3) :: a,b,c                    !primitive translation vectors
  	real*8, dimension (2) :: ka,kb                    !G vectors in BZ
    character*3 group                                 !point group of the model
    integer ng                                        !number of group elements 
	  real*8 U,Vnn,VL2L,Vbonding,Vantibonding,Jnn,mu,la,t1 !parameters of the system
	  real*8, allocatable, dimension (:) :: detune      !energy detune for unequivalent atoms
    character*10 sample                               !sample = grphn, grpht, trngl, kagme, dice
	  logical propergauge
  end type modelblock

  type edgemodelblock
    integer natom                                !number of atoms
	real*8, allocatable, dimension (:,:) :: ra   !atom positions within a unit cell
	real*8, dimension (3) :: a,b,c               !primitive translation vectors
	real*8, dimension (2) :: ka,kb               !G vectors in BZ

	character*10 sample
	logical periodicslab              !periodic slab or open-boundary slab
    integer nb                        !direction of the boundary normal
    integer Nq                        !number of conserved transverse momenta 
    integer Lb                        !number of sites along an open-boundary slice
    integer norb                      !number of effective orbits in the model (including spin)
    integer nambu                     !dimension of nambu space (nambu=norb for normal state)
    integer ndim                      !hilbert space dimension of the open-boundary slice with a conserved transverse momentum
  end type edgemodelblock

  type outputblock
	real*8 Vmax(3),q(2,2)
	real*8 dope
  real*8 fm
  integer sc
  end type outputblock


end module standard_derived_types


module workplace
  use standard_derived_types

  integer norb,natom

  integer nLform
  character*3, allocatable, dimension (:) :: Lformtable
  real*8, allocatable, dimension (:,:) :: LformBasis

  integer ndim
  type (Mformconfig), allocatable, target, dimension (:) :: Mformfunc

  type (modelblock), target :: model
  type (edgemodelblock), target :: edgemodel
  type (outputblock) :: output

  integer nq
  real*8, allocatable, dimension (:,:) :: qv 
  type (meshconfig) :: qmesh

  integer nk
  real*8, allocatable, dimension (:,:) :: kv
  type (meshconfig) :: kmesh
  logical appendkmesh
  character*20 kmeshfile

  integer nkx,nky
  real*8, allocatable, dimension (:,:,:) :: ek

  integer nR,Lcontact
  integer, allocatable, dimension (:,:) :: R

  character*3 group
  integer ng
  logical usegroup
  integer, allocatable, dimension (:,:) :: indexgroup
  integer, allocatable, dimension (:,:) :: indexqg
  complex*16, allocatable, dimension (:,:) :: ppXg,phXg
 
  real*8 pi,twopi
  complex*16 one

  integer nw
  real*8 wmax,wmin,wir,diverge   !energy scales and divergence criterion

  complex*16, allocatable, dimension (:,:,:) :: P,C,D

  logical DataReady,AppendBreakPoint,QuickSearch,bcsflow

  logical useprojectors,AppendProjectors
  character*20 Projectorfile
  integer nP2C,nP2D,nC2P,nC2D,nD2C,nD2P
  complex*16, allocatable, dimension (:) :: P2C,P2D,C2P,C2D,D2C,D2P
  integer, allocatable, dimension (:,:) :: infoP2C,infoP2D,infoC2P,infoC2D,infoD2C,infoD2P

end module workplace


!=======================================
program main
use workplace

call definejob()
call setupworkplace()
call runfrg()

stop
end
