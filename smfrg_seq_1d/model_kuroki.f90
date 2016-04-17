subroutine Kuroki()
  use workplace
  implicit none

  real*8 Q,cutoff,factor
  integer nband,np
  logical dorpa,appendrpax0

  dorpa=.false.
  appendrpax0=.false.

  if(.true.)then
    appendkmesh=.true.
    appendfs=.true.
    appendX0=.true.
  else
    appendkmesh=.false.
    appendfs=.false.
    appendx0=.false.
  end if

  useX0=.true.
  BCSflow=.false.
  QuickSearch=.false.
  outputXq=.false.
  DataReady=.false.
  AppendBreakPoint=.false.
  SkipInterOrbitPair=.false.
  kmeshfile='KurokiKmesh.dat'
  useprojectors=.false.; Appendprojectors=.true.; projectorfile='KurokiPrj.dat'
  useprojectors=useprojectors.and.(.not.QuickSearch)   !projectors are not used in the QuickSearch mode.

    norb=5; natom=1; nambu=10; model%n4logqmesh=1
    model%norb=norb; model%natom=natom; model%sample='Kuroki'
    print*,'The sample chosen is ',model%sample 
	allocate(model%ra(3,natom)); model%ra=0   
	allocate(model%orbits(norb))
  model%orbits=(/'A1g','E2g','E1g','B1g','B2g'/)     
	model%periodiclayers=.false.
  
  model%sample='Kuroki'
  open(10,file='frg.input')
  read(10,*)model%U,model%JH,model%mu,np
  close(10)
	!model%U=2.0
	!model%JH=0.4
	if(dorpa)then
    print*,'U='
	  read(*,*)model%U
	  print*,'JH='
	  read(*,*)model%JH
    end if
	model%Uab=model%U-2*model%JH  
  
  
  
  allocate(model%t(5,5,-4:4,-4:4,1))     
  open(10,file='Kuroki5d.model')
  read(10,'(2e16.8)')model%t(:,:,:,:,1); close(10)

  nkx=100
  nky=100

  usegroup=.true.  
  group='D4h'; ng=8
	model%group=group; model%ng=ng
	Q=0; call checksymmetry(Q)
	call FermiSurf(norb,Q,model)
  call BandStructure(norb,Q,model)

	nband=3; model%nband=nband
	allocate(model%npatch(nband),model%nativeband(nband),model%pocket(nband),model%mushift(nband))
	model%nativeband=(/2,3,4/)
  if(np==4)then
    model%pocket=(/'G_','G_','XY'/); model%npatch=(/32,32,64/)
  end if
  if(np==5)then
    model%pocket=(/'G_','GM','XY'/); model%npatch=(/32,64,64/)
  end if
  model%mushift=0.
	call meshfermisurface(model)

  
	nLform=9; Lcontact=2
  allocate(Lformtable(nLform)); Lformtable='A1g' !(/'A1g','A1g','B1g','A1g','B2g'/)
  allocate(Lformbasis(3,nLform)); 
  Lformbasis(:,1)=0
	Lformbasis(:,2)=(/1,0,0/)
  Lformbasis(:,3)=(/-1,0,0/)
  Lformbasis(:,4)=(/0,1,0/)
  Lformbasis(:,5)=(/0,-1,0/)
  Lformbasis(:,6)=(/1,1,0/)
  Lformbasis(:,7)=(/-1,-1,0/)
  Lformbasis(:,8)=(/1,-1,0/)
  Lformbasis(:,9)=(/-1,1,0/)

	call setLformfunc(1)
	call setOformfunc()
	guess='***';  call setMolecules()

    nw=400; wmax=100
    wmin=1.e-4; wir=1.e-5
    diverge=40

    if(.not.dorpa)return

    call writeinputs()
    call setupworkplace()
    cutoff=1.e-3; call rpa(cutoff,appendrpax0,'kurokiX0.dat'); stop

  return
end subroutine Kuroki

subroutine Kuroki5dhk(norb,kv,Q,hk,model,mf,slvbsn)
  use standard_derived_types
  implicit none
  integer norb
  real*8 kv(2),Q
  complex*16 hk(norb,norb)
  type (modelblock), target :: model
  type (mfblock) :: mf
	logical slvbsn

  complex*16, pointer, dimension (:,:) :: t
  integer ix,iy,iorb,bondmax
  complex*16 one
  real*8 pi,vk(2),error

  one=cmplx(0,1); pi=asin(1.d0)*2
  vk=kv*pi

  hk=0
  do iorb=1,5; hk(iorb,iorb)=-model%mu; end do

  bondmax=2
  do iy=-bondmax,bondmax; do ix=-bondmax,bondmax
     t=>model%t(:,:,ix,iy,1)
	 hk=hk+t*exp( one*( ix*vk(1)+iy*vk(2) ) )
  end do; end do
 
  return
  error=sum( abs( aimag(hk) ) )
  if(error>1.e-4)pause 'imag(hk)/=0?'

  return
end subroutine Kuroki5dhk
