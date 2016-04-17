subroutine KFeSe1Lav()
  use workplace
  implicit none

  real*8 Q,factor
  integer iz

  DataReady=.false.
  useX0=.true.
  appendX0=.true.
  appendkmesh=.true.
  appendfs=.true.
  QuickSearch=.false.
  SkipInterOrbitPair=.false.
  AppendBreakPoint=.false.
  kmeshfile='KFSLavkmesh.dat'
  useprojectors=.false.;  Appendprojectors=.false.; projectorfile='KFS1LavPrj.dat'
  useprojectors=useprojectors.and.(.not.QuickSearch)

    norb=5; natom=1    
	model%norb=norb; model%natom=natom; model%sample='KFS1Lav'
    print*,'The sample chosen is ',model%sample 

	allocate(model%ra(3,natom)); model%ra(:,1)=0
	model%a=(/1,0,0/); model%b=(/0,1,0/); model%c=(/0,0,1/)

    allocate(model%orbits(norb)); model%orbits(1:5)=(/'A1g','E1u','E2u','B1g','B2g'/)
    
	usegroup=.true.;  group='D4h';  ng=8
    model%group=group; model%ng=ng

    model%U=2; model%JH=0.245; model%Uab=model%U-2*model%JH; model%mu=9.223

	model%ephcoupled=.true.
	model%Veph=2
	if(model%ephcoupled.and.model%Veph<1.e-6)pause 'Error: Veph=0 while ephcoupled @ KFeSe1Lav'


    allocate(model%t(5,5,-2:2,-2:2,1))
    open(10,file='KFeSe1Lav.model'); read(10,*)model%t; close(10)

	
    factor=1.0
	model%t(2:3,2:3,1,0,1)=model%t(2:3,2:3,1,0,1)*factor
	model%t(2:3,2:3,-1,0,1)=model%t(2:3,2:3,-1,0,1)*factor
	model%t(2:3,2:3,0,1,1)=model%t(2:3,2:3,0,1,1)*factor
	model%t(2:3,2:3,0,-1,1)=model%t(2:3,2:3,0,-1,1)*factor

    factor=0.8
	model%t(5,5,1,0,1)=model%t(5,5,1,0,1)*factor
	model%t(5,5,-1,0,1)=model%t(5,5,-1,0,1)*factor
	model%t(5,5,0,1,1)=model%t(5,5,0,1,1)*factor
	model%t(5,5,0,-1,1)=model%t(5,5,0,-1,1)*factor
    
	factor=1.
	model%t(2:3,2:3,1,1,1)=model%t(2:3,2:3,1,1,1)*factor
	model%t(2:3,2:3,-1,1,1)=model%t(2:3,2:3,-1,1,1)*factor
	model%t(2:3,2:3,1,-1,1)=model%t(2:3,2:3,1,-1,1)*factor
	model%t(2:3,2:3,-1,-1,1)=model%t(2:3,2:3,-1,-1,1)*factor

	Q=0
	call checksymmetry(Q)
	call bandstructure(norb,Q,model)
    call FermiSurf(norb,Q,model)
	
    nkx=40; nky=40      
    nLform=5; Lcontact=2
    allocate(Lformtable(nLform)); Lformtable=(/'A1g','A1g','B1g','A1g','B2g'/)
    allocate(Lformbasis(3,nLform)); Lformbasis(:,1)=0
	Lformbasis(:,2)=(/1,0,0/); Lformbasis(:,3)=(/1,0,0/)
	Lformbasis(:,4)=(/1,1,0/); Lformbasis(:,5)=(/1,1,0/)

	call setLformfunc(ng)
	call setOformfunc()

    SkipInterOrbitPair=.false.
	guess='***';  call setMolecules()
	
    nw=400; wmax=100
    wmin=1.e-4; wir=1.e-3
    diverge=40

  return
end subroutine KFeSe1Lav

subroutine KFeSe1Lavhk(norb,kv,Q,hk,model,mf,slvbsn)
  use standard_derived_types
  implicit none
  integer norb
  real*8 kv(2),Q
  type (modelblock), target :: model
  type (mfblock) :: mf
	logical slvbsn
  complex*16, pointer, dimension (:,:) :: t
  complex*16 hk(norb,norb),one,phase,error

  integer ix,iy,iz,iorb
  real*8 bond(2),pi

  pi=asin(1.d0)*2; one=cmplx(0,1)

  hk=0; do iorb=1,5; hk(iorb,iorb)=-model%mu; end do
  do iy=-2,2; do ix=-2,2; bond=(/ix,iy/)
     t=>model%t(:,:,ix,iy,1)
	 phase=exp(-one*pi*sum(kv*bond))
	 hk(1:5,1:5)=hk(1:5,1:5)+t*phase
  end do; end do

  return
end subroutine KFeSe1Lavhk
