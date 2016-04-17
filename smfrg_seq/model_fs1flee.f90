subroutine FS1FLee()
  use workplace
  implicit none
  real*8 Q,tt
  integer ix,iy,a,b,nkf,ik,nband

  usex0=.true.
  appendx0=.true.
  appendfs=.true.
  appendkmesh=.true.

  DataReady=.false.
  QuickSearch=.false.
  SkipInterOrbitPair=.false.
  AppendBreakPoint=.false.
  kmeshfile='FS1FLeekmesh.dat'
  useprojectors=.false.;  Appendprojectors=.false.; projectorfile='FS1FLeePrj.dat'
  useprojectors=useprojectors.and.(.not.QuickSearch)

  norb=5; natom=1
	model%norb=norb; model%natom=natom; model%sample='FS1FLee'
  print*,'The sample chosen is ',model%sample

	allocate(model%ra(3,natom)); model%ra=0   
	allocate(model%orbits(norb))
  model%orbits=(/'A1g','E2g','E1g','B1g','B2g'/)     

	open(10,file='frg.input')
  read(10,*)model%U,model%JH,model%mu,model%Veph
	model%Uab=model%U-2*model%JH
	
	model%ephcoupled=.true.
	!model%Veph=1.0
	if(model%ephcoupled.and.model%Veph<1.e-6)model%ephcoupled=.false.

  allocate(model%t(5,5,-2:2,-2:2,1))
  model%t=0
  open(10,file='FS1FLee_new.model')
    read(10,*)model%t
	close(10)
	call basisrotate(model%t(:,:,:,:,1))

	model%nband=1 
  allocate(model%nativeband(1),model%npatch(1),model%pocket(1),model%mushift(1))
	model%nativeband=4; model%npatch=120; model%pocket='XY'; model%mushift=0
  call meshfermisurface(model)

  usegroup=.true.  
  group='D4h'; ng=8
  model%group=group; model%ng=ng
  Q=0; call checksymmetry(Q)
	call FermiSurf(norb,Q,model); 
  call bandstructure(norb,Q,model)


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
    wmin=1.e-4; wir=1.e-5
    diverge=100

  return
end subroutine FS1FLee

subroutine FS1FLeehk(norb,kv,Q,hk,model,mf,slvbsn)
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
  real*8 pi,vk(2),error,delta

  one=cmplx(0,1); pi=asin(1.d0)*2
  vk=kv*pi

  hk=0; delta=0
  do iorb=1,5; hk(iorb,iorb)=-model%mu; end do

  bondmax=2
  do iy=-bondmax,bondmax; do ix=-bondmax,bondmax
     t=>model%t(:,:,ix,iy,1)
	 hk=hk+t*exp( one*( ix*vk(1)+iy*vk(2) ) )
  end do; end do
  
  hk(2,2)=hk(2,2)+delta; hk(3,3)=hk(3,3)+delta
  return
  error=sum( abs( aimag(hk) ) )
  if(error>1.e-4)pause 'imag(hk)/=0?'

  return
end subroutine 

 
