subroutine KFeSe2vL()
  use workplace
  implicit none

  real*8 Q
  integer iz

  DataReady=.false.
  QuickSearch=.false.
  SkipInterOrbitPair=.false.
  AppendBreakPoint=.false.
  Appendkmesh=.true.; kmeshfile='KFS2vLkmesh.dat'
  useprojectors=.true.;  Appendprojectors=.false.; projectorfile='KFS1LPrj.dat'
  useprojectors=useprojectors.and.(.not.QuickSearch)

    norb=10; natom=2    
	model%norb=norb; model%natom=natom; model%sample='KFS2vL'
    print*,'The sample chosen is ',model%sample 

	allocate(model%ra(3,natom)); model%ra(:,1)=0; model%ra(:,2)=(/0,0,1/)*0.5
	model%a=(/1,0,0/); model%b=(/0,1,0/); model%c=(/0,0,1/)
	model%periodiclayers=.true.

    allocate(model%orbits(norb)); model%orbits(1:5)=(/'A1g','E1u','E2u','B1g','B2g'/)
	model%orbits(6:10)=model%orbits(1:5)	
    
	usegroup=.true.;  group='D4h';  ng=8
    model%group=group; model%ng=ng

    model%U=4/2; model%JH=0.4/2; model%Uab=2/2 !model%U-2*model%JH; 
	model%mu=9.27

	model%zerotz=.false.
    allocate(model%t(5,5,-2:2,-2:2,-1:1))
    open(10,file='KFeSe2vL.model')
    read(10,*)model%t; close(10)
    do iz=-1,1; call basisrotate(model%t(:,:,:,:,iz)); end do
	Q=0; call checksymmetry(Q);	!call FermiSurf(norb,Q,model)

    nkx=40; nky=40      
    nLform=6; Lcontact=2
    allocate(Lformtable(nLform)); Lformtable=(/'A1g','A1g','B1g','A1g','B2g','A1g'/)
    allocate(Lformbasis(3,nLform)); Lformbasis(:,1)=0
	Lformbasis(:,2)=(/1,0,0/); Lformbasis(:,3)=(/1,0,0/)
	Lformbasis(:,4)=(/1,1,0/); Lformbasis(:,5)=(/1,1,0/)
	Lformbasis(:,6)=(/0.,0.,0.5/)

	call setLformfunc(ng)
	call setOformfunc()

    SkipInterOrbitPair=.false.
	guess='A1g';  call setMolecules()
	
    nw=400; wmax=100
    wmin=1.e-4; wir=1.e-3
    diverge=40

  return
end subroutine KFeSe2vL

subroutine KFeSe2vLhk(norb,kv,hk,model,mf,slvbsn)
  use standard_derived_types
  implicit none
  integer norb
  real*8 kv(2)
  type (modelblock), target :: model
  type (mfblock) :: mf
	logical slvbsn
  complex*16, pointer, dimension (:,:) :: t
  complex*16 hk(norb,norb),one,phase,error

  integer ix,iy,iz,iorb
  real*8 bond(2),pi

  pi=asin(1.d0)*2; one=cmplx(0,1)

  hk=0; do iorb=1,10; hk(iorb,iorb)=-model%mu; end do
  do iz=-1,1; do iy=-2,2; do ix=-2,2; bond=(/ix,iy/)
     t=>model%t(:,:,ix,iy,iz)
	 phase=exp(-one*pi*sum(kv*bond))
     if(iz==0)then
	   hk(1:5,1:5)=hk(1:5,1:5)+t*phase
	   hk(6:10,6:10)=hk(6:10,6:10)+t*phase
	 else 
	   hk(1:5,6:10)=hk(1:5,6:10)+t*phase
	   hk(6:10,1:5)=hk(6:10,1:5)+t*phase
	 end if
  end do; end do; end do

  !do ix=1,10; do iy=ix+1,10
  !   error=hk(ix,iy)-conjg(hk(iy,ix))
  !	 if(abs(error)>1.e-5)then
  !	   print*,ix,iy
  !	   print*,hk(ix,iy)
  !	   print*,conjg(hk(iy,ix))	   
  !	   pause 'hermitian error'
  !	 end if
  !end do; end do

  return
end subroutine KFeSe2vLhk
