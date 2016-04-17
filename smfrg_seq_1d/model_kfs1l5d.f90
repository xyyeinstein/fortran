subroutine KFeSe1L5d()
  use workplace
  implicit none
  real*8 Q,tt
  integer ix,iy,a,b,error

  DataReady=.true.
  QuickSearch=.false.
  SkipInterOrbitPair=.false.
  AppendBreakPoint=.false.
  Appendkmesh=.true.; kmeshfile='KFS1Lkmesh.dat'
  useprojectors=.true.;  Appendprojectors=.true.; projectorfile='KFS1LPrj.dat'
  useprojectors=useprojectors.and.(.not.QuickSearch)

    norb=5; natom=1
	model%norb=norb; model%natom=natom; model%sample='KFS1L'
    print*,'The sample chosen is ',model%sample 
    
	allocate(model%ra(3,natom)); model%ra=0   
	allocate(model%orbits(norb))
    model%orbits=(/'A1g','E1u','E2u','B1g','B2g'/)
	model%a=(/1,0,0/); model%b=(/0,1,0/); model%c=(/0,0,1/)

	model%U=2; model%JH=0.245; model%Uab=model%U-2*model%JH; model%mu=9.22-0.25
  
    allocate(model%t(5,5,-2:2,-2:2,1)); model%t=0
    open(10,file='KFeSe5d0tz.model')
    do while ( .true. )
       read(10,*,iostat=error)a,b,ix,iy,tt
       model%t(a,b,-ix,-iy,1)=tt
       if(error/=0)exit
    end do; close(10)
    call basisrotate(model%t(:,:,:,:,1))

    usegroup=.true.  
    group='C4v'; ng=8
	model%group=group; model%ng=ng

	Q=0; call checksymmetry(Q);	call FermiSurf(norb,Q,model); call bandstructure(norb,Q,model)  

    nkx=40; nky=40  
	nLform=5; Lcontact=2
    allocate(Lformtable(nLform)); Lformtable=(/'A1g','A1g','B1g','A1g','B2g'/)
    allocate(Lformbasis(3,nLform)); Lformbasis(:,1)=0
	Lformbasis(:,2)=(/1,0,0/); Lformbasis(:,3)=(/1,0,0/)
	Lformbasis(:,4)=(/1,1,0/); Lformbasis(:,5)=(/1,1,0/)
	call setLformfunc(ng)

	call setOformfunc()

	guess='A1g';  call setMolecules()

    nw=400; wmax=100
    wmin=1.e-4; wir=1.e-3
    diverge=100

  return
end subroutine KFeSe1L5d

subroutine KFeSe1L5dhk(norb,kv,Q,hk,model,mf,slvbsn)
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
end subroutine KFeSe1L5dhk
