subroutine WangFa5d()
  use workplace
  implicit none

  real*8 Q

  DataReady=.false.
  AppendBreakPoint=.false.
  Appendkmesh=.true.; kmeshfile='WangFakmesh.dat'

    norb=5; natom=1
	model%norb=norb; model%natom=natom; model%sample='WangFa'
    print*,'The sample chosen is ',model%sample

	allocate(model%ra(3,natom)); model%ra=0   
	allocate(model%orbits(norb))
    model%orbits=(/'A1g','E2g','E1g','B1g','B2g'/)     
	model%periodiclayers=.false.

	model%sample='WangFa'; model%U=4/2; model%JH=0.54/2; model%Uab=model%U-2*model%JH;  !B1g superconducting at 0.017 eV    
    model%mu=7.5

    allocate(model%t(5,5,-2:2,-2:2,1))
    open(10,file='WangFa5d.dat')
    read(10,*)model%t(:,:,:,:,1); close(10)
    call basisrotate(model%t(:,:,:,:,1))

    usegroup=.true.  
    group='D4h'; ng=8
    model%group=group; model%ng=ng
    Q=0; call checksymmetry(Q);	call FermiSurf(norb,Q,model); !pause 'continue?'

    nkx=40; nky=40  
	nLform=5; Lcontact=2
    allocate(Lformtable(nLform)); Lformtable=(/'A1g','A1g','B1g','A1g','B2g'/)
    allocate(Lformbasis(3,nLform)); Lformbasis(:,1)=0
	Lformbasis(:,2)=(/1,0,0/); Lformbasis(:,3)=(/1,0,0/)
	Lformbasis(:,4)=(/1,1,0/); Lformbasis(:,5)=(/1,1,0/)

	call setLformfunc(ng)
	call setOformfunc()

    SkipInterOrbitPair=.false.
	guess='A1g';  call setMolecules()

    nw=400; wmax=100
    wmin=1.e-4; wir=1.e-3
    diverge=40
  
  return
end subroutine WangFa5d

subroutine WangFa5dhk(norb,kv,Q,hk,model,mf,slvbsn)
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

  bondmax=1
  do iy=-bondmax,bondmax; do ix=-bondmax,bondmax
     t=>model%t(:,:,ix,iy,1)
	 hk=hk+t*exp( one*( ix*vk(1)+iy*vk(2) ) )
  end do; end do
 
  return
  error=sum( abs( aimag(hk) ) )
  if(error>1.e-4)pause 'imag(hk)/=0?'

  return
end subroutine WangFa5dhk
