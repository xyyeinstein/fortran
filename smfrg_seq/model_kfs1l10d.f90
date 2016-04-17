subroutine KFeSe1L10d()
  use workplace
  implicit none

  real*8 Q,x,y,z,tt
  integer ix,iy,iz,a,b,error

  ! this is a test of supercell method, and should be compared to the 5-orbit model

  DataReady=.false.
  AppendBreakPoint=.false.
  Appendkmesh=.true.; kmeshfile='K1L10dkmesh.dat'

    norb=10; natom=2
    model%norb=norb; model%natom=natom; model%sample='K1L10d'
    print*,'The sample chosen is ',model%sample 

	allocate(model%ra(3,natom)); model%ra(:,1)=0; model%ra(:,2)=(/1,1,0/)*0.5
    model%a=(/1,0,0/); model%b=(/0,1,0/); model%c=(/0,0,1/)

    allocate(model%orbits(norb)); model%orbits(1:5)=(/'A1g','E2g','E1g','B1g','B2g'/)
	model%orbits(6:10)=model%orbits(1:5)
	model%periodiclayers=.false.

    usegroup=.true.;  group='D4h';  ng=8    
	model%group=group; model%ng=ng

    model%U=4/2; model%JH=0.54/2; model%Uab=model%U-2*model%JH; model%mu=9.27

	model%zerotz=.true.
    allocate(model%t(10,10,-4:4,-4:4,-1:1))
    open(10,file='KFeSe10d.model')
    model%t=0
    do while (.true.)
       read(10,*,iostat=error)a,b,x,y,z,tt
       ix=nint(x*2); iy=nint(y*2); iz=nint(z*2)
       model%t(a,b,ix,iy,iz)=tt
       if(error/=0)exit
    end do; close(10)
	Q=0; call checksymmetry(Q);	!call FermiSurf(norb,Q,model)

    nkx=40; nky=40     
    nLform=3; Lcontact=2
    allocate(Lformtable(nLform)); Lformtable=(/'A1g','A1g','B2g'/)
    allocate(Lformbasis(3,nLform));	Lformbasis(:,1)=0
	Lformbasis(:,2)=(/1,1,0/)*0.5; Lformbasis(:,3)=(/1,1,0/)*0.5

	call setLformfunc(ng)
	call setOformfunc()

    SkipInterOrbitPair=.false.
	guess='A1g';  call setMolecules()

    nw=400; wmax=100
    wmin=1.e-4; wir=1.e-3
    diverge=40

  return
end subroutine KFeSe1L10d

subroutine KFeSe1L10dhk(norb,kv,Q,hk,model,mf,slvbsn)
  use standard_derived_types
  implicit none
  integer norb
  real*8 kv(2),Q
  complex*16 hk(norb,norb)
  type (modelblock), target :: model
  type (mfblock) :: mf
	logical slvbsn

  integer ix,iy,iz,r(3),iorb
  real*8 pi,p(3)
  complex*16, pointer, dimension (:,:) :: t
  complex*16 one,phase

  one=cmplx(0,1); pi=asin(1.d0)*2
  
  p(1:2)=kv*pi; p(3)=2*Q*pi-sum(p(1:2)); p=p/2
 
  hk=0
  do iz=-1,1; if(model%zerotz.and.iz/=0)cycle; do iy=-4,4; do ix=-4,4
     t=>model%t(:,:,ix,iy,iz)
     r=(/ix,iy,iz/); phase=exp(-one*sum(p*r))
	 hk=hk+t*phase
  end do; end do; end do

  do iorb=1,norb; hk(iorb,iorb)=hk(iorb,iorb)-model%mu; end do
  
  return
end subroutine KFeSe1L10dhk
