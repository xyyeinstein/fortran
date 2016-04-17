subroutine FS2F()
  use workplace
  implicit none

  real*8 Q,x,y,z,tt
  integer ix,iy,iz,a,b,error
  complex*16 ttt

  ! this is a test of supercell method, and should be compared to the 5-orbit model

  DataReady=.false.
  AppendBreakPoint=.false.
  BCSflow=.false.
  Appendkmesh=.true.; kmeshfile='FS2Fkmesh.dat'

  norb=10; natom=2
  model%norb=norb; model%natom=natom; model%sample='FS2F'
  print*,'The sample chosen is ',model%sample

  allocate(model%ra(3,natom)); model%ra(:,1)=0; model%ra(:,2)=(/0.5,0.5,0./)  !(/0.5,0.,0./); model%ra(:,2)=(/0.,0.5,0./)
  model%a=(/1,0,0/); model%b=(/0,1,0/); model%c=(/0,0,1/)

  allocate(model%orbits(norb)); model%orbits(1:5)=(/'A1g','E2g','E1g','B1g','B2g'/)
  model%orbits(6:10)=model%orbits(1:5)
  model%periodiclayers=.false.

  usegroup=.true.;  group='D2h';  ng=4   
  model%group=group; model%ng=ng

  model%U=4/2; model%JH=0.32/2; model%Uab=model%U-2*model%JH; 

  model%zerotz=.true.
  allocate(model%t(10,10,-4:4,-4:4,0:0))
  open(10,file='FS2F.model')
  model%t=0; model%mu=7.05+0.0135
  do while (.true.)
     read(10,*,iostat=error)a,b,x,y,z,tt
     ix=nint(x*2); iy=nint(y*2); iz=nint(z*2)
     model%t(a,b,ix,iy,iz)=tt
     if(error/=0)exit
  end do; close(10)
  
  Q=0; call checksymmetry(Q); 
  call FermiSurf(norb,Q,model)
  call BandStructure(norb,Q,model)
  

  model%ephcoupled=.true.
  print*,'Veph='; read*,model%Veph
  !model%Veph=2.0
  if(model%ephcoupled.and.model%Veph<1.e-6)model%ephcoupled=.false.



  nkx=40; nky=40     
  nLform=3; Lcontact=2
  allocate(Lformtable(nLform)); Lformtable=(/'A1g','A1g','B2g'/) !,'A1g','B1g'/)
  allocate(Lformbasis(3,nLform));	Lformbasis(:,1)=0
  Lformbasis(:,2)=(/1,1,0/)*0.5; Lformbasis(:,3)=(/1,1,0/)*0.5
 ! Lformbasis(:,4)=(/1,0,0/); Lformbasis(:,5)=(/0,1,0/)

  call setLformfunc(ng)
  call setOformfunc()

  SkipInterOrbitPair=.false.
  guess='A1g';  call setMolecules()

  nw=400; wmax=100
  wmin=1.e-4; wir=1.e-3
  diverge=40

  return
end subroutine FS2F

subroutine FeSe2Fehk(norb,kv,Q,hk,model,mf,slvbsn)
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
  
  p(1:2)=kv*pi;p(3)=2*Q*pi-sum(p(1:2)); p=p/2
 
  hk=0
  do iz=-1,1; if(model%zerotz.and.iz/=0)cycle; do iy=-4,4; do ix=-4,4
     t=>model%t(:,:,ix,iy,iz)
     r=(/ix,iy,iz/); phase=exp(-one*sum(p*r))
	 hk=hk+t*phase
  end do; end do; end do

  do iorb=1,norb; hk(iorb,iorb)=hk(iorb,iorb)-model%mu; end do
  do ix=1,norb; do iy=1,norb
    if(abs(hk(ix,iy)-conjg(hk(iy,ix)))>1.e-5)then
      pause 'hk not Hermitian'
    end if
  end do; end do
  
  return
end subroutine FeSe2Fehk
