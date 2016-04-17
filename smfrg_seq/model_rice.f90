subroutine rice()
  use workplace
  implicit none

  real*8 Q,x,y,z,tt
  integer ix,iy,iz,a,b,nband

  BCSflow=.false.
  DataReady=.false.
  AppendBreakPoint=.false.
  SkipInterOrbitPair=.false.
  skipOddOrbPair=.false.
  Appendkmesh=.false.; kmeshfile='ricekmesh.dat'
  useprojectors=.false.;  Appendprojectors=.false.; projectorfile='ricePrj.dat'

    norb=2; natom=1
    model%norb=norb; model%natom=natom; model%sample='rice'
    print*,'The sample chosen is ',model%sample 

	allocate(model%ra(3,natom)); model%ra(:,1)=0
    model%a=(/1,0,0/); model%b=(/0,1,0/); model%c=(/0,0,1/)

    allocate(model%orbits(norb)); model%orbits=(/'B1g','A1g'/)
	model%periodiclayers=.false.

    usegroup=.true.;  group='C2v';  ng=4
	model%group='C2v'; model%ng=4

    model%U=4; model%JH=0; model%Uab=0; model%Vnn=0.4; model%mu=0 
	model%nest=8
	model%qnest(:,1)=(/0.25,0./); model%qnest(:,2)=(/0.75,1./)
	model%qnest(:,3)=(/0.,0.25/); model%qnest(:,4)=(/1.,0.75/)
	model%qnest(:,5)=(/0.5,0./); model%qnest(:,6)=(/0.5,1./)
	model%qnest(:,5)=(/0.,0.5/); model%qnest(:,6)=(/1.,0.5/)

	Q=0; call checksymmetry(Q);	
	call FermiSurf(norb,Q,model); 
	call fermisurf_rice(norb,model);
	call bandstructure(norb,Q,model) 

    nkx=40; nky=40     
    nLform=9; Lcontact=2
    !nLform=1; Lcontact=0
	allocate(Lformtable(nLform)); !Lformtable=(/'A1g'/) !,'A1g','B1g','A1g','B2g'/)
    allocate(Lformbasis(3,nLform))
	Lformbasis(:,1)=0;  Lformtable(1)='A1g'
	!Lformbasis(:,2)=(/1,0,0/); Lformbasis(:,3)=(/1,0,0/)
	!Lformbasis(:,4)=(/1,1,0/); Lformbasis(:,5)=(/1,1,0/)
	Lformbasis(:,2)=(/1,0,0/); Lformtable(2)='A1g'
	Lformbasis(:,3)=(/1,0,0/); Lformtable(3)='E1u'
	Lformbasis(:,4)=(/0,1,0/); Lformtable(4)='A1g'
	Lformbasis(:,5)=(/0,1,0/); Lformtable(5)='E2u'
	Lformbasis(:,6)=(/1,1,0/); Lformtable(6)='A1g'
	Lformbasis(:,7)=(/1,1,0/); Lformtable(7)='B2g'
	Lformbasis(:,8)=(/1,1,0/); Lformtable(8)='E1u'
	Lformbasis(:,9)=(/1,1,0/); Lformtable(9)='E2u'

	call setLformfunc(ng)
	call setOformfunc()

	guess='***';  call setMolecules()

	!use higher symmetry point group to speed up the calculation
    group='C4v';  ng=8   
	model%group='C4v'; model%ng=8

    nw=400; wmax=100
    wmin=1.e-4; wir=1.e-4
    diverge=40

  return
end subroutine rice


subroutine Ricehk(norb,kv,hk,model,mf,slvbsn)
  use standard_derived_types
  implicit none

  integer norb
  real*8 kv(2),z
  complex*16 hk(norb,norb)
  type (modelblock) :: model
  type (mfblock) :: mf
	logical slvbsn

  real*8 pi,ck(2),t1,t2,mu,ek,xk,dk,spingap
  real*8 normanhk

  pi=asin(1.d0)*2
  
  ck=cos(kv*pi); spingap=0.1; z=1
  
  ek=normanhk(kv,model%mu)/z
  xk=-2*sum(ck)
  dk=spingap*0.5*(ck(1)-ck(2))

  hk=0
  hk(1,1)=ek; hk(2,2)=-xk; hk(1,2)=dk; hk(2,1)=dk

  return
end subroutine Ricehk

function normanhk(kv,mu)
  implicit none
  real*8 kv(2),normanhk,mu

  real*8 pi,ck(2),ck2(2),t1,t2,t3,t4,t5,mueff,ek

  pi=asin(1.d0)*2
  
  t1=-0.5911; t2=0.1636; t3=-0.0519; t4=-0.1117; t5=0.051; mueff=-0.12+mu*abs(t1)
  ck=cos(kv*pi); ck2=cos(kv*pi*2)
  
  ek=t1*sum(ck)/2+t2*ck(1)*ck(2)+t3*sum(ck2)/2+t4*(ck(1)*ck2(2)+ck(2)*ck2(1))/2+t5*ck2(1)*ck2(2)-mueff
  ek=ek*4./0.5911

  normanhk=ek

  return
end function normanhk

subroutine FermiSurf_rice(norb,model)
  use standard_derived_types
  integer norb
  type (modelblock) :: model

  integer ix,iy
  real*8 pi,kv(2)
  complex*16 hk(norb,norb),gk(norb,norb),one
  real*8 eta,dos

  pi=asin(1.d0)*2; one=cmplx(0,1); eta=0.1
  factor=1
 
  open(10,file='Akdos.dat')
  do iy=-40,40;  do ix=-40,40; kv=(/ix,iy/)*factor/40
     call ricehk(norb,kv,hk,model)
	 gk=-hk; gk(1,1)=gk(1,1)+eta*one; gk(2,2)=gk(2,2)+eta*one
	 call ZINVERT(norb,gk)
	 dos=-aimag( gk(1,1) )/pi
	 write(10,*)dos
  end do; end do
  close(10)
  print*,'Ak at the fermi level output to file'
  return
end subroutine FermiSurf_rice
