subroutine CuO2()
  use workplace
  implicit none

  real*8 Q
  
  useX0=.true.; appendX0=.true.
  BCSflow=.false.
  DataReady=.false.
  QuickSearch=.false.
  SkipInterOrbitPair=.false.
  AppendBreakPoint=.false.
  Appendkmesh=.false.; kmeshfile='CuO2kmesh.dat'
  useprojectors=.false.;  Appendprojectors=.true.; projectorfile='CuO2Prj.dat'
  useprojectors=useprojectors.and.(.not.QuickSearch)

  norb=3; natom=3
  model%norb=3; model%natom=3; model%sample='CuO2'
  print*,'The sample chosen is ',model%sample

  
  allocate(model%ra(3,3),model%orbits(3))
  model%orbits=(/'A1g','E1u','E2u'/)
  model%ra(:,1)=(/0.,0.,0./); model%ra(:,2)=(/0.5,0.,0./)
  model%ra(:,3)=(/0.,0.5,0./)
  model%a=(/1,0,0/); model%b=(/0,1,0/); model%c=(/0,0,1/)
  model%propergauge=.false.

  model%Udd=10; model%Upp=1; model%Vpp=0.2; model%Vdd=1; model%Vpd=3; model%mu=-2
  call bandstructure(model%norb,Q,model); 

  model%nband=1; allocate(model%nativeband(1),model%npatch(1),model%pocket(1),model%mushift(1))
  model%nativeband=3; model%npatch=96; model%pocket='G_'; model%mushift=0
  call meshfermisurface(model); call varma()

    usegroup=.true.  
    group='C2v'; ng=4
    model%group='C4v'; model%ng=8
    !Q=0; call checksymmetry(Q)
	!call FermiSurf(norb,Q,model); !pause 'continue?'

    nkx=40; nky=40  
	nLform=13; Lcontact=2
	!nLform=1; Lcontact=0
    allocate(Lformtable(nLform),Lformbasis(3,nLform))
    Lformbasis(:,1)=0; Lformtable(1)='A1g'
	Lformbasis(:,2)=(/1,0,0/); Lformbasis(:,3)=(/1,0,0/); Lformtable(2:3)=(/'A1g','E1u'/)
	Lformbasis(:,4)=(/0,1,0/); Lformbasis(:,5)=(/0,1,0/); Lformtable(4:5)=(/'A1g','E2u'/)
	Lformbasis(:,6)=(/0.5,0.,0./); Lformbasis(:,7)=(/0.5,0.,0./); Lformtable(6:7)=(/'A1g','E1u'/)
	Lformbasis(:,8)=(/0.,0.5,0./); Lformbasis(:,9)=(/0.,0.5,0./); Lformtable(8:9)=(/'A1g','E2u'/)
	Lformbasis(:,10)=(/0.5,0.5,0./); Lformbasis(:,11)=(/0.5,0.5,0./); Lformtable(10:11)=(/'A1g','B2g'/)
	Lformbasis(:,12)=(/0.5,0.5,0./); Lformbasis(:,13)=(/0.5,0.5,0./); Lformtable(12:13)=(/'E1u','E2u'/)

	call setLformfunc(ng)
	call setOformfunc()

    SkipInterOrbitPair=.false.
	guess='***';  call setMolecules()

    nw=400; wmax=100
    wmin=1.e-4; wir=2.5e-4
    diverge=80

    if(useX0.and.appendX0)then
	  print*,'Udd,Upp,Vdd,Vpp,Vpd='
	  read*,model%Udd,model%Upp,model%Vdd,model%Vpp,model%Vpd
	end if

  return
end subroutine CuO2

subroutine CuO2Hk(norb,vk,hk,model,mf,slvbsn)
  use standard_derived_types
  implicit none

  integer norb
  real*8 vk(2)
  complex*16 hk(norb,norb)
  type(modelblock) :: model
	type(mfblock) :: mf
	logical slvbsn

  real*8 tpd,tpp,tdd,ed,ep,mu
  real*8 ckx,cky,sx,sy,cxy1,cxy2

  complex*16 one
  real*8 pi,dr(2)
  integer i,j

  one=cmplx(0,1); pi=asin(1.d0)*2
  !tpd=1 ; tpp=0.8 ; tdd=0 ; ed=0 ; ep=-3 ; mu=1.2
  tpd=1.3 ; tpp=0.6 ; tdd=0. ; ed=0 ; ep=0 ; mu=4.71+model%mu
  !tpd=0. ; tpp=0. ; tdd=1. ; ed=0 ; ep=1 ; mu=0

  pi = asin(1.d0)*2
  one = cmplx(0,1)

  ckx=cos(vk(1)*pi) ;  cky=cos(vk(2)*pi) 
  sx=sin(pi*vk(1)/2) ; sy=sin(pi*vk(2)/2)
  cxy1=cos((vk(1)+vk(2))*pi/2) ; cxy2=cos((vk(1)-vk(2))*pi/2)

  hk(1,1)=-2*tdd*(ckx+cky)+ed-mu
  hk(2,2)=ep-mu
  hk(3,3)=ep-mu
  hk(1,2)=2*one*tpd*sx ; hk(2,1)=conjg(hk(1,2))
  hk(1,3)=-2*one*tpd*sy ; hk(3,1)=conjg(hk(1,3))
  hk(2,3)=2*tpp*(cxy1-cxy2) ; hk(3,2)=conjg(hk(2,3))
  
  if(model%propergauge)then
    do i=1,3; do j=1,3; dr=model%ra(1:2,j)-model%ra(1:2,i)
       hk(i,j)=hk(i,j)*exp(-one*sum(vk*dr)*pi )
    end do; end do
  end if

  return
end subroutine CuO2Hk

subroutine varma()
  use workplace
  implicit none

  integer ik,ib
  real*8 k(2),eval(norb),fk,gk,Q
  complex*16 hk(norb,norb)

  open(10,file='varmafk.dat')
  Q=0
  do ik=1,model%nkf
     k=model%kf(:,ik)
	 call gethk(norb,k,Q,hk,model,MF,.false.)
	 call zheigen(norb,hk,eval)
	 ib=model%band(ik)
	 fk = cos(k(1)*pi/2)*aimag( conjg(hk(1,ib))*hk(2,ib) ) 
	 gk = cos(k(2)*pi/2)*aimag( conjg(hk(1,ib))*hk(3,ib) ) 
	 write(10,*)fk+gk,fk-gk
  end do

  return
end subroutine varma
