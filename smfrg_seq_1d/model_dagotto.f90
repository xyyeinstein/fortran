subroutine Dagotto()
  use workplace
  implicit none

  real*8 Q
  integer nband

  usex0=.true.
  appendx0=.true.
  appendkmesh=.true.
  appendfs=.true.
  DataReady=.false.
  uniformqmesh=.true.
  AppendBreakPoint=.false.
  kmeshfile='Dagottokmesh.dat'

  norb=3; natom=1; model%n4logqmesh=2
  model%norb=norb; model%natom=natom; model%sample='Dagott'
  print*,'The sample chosen is ',model%sample

  allocate(model%ra(3,natom)); model%ra=0
  allocate(model%orbits(norb)); model%orbits=(/'E2u','E1u','B2g'/)
  model%periodiclayers=.false.

  usegroup=.true.;  group='D4h';  ng=8
  model%group=group; model%ng=ng

  open(10,file='frg.input')
  read(10,*)model%U,model%JH,model%mu
  close(10)
  !model%U=8; model%JH=3
  model%Uab=model%U-2*model%JH

  nband=3; model%nband=nband
  allocate(model%npatch(nband),model%nativeband(nband),model%pocket(nband),model%mushift(nband))
  model%nativeband=(/1,2,3/); 
  model%pocket=(/'G_','G_','XY'/); 
  model%npatch=(/64,64,128/); 
  model%mushift=(/0.d0,0.d0,0.d0/)
  Q=0; call checksymmetry(Q)
  call bandstructure(norb,Q,model)
  call fermisurf(norb,Q,model)
  call meshfermisurface(model)
  

  nkx=40; nky=40    
  nLform=5; Lcontact=2
  allocate(Lformtable(nLform)); Lformtable=(/'A1g','A1g','B2g','A1g','B1g'/)
  allocate(Lformbasis(3,nLform)); Lformbasis(:,1)=0
  Lformbasis(:,2)=(/1,1,0/); Lformbasis(:,3)=(/1,1,0/)
  Lformbasis(:,4)=(/1,0,0/); Lformbasis(:,5)=(/1,0,0/)

  call setLformfunc(ng)
  call setOformfunc()

  SkipInterOrbitPair=.false.
	guess='***';  call setMolecules()

  nw=400; wmax=100
  wmin=1.e-4; wir=1.e-5
  diverge=80
  return
end subroutine Dagotto

subroutine Dagottohk(norb,kv,Q,hk,model,mf,slvbsn)
  use standard_derived_types
  implicit none
  integer norb
  real*8 kv(2),Q
  complex*16 hk(3,3)
  type (modelblock) :: model
  type (mfblock) :: mf
	logical slvbsn

  real*8 t1,t2,t3,t4,t5,t6,t7,t8,delta,mu
  real*8 ck(2),sk(2),vk(2)
  complex*16 one

  t1=0.2; t2=0.6; t3=0.3; t4=-0.1; t5=2; t6=3; t7=-2; t8=1
  delta=4; one=cmplx(0,1); mu=model%mu

  vk=kv*3.1415926; ck=cos(vk); sk=sin(vk)

  hk(1,1)=2*t2*ck(1)+2*t1*ck(2)+4*t3*ck(1)*ck(2)-mu
  hk(2,2)=2*t1*ck(1)+2*t2*ck(2)+4*t3*ck(1)*ck(2)-mu
  hk(3,3)=2*t5*(sum(ck))+4*t6*ck(1)*ck(2)-mu+delta
  hk(1,2)=4*t4*sk(1)*sk(2); hk(2,1)=hk(1,2)
  hk(1,3)=2*one*t7*sk(1)+4*one*t8*sk(1)*ck(2); hk(3,1)=conjg( hk(1,3) )
  hk(2,3)=2*one*t7*sk(2)+4*one*t8*sk(2)*ck(1); hk(3,2)=conjg( hk(2,3) )

  return
end subroutine Dagottohk
