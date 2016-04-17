subroutine Raghu()
  use workplace
  implicit none

  integer nband
  real*8 Q
  
  useX0=.true.
  appendX0=.false.
  Appendkmesh=.false.
  appendfs=.false.
  if(.true.)then
    appendx0=.true.
    appendkmesh=.true.
    appendfs=.true.
  end if
  DataReady=.false.
  uniformkmesh=.false.
  AppendBreakPoint=.false.
  kmeshfile='Raghukmesh.dat'

  norb=2; natom=1; nambu=4
  model%norb=norb; model%natom=natom; model%sample='Raghu'
  print*,'The sample chosen is ',model%sample 

	allocate(model%ra(3,natom));  model%ra=0
  allocate(model%orbits(norb)); model%orbits=(/'E1u','E2u'/)
	model%periodiclayers=.false.
	
	usegroup=.true.;  group='D4h';  ng=8
  model%group=group; model%ng=ng

  !open(10,file='frg.input')
  model%U=1.0
  model%JH=0.2
  model%mu=1.55
  model%Uab=model%U-2*model%JH
  nband=2; model%nband=nband
  allocate(model%npatch(nband),model%nativeband(nband),model%pocket(nband),model%mushift(nband))
  model%nativeband=(/1,2/);	model%pocket=(/'GM','XY'/); model%npatch=120; model%mushift=0.

  Q=0; call checksymmetry(Q)
  call fermisurf(norb,Q,model)
  call bandstructure(norb,Q,model) 
  stop
  call meshfermisurface(model)
  nkx=400; nky=400  !used presently only to calculate doping level  
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
    diverge=100

  return
end subroutine Raghu

subroutine Raghuhk(norb,kv,Q,hk,model,mf,slvbsn)
  use standard_derived_types
  implicit none
  integer norb
  real*8 kv(2),Q
  complex*16 hk(norb,norb)
  type (modelblock) :: model
  type (mfblock) :: mf
	logical slvbsn

  real*8 t1,t2,t3,t4,ek,xk,zk
  real*8 ck(2),sk(2),vk(2)

  t1=-1.; t2=1.3; t3=-0.85; t4=-0.85

  t1=t1*1.
  t2=t2*1.
  t3=t3*1.
  t4=t4*1.

  vk=kv*3.1415926;  ck=cos(vk); sk=sin(vk)
  ek=-(t1+t2)*sum(ck)-4*t3*ck(1)*ck(2)-model%mu
  xk=-4*t4*sk(1)*sk(2)
  zk=-(t1-t2)*(ck(1)-ck(2))
  hk(1,1)=ek+zk; hk(2,2)=ek-zk; hk(1,2)=xk; hk(2,1)=xk

  return
end subroutine Raghuhk
