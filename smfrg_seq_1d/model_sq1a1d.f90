subroutine square1atom1orb()
  use workplace
  implicit none

  real*8 Q,cutoff
  integer nband,iform
  logical dorpa
  usex0=.false. 
  dorpa=.false.
  doMF=.false.
  dataready=.false.
  uniformkmesh=.false.
  appendkmesh=.false.
  appendfs=.false.
  QuickSearch=.false.
  outputXq=.false.
  AppendBreakPoint=.false.
  SkipInterOrbitPair=.false.
  kmeshfile='sq11kmesh.dat'
  useprojectors=.false.; Appendprojectors=.false.; projectorfile='sq1a1dPrj.dat'
  useprojectors=useprojectors.and.(.not.QuickSearch)   !projectors are not used in the QuickSearch mode.

  norb=1; natom=1; nambu=2
  model%norb=norb; model%natom=natom; model%sample='sq1a1d'
  print*,'The sample chosen is ',model%sample 

  allocate(model%ra(3,natom)); model%ra=0
  allocate(model%orbits(norb)); model%orbits='A1g'
  model%a=(/1,0,0/); model%b=(/0,1,0/); model%c=(/0,0,1/)
  !model%ladder=.true.; model%nleg=8
    
  usegroup=.true.;  group='C4v';  ng=8
  model%group=group; model%ng=ng
    
  model%U=2.5
  model%Vnn=0.0
  model%Jnn=0
  model%mu=-0.  
  model%nest=0; 
  !model%qnest(:,1)=(/1,0/)
  !model%qnest(:,1)=(/0.25,0./); model%qnest(:,2)=(/0.75,1./)
  !model%qnest(:,3)=(/0.,0.25/); model%qnest(:,4)=(/1.,0.75/)
  !model%qnest(:,5)=(/0.5,0./); model%qnest(:,6)=(/0.5,1./)
  !model%qnest(:,7)=(/0.,0.5/); model%qnest(:,8)=(/1.,0.5/)
  !model%qnest(:,9)=(/1,1/)
    
  Q=0
  call checksymmetry(Q)
  call FermiSurf(norb,Q,model)
  call bandstructure(norb,Q,model)  
  nband=1
  model%nband=nband
  allocate(model%npatch(nband),model%nativeband(nband),model%pocket(nband),model%mushift(nband))
  model%nativeband=1
  model%pocket='G_'
  model%npatch=96
  model%mushift=0.
  call meshfermisurface(model)
  print*,'meshfermisurface OK!'
  

!============================Set Form Factor========================================
  nkx=600
  nky=nkx  
  nLform=3; Lcontact=2
  allocate(Lformtable(nLform))
  allocate(Lformbasis(3,nLform))
  iform=0
  call addform((/0,0,0/)*0d0,'A1g')
  call addform((/1,0,0/)*1.d0,'A1g')
  call addform((/1,0,0/)*1.d0,'B1g')
  call setLformfunc(ng);print*,'Lform ok'
  call setOformfunc();print*,'Oform OK'
  guess='***';  
  call setMolecules()
  print*,'Mform ok'
  !===================================================================================
	
  nw=400; wmax=100
  wmin=1.e-4; wir=1.e-4
  diverge=40

  if(dorpa)then
    call writeinputs()
    call setupworkplace()
    cutoff=1.e-3
    call rpa(cutoff,0,'lhqX0.dat')
    stop
  end if

  return
  contains
    subroutine addform(basis,form)
      real*8 basis(3)
      character*3 form
      if(iform+1>nLform)return
      iform=iform+1
      Lformbasis(:,iform)=basis
      Lformtable(iform)=form
    end subroutine addform

end subroutine square1atom1orb

subroutine square1a1dhk(norb,kv,Q,hk,model,mf,slvbsn)
  use standard_derived_types
  implicit none
  integer norb
  real*8 kv(2),Q
  complex*16 hk(norb,norb)
  type (modelblock) :: model
  type (mfblock) :: mf
logical slvbsn

  real*8 pi,ck(2),t1,t2,mu,delt1,QM(2),kvv(2)

  !call square1a1dhk_(norb,kv,Q,hk,model); return

  pi=asin(1.d0)*2; delt1=0
  QM=(/1,1/)*1.d0
  kvv=kv !+QM
    
  t1=1; t2=-0.3; mu=4*t2+model%mu
  ck=cos(kvv*pi); hk=-2*t1*ck(1)-2*t1*ck(2)-4*t2*ck(1)*ck(2)

  if(slvbsn)then
	  hk=hk*mf%nhole
	  mf%J(1,1,1)=0
	  mf%J(1,1,2:5)=4*t1*t1/model%U
	  mf%J(1,1,6:9)=0 !4*t2*t2/model%U
  end if
	hk=hk-mu

  return
end subroutine square1a1dhk

subroutine square1a1dhk_(norb,kv,Q,hk,model)
  use standard_derived_types
  implicit none
  integer norb
  real*8 kv(2),Q
  complex*16 hk(norb,norb)
  type (modelblock) :: model

  real*8 pi,ck(2),ck2(2),t1,t2,t3,t4,t5,mu,delt1

  pi=asin(1.d0)*2; delt1=0.
  
  t1=-0.5911; t2=0.1636; t3=-0.0519; t4=-0.1117; t5=0.051; mu=-0.12+model%mu
  ck=cos(kv*pi); ck2=cos(kv*pi*2)
  
  hk=( t1*ck(1)+t1*(1-delt1)*ck(2) )/2+t2*ck(1)*ck(2)+t3*sum(ck2)/2+t4*(ck(1)*ck2(2)+ck(2)*ck2(1))/2+t5*ck2(1)*ck2(2)-mu
  hk=hk*4./0.5911

  return
end subroutine square1a1dhk_
