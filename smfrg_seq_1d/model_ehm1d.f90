subroutine ehm1d()
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
  uniformqmesh=.true.
  appendkmesh=.false.
  appendfs=.false.
  QuickSearch=.false.
  outputXq=.false.
  AppendBreakPoint=.false.
  SkipInterOrbitPair=.false.
  kmeshfile='ehm1dkmesh.dat'
  useprojectors=.false.; Appendprojectors=.false.; projectorfile='ehm1dPrj.dat'
  useprojectors=useprojectors.and.(.not.QuickSearch)   !projectors are not used in the QuickSearch mode.

  norb=1; natom=1; nambu=2
  model%norb=norb; model%natom=natom; model%sample='ehm1d'
  print*,'The sample chosen is ',model%sample 

  model%nest=4
  model%qnest(:,1)=(/-1.d0,0.d0/)
  model%qnest(:,2)=(/1.d0,0.d0/)
  model%qnest(:,3)=(/-0.5d0,0.d0/)
  model%qnest(:,4)=(/0.5d0,0.d0/)
  allocate(model%ra(3,natom)); model%ra=0
  allocate(model%orbits(norb)); model%orbits='A1g'
  model%a=(/1,0,0/); model%b=(/0,1,0/); model%c=(/0,0,1/)
    
  usegroup=.true.;  group='C2_';  ng=2
  model%group=group; model%ng=ng
    
  open(10,file='U.input'); read(10,*)model%U; close(10)
  open(10,file='Vnn.input'); read(10,*)model%Vnn; close(10)
  open(10,file='mu.input'); read(10,*)model%mu; close(10)
  
      
  Q=0
  call bandstructure(norb,Q,model)  
  nband=1
  model%nband=nband
  allocate(model%npatch(nband),model%nativeband(nband),model%pocket(nband),model%mushift(nband))
  model%nativeband=1
  model%pocket='G_'
  model%npatch=1
  model%mushift=0.
  !call meshfermisurface(model)
  !print*,'meshfermisurface OK!'
  

  !============================Set Form Factor========================================
  nkx=600
  nky=nkx  
  nLform=5; Lcontact=4
  allocate(Lformtable(nLform))
  allocate(Lformbasis(3,nLform))
  iform=0
  call addform((/0,0,0/)*0d0,'A1g')
  call addform((/1,0,0/)*1.d0,'A1g')
  call addform((/-1,0,0/)*1.d0,'A1g')
  call addform((/2,0,0/)*1.d0,'A1g')
  call addform((/-2,0,0/)*1.d0,'A1g')
  call setLformfunc(1)
  print*,'Lform ok'
  call setOformfunc()
  print*,'Oform OK'
  guess='***'
  call setMolecules()
  print*,'Mform ok'
  !===================================================================================
	
  nw=800; wmax=100
  wmin=1.e-6; wir=1.e-5
  diverge=40

  if(dorpa)then
	  call writeinputs()
		call setupworkplace()
		cutoff=1.e-3; call rpa(cutoff,0,'lhqX0.dat'); stop
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

end subroutine

subroutine ehm1dhk(norb,kv,Q,hk,model,mf,slvbsn)
  use standard_derived_types
  implicit none
  integer norb
  real*8 kv(2),Q
  complex*16 hk(norb,norb)
  type (modelblock) :: model
  type (mfblock) :: mf
  logical slvbsn

  real*8 pi,ck(2),c2k(2),t1,t2,mu,delt1,QM(2),kvv(2)

  pi=asin(1.d0)*2
  kvv=kv
  hk=0

  t1=1
  t2=-0.0
  mu=model%mu
  ck=cos(kvv*pi)
  c2k=cos(2*kvv*pi)
  hk=-2*t1*ck(1)-2*t2*c2k(1)
  hk=hk-mu

  return
end subroutine ehm1dhk

