subroutine lhq()
  use workplace
  implicit none

  real*8 Q,cutoff,nhole
  integer iform,nband
  logical dorpa

  dorpa=.false.
  Appendkmesh=.false.
  bcsflow=.false.
  useX0=.false.
  appendX0=.false.
  uniformqmesh=.false.
  DataReady=.false.
  QuickSearch=.false.
  AppendBreakPoint=.false.
  kmeshfile='lhqkmesh.dat'
  useprojectors=.false.; appendprojectors=.false.; projectorfile='lhqdproj.dat'
  useprojectors=useprojectors.and.(.not.QuickSearch)

  !open(10,file='ctrl.input')
  !  read(10,*)Appendkmesh
  !  read(10,*)AppendX0
  !close(10)
  !print*,appendkmesh,appendx0; stop
  norb=2; natom=2; model%n4logqmesh=2
	model%norb=norb; model%natom=natom; model%sample='lhq'
  print*,'The sample chosen is ',model%sample 

	allocate(model%ra(3,natom)); model%ra(:,1)=(/0.5,0.,0./); model%ra(:,2)=(/0.,0.5,0./)
	model%a=(/1,0,0/); model%b=(/0,1,0/); model%c=(/0,0,1/)
  allocate(model%orbits(norb)); model%orbits='A1g'
	model%periodiclayers=.false.
    
	usegroup=.true.;  group='C4v';  ng=8
	model%group=group; model%ng=ng
    
	model%propergauge=.false.

  if(.true.)then
	  !for KFe2Se2
	  nband=1; model%nband=nband; model%mu=0.0 
	  allocate(model%npatch(nband),model%nativeband(nband),model%pocket(nband),model%mushift(nband))
	  model%nativeband=2;	model%pocket='M_'; model%npatch=120; model%mushift=0.
  else
    !for KFe2As2
    nband=2; model%nband=nband; 
	  allocate(model%npatch(nband),model%nativeband(nband),model%pocket(nband),model%mushift(nband))
	  model%nativeband=(/1,2/);	model%pocket=(/'G_','M_'/); model%npatch=120; model%mushift=0.
  end if
  
  if(.false.)then
    open(10,file='frg.input')
    read(10,*)model%t1
    read(10,*)model%t2
    read(10,*)model%t3
    read(10,*)model%filling
    read(10,*)model%U
    read(10,*)model%Vnn
    close(10)
  end if

  model%U=3.0
  model%filling=0
  model%filling=-model%filling+1
  print*,model%filling
  call filling2mu(.false.)
  !model%mu=-0.30
  Q=0; call checksymmetry(Q); call fermisurf(norb,Q,model)
	call bandstructure(norb,Q,model) 
	call meshfermisurface(model)


  nkx=100
	nky=100     
    	
	if(dorpa)then
	  nLform=1
	  Lcontact=0
	else
	  nLform=5
	  Lcontact=1
	end if
	allocate(Lformbasis(3,nLform),Lformtable(nLform))

  iform=0
	call addform( (/0,0,0/)*1.d0,'A1g')
	call addform( (/0.5,0.5,0./)*1.d0,'A1g')
	!call addform( (/0.5,0.5,0./)*1.d0,'B2g')
	call addform( (/-0.5,0.5,0./)*1.d0,'A1g')
	call addform( (/-0.5,-0.5,0./)*1.d0,'A1g')
	call addform( (/0.5,-0.5,0./)*1.d0,'A1g')
	call addform( (/1,0,0/)*1.d0,'A1g')
	!call addform( (/1,0,0/)*1.d0,'B1g')
	call addform( (/-1,0,0/)*1.d0,'A1g')
  call addform( (/0,1,0/)*1.d0,'A1g')
	call addform( (/0,-1,0/)*1.d0,'A1g')

	call setLformfunc(1)  !set fundamental lattice forms instead of irreducible harmonics
	call setOformfunc()

	guess='***';  call setMolecules()

  nw=400
	wmax=1.e2
  wmin=1.e-4
	wir=1.e-5
  wirx0=1.e-5
  diverge=200

	if(.not.dorpa)return

	call writeinputs()
	call setupworkplace()
	cutoff=1.e-3; call rpa(cutoff,appendkmesh,'lhqX0.dat'); stop

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
end subroutine lhq

subroutine lhqhk(norb,kv,hk,model,mf,slvbsn)
  use standard_derived_types
  implicit none
  integer norb
  real*8 kv(2),vk(2)
  complex*16 hk(2,2),one
  real*8 pi,t1,t2,t3,mu
  type (modelblock) :: model
  type (mfblock) :: mf
  logical slvbsn

  one=cmplx(0d0,1d0)
  pi=asin(1.d0)*2; vk=kv*pi
  !t1=-model%t1
  !t2=-model%t2
  !t3=-model%t3

  t1=-.3; t2=-1.2; t3=0.8       !for KFe2Se2
  !t1=-0.3; t2=-1.4; t3=0.6       !for BaFe2As2
  mu=model%mu

  !intracell hopping between atom 1 and 2
  hk(1,2)=-4*t1*cos(vk(1)*0.5)*cos(vk(2)*0.5)
  hk(2,1)=hk(1,2)  

  !intercell hopping between like atoms
  hk(1,1)=-2*t2*cos(vk(1))-2*t3*cos(vk(2))-mu
  hk(2,2)=-2*t3*cos(vk(1))-2*t2*cos(vk(2))-mu
  if(model%propergauge)then
    hk(1,2)=hk(1,2)*exp(-one*sum(vk*(model%ra(1:2,1)-model%ra(1:2,2))))
    hk(2,1)=conjg(hk(1,2))
  end if
  
  return
end subroutine lhqhk

