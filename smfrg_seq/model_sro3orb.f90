subroutine SrRuO3orb()
  use workplace
  implicit none

  real*8 Q,cutoff
  integer iform,nband
  logical dorpa

  dorpa=.true.
  !useX0=.true.; appendX0=.false.
  QuickSearch=.false.
  DataReady=.false.
  AppendBreakPoint=.false.
  BCSFlow=.false.
  Appendkmesh=.true.; kmeshfile='SRO3orbkmesh.dat'
  useprojectors=.false.

    norb=3; natom=1
    model%norb=norb; model%natom=natom; model%sample='SRO3orb'
    print*,'The sample chosen is ',model%sample 

	allocate(model%ra(3,natom));  model%ra=0
    allocate(model%orbits(norb)); model%orbits=(/'E1g','E2g','B2g'/)
	model%periodiclayers=.false.
	
	usegroup=.true.;  group='C4v';  ng=8
    model%group=group; model%ng=ng

    model%U=3; model%JH=0.5; model%Uab=model%U-2*model%JH; model%mu=1.1 !to be supplemented case: mu=1
	nband=3; model%nband=nband
	allocate(model%npatch(nband),model%nativeband(nband),model%pocket(nband),model%mushift(nband))
	model%nativeband=(/1,2,3/); model%pocket=(/'G_','G_','G_'/); model%npatch=(/96,96,96/); model%mushift=0.


    !---------------------------------tune interactions interactively----------------------
	if(dorpa.or.(useX0.and.appendX0))then
	  print*,'you are now free to tune interactions:'
	  print*,'U='; read*,model%U
	  print*,'JH='; read*,model%JH; model%Uab=model%U-2*model%JH
	  print*,'Vnn='; read*,model%Vnn
	end if
    !---------------------------------------------------------------------------------------

	Q=0; call checksymmetry(Q);	
	call FermiSurf(norb,Q,model); call bandstructure(norb,Q,model)
	call meshfermisurface(model)

    nkx=40; nky=40     
    nLform=9; Lcontact=2
	allocate(Lformbasis(3,nLform),Lformtable(nLform))

    iform=0
	call addform( (/0,0,0/)*1.d0,'A1g')
	call addform( (/1,0,0/)*1.d0,'A1g')
	call addform( (/-1,0,0/)*1.d0,'A1g')
	call addform( (/0,1,0/)*1.d0,'A1g')
	call addform( (/0,-1,0/)*1.d0,'A1g')
	call addform( (/1,1,0/)*1.d0,'A1g')
	call addform( (/1,-1,0/)*1.d0,'A1g')
	call addform( (/-1,1,0/)*1.d0,'A1g')
	call addform( (/-1,-1,0/)*1.d0,'A1g')
	!call addform( (/2,0,0/)*1.d0,'A1g')
	!call addform( (/0,2,0/)*1.d0,'A1g')
	!call addform( (/2,0,0/)*1.d0,'E1u')
	!call addform( (/0,2,0/)*1.d0,'E2u')
	if(iform/=nLform)stop 'nform error @ SrRuO1orb'


	call setLformfunc(1)   !use identity group (by setting ng=1) to setup fundamental lattice forms

	group='C2v'; ng=4; model%group=group; model%ng=ng
	call setOformfunc()
	guess='***';  call setMolecules()

	group='C4v'; ng=8; model%group=group; model%ng=ng

    nw=400; wmax=1.e2
    wmin=1.e-4; wir=1.e-4
    diverge=40

    if(.not.dorpa)return

	call writeinputs()
	call setupworkplace()
	cutoff=1.e-2; call rpa(cutoff,0,'SRO3dX0.dat'); stop

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
end subroutine SrRuO3orb


subroutine SrRuO3orbhk(norb,kv,hk,model,mf,slvbsn)
  use standard_derived_types
  implicit none
  integer norb
  real*8 kv(2)
  complex*16 hk(norb,norb)
  type (modelblock) :: model
  type (mfblock) :: mf
	logical slvbsn

  real*8 t1,t2,t3,t4,gap,mu
  real*8 ck(2),sk(2),vk(2)

  t1=1; t2=0.1; t3=0.8; t4=0.35; gap=-0.2; mu=model%mu

  vk=kv*3.1415926;  ck=cos(vk); sk=sin(vk)

  hk=0

  hk(1,1)=-2*t1*ck(1)-mu; hk(2,2)=-2*t1*ck(2)-mu
  hk(1,2)=-4*t2*sk(1)*sk(2); hk(2,1)=hk(1,2)
  hk(3,3)=-2*t3*sum(ck)-4*t4*ck(1)*ck(2)-mu+gap

  return
end subroutine SrRuO3orbhk
