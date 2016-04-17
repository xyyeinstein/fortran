subroutine SrRuO1orb()
  use workplace
  implicit none

  real*8 Q,cutoff
  integer iform,nband
  logical dorpa

  dorpa=.false.
  useX0=.true.; appendX0=.true.
  QuickSearch=.false.
  DataReady=.false.
  AppendBreakPoint=.false.
  BCSFlow=.false.
  Appendkmesh=.true.; kmeshfile='SRO1orbkmesh.dat'
  useprojectors=.false.

    norb=1; natom=1
    model%norb=norb; model%natom=natom; model%sample='SRO1orb'
    print*,'The sample chosen is ',model%sample 

	allocate(model%ra(3,natom));  model%ra=0
    allocate(model%orbits(norb)); model%orbits=(/'A1g'/)
	model%periodiclayers=.false.
	
	usegroup=.true.;  group='C4v';  ng=8
    model%group=group; model%ng=ng

    model%U=4; model%mu=1.2+0.1
	nband=1; model%nband=nband
	allocate(model%npatch(nband),model%nativeband(nband),model%pocket(nband),model%mushift(nband))
	model%nativeband=(/1/); model%pocket=(/'G_'/); model%npatch=(/96/); model%mushift=0.


    !---------------------------------tune interactions interactively----------------------
	if(dorpa.or.(useX0.and.appendX0))then
	  print*,'you are now free to tune interactions:'
	  print*,'U='; read*,model%U
	  print*,'Vnn='; read*,model%Vnn
	end if
    !---------------------------------------------------------------------------------------

	Q=0; call checksymmetry(Q);	
	call FermiSurf(norb,Q,model); call bandstructure(norb,Q,model) 
	call meshfermisurface(model)

    nkx=40; nky=40     
    nLform=13; Lcontact=3
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
	call addform( (/2,0,0/)*1.d0,'A1g')
	call addform( (/-2,0,0/)*1.d0,'A1g')
	call addform( (/0,2,0/)*1.d0,'A1g')
	call addform( (/0,-2,0/)*1.d0,'A1g')
	if(iform/=nLform)stop 'nform error @ SrRuO1orb'


	call setLformfunc(1)   !use identity group (by setting ng=1) to setup fundamental lattice forms
	call setOformfunc()
	guess='***';  call setMolecules()

    nw=400; wmax=1.e2
    wmin=1.e-4; wir=1.e-4
    diverge=40

    if(.not.dorpa)return

	call writeinputs()
	call setupworkplace()
	cutoff=1.e-2; call rpa(cutoff,0,'SRO1dX0.dat'); stop

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
end subroutine SrRuO1orb


subroutine SrRuO1orbhk(norb,kv,hk,model,mf,slvbsn)
  use standard_derived_types
  implicit none
  integer norb
  real*8 kv(2)
  complex*16 hk(norb,norb)
  type (modelblock) :: model
  type (mfblock) :: mf
	logical slvbsn

  real*8 t1,t2,mu
  real*8 ck(2),sk(2),vk(2)

  t1=0.8; t2=0.35; mu=model%mu

  vk=kv*3.1415926;  ck=cos(vk)

  hk=-2*t1*sum(ck)-4*t2*ck(1)*ck(2)-mu

  return
end subroutine SrRuO1orbhk
