subroutine square2atom2orb()
  use workplace
  implicit none

  real*8 Q

  DataReady=.false.
  QuickSearch=.false.
  AppendBreakPoint=.false.
  Appendkmesh=.false.; kmeshfile='sq22kmesh.dat'
  useprojectors=.true.; appendprojectors=.false.; projectorfile='sq2a2dproj.dat'
  useprojectors=useprojectors.and.(.not.QuickSearch)

    norb=2; natom=2
	model%norb=norb; model%natom=natom; model%sample='sq2a2d'
    print*,'The sample chosen is ',model%sample; 

	allocate(model%ra(3,natom)); model%ra(:,1)=0; model%ra(:,2)=(/1,1,0/)*0.5
	model%a=(/1,0,0/); model%b=(/0,1,0/); model%c=(/0,0,1/)
    allocate(model%orbits(norb)); model%orbits='A1g'
	model%periodiclayers=.false.
    
	usegroup=.true.;  group='D4h';  ng=8
	model%group=group; model%ng=ng
    
	model%U=3; model%mu=0

    Q=0; call checksymmetry(Q); call fermisurf(norb,Q,model); 

    nkx=40; nky=40  
    nLform=2; Lcontact=2
    allocate(Lformtable(nLform)); Lformtable=(/'A1g','B2g'/)
    allocate(Lformbasis(3,nLform)); Lformbasis(:,1)=0; Lformbasis(:,2)=(/1,1,0/)*0.5
	!Lformbasis(:,3)=(/1,0,0/)

	call setLformfunc(ng)
	call setOformfunc()
    SkipInterOrbitPair=.false.
	guess='B2g';  call setMolecules()

    nw=400; wmax=100
    wmin=1.e-4; wir=1.e-3
    diverge=40

  return
end subroutine square2atom2orb

subroutine square2a2dhk(kv,hk,model,mf,slvbsn)
  use standard_derived_types
  implicit none
  real*8 kv(2),vk(2)
  complex*16 hk(2,2)
  real*8 pi,t1,t2,mu
  type (modelblock) :: model
  type (mfblock) :: mf
	logical slvbsn

  pi=asin(1.d0)*2; vk=kv*pi
  t1=1; t2=-0.3; mu=4*t2+model%mu

  !intracell hopping between atom 1 and 2
  hk(1,2)=-4*t1*cos(vk(1)*0.5)*cos(vk(2)*0.5)
  hk(2,1)=hk(1,2)  

  !intercell hopping between like atoms
  hk(1,1)=-2*t2*sum(cos(vk))-mu; hk(2,2)=hk(1,1)
  
  return
end subroutine square2a2dhk
