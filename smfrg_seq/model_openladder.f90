subroutine openladder()
  use workplace
  implicit none

  !notice: 1) this subprogram applies for open ladders, a realistic 1d system. in this case model%nleg must be set to 1 since the ladder is open in the rung direction,
  !           although there may be  many legs on the ladder (they are reflected in the hamiltonian oppenladderhk).
  !        2) the y-coordinates of the atoms are set in such a way that no atoms can be found at (x,y)+/-b. on the other hand, the positions are symmetric wrt x-axis.
  !        3) the basis bonds for form factors are given in such a way that under E and mx of C2v the bonds always connect atoms if the starting atom is appropriate.
  !        4) during the calculation the full C2v group is applied if 2) is implemented correctly.  

  real*8 Q
  integer nband

  QuickSearch=.false.
  outputXq=.false.
  DataReady=.false.
  AppendBreakPoint=.false.
  SkipInterOrbitPair=.false.
  Appendkmesh=.false.; kmeshfile='ladderkmesh.dat'
  useprojectors=.false.; Appendprojectors=.true.; projectorfile='sq1a1dPrj.dat'
  useprojectors=useprojectors.and.(.not.QuickSearch)   !projectors are not used in the QuickSearch mode.

    norb=2; natom=2
	model%norb=norb; model%natom=natom; model%sample='ladder'
    print*,'The sample chosen is ',model%sample 
	
	allocate(model%ra(3,natom)); model%ra(:,1)=(/0.,-0.2,0./); model%ra(:,2)=(/0.,0.2,0./)
    allocate(model%orbits(norb)); model%orbits='A1g'
	model%a=(/1,0,0/); model%b=(/0,1,0/); model%c=(/0,0,1/)
	model%ladder=.true.; model%nleg=1   
    
	usegroup=.true.;  group='C2v';  ng=2    !set ng=2 so that only E and mx operations are used to setup the lattice form factors
    model%group=group; model%ng=ng
    
	model%U=4; model%Vnn=0.; model%Jnn=0; model%mu=0.  
	model%nest=1; model%qnest(:,1)=(/1,0/)
    
	Q=0;  call checksymmetry(Q); call FermiSurf(norb,Q,model); call bandstructure(norb,Q,model)  

    nLform=9; Lcontact=2
    allocate(Lformtable(nLform)); allocate(Lformbasis(3,nLform))
	Lformbasis(:,1)=0; Lformtable(1)='A1g'                   !onsite 
	Lformbasis(:,2)=(/1,0,0/); Lformtable(2)='A1g'           !x-bonds connecting 1-1 or 2-2 atoms
	Lformbasis(:,3)=(/1,0,0/); Lformtable(3)='E1u'           !x-bonds connecting 1-1 or 2-2 atoms
	Lformbasis(:,4)=(/0.,0.4,0./); Lformtable(4)='A1g'       !y-bond connecting 1-2 atoms
	Lformbasis(:,5)=(/0.,-0.4,0./); Lformtable(5)='A1g'      !y-bond connecting 2-1 atoms
	Lformbasis(:,6)=(/1.,0.4,0./); Lformtable(6)='A1g'       !diag-bonds connecting 1-2 atoms
	Lformbasis(:,7)=(/1.,-0.4,0./); Lformtable(7)='A1g'      !diag-bonds connecting 2-1 atoms
	Lformbasis(:,8)=(/1.,0.4,0./); Lformtable(8)='E1u'       !diag-bonds connecting 1-2 atoms
	Lformbasis(:,9)=(/1.,-0.4,0./); Lformtable(9)='E1u'      !diag-bonds connecting 2-1 atoms

	call setLformfunc(ng)

	call setOformfunc()

	guess='***';  call setMolecules()

    group='C2v';  ng=4                   !reset ng=4 so that the full C2v operations are applied during the calculations
    model%group=group; model%ng=ng
    
	nw=400; wmax=100
    wmin=1.e-4; wir=1.e-4
    diverge=40
  return
end subroutine openladder

subroutine openladderhk(norb,kv,Q,hk,model,mf,slvbsn)
  use standard_derived_types
  implicit none
  integer norb
  real*8 kv(2),Q
  complex*16 hk(norb,norb),one
  type (modelblock) :: model
  type (mfblock) :: mf
	logical slvbsn

  real*8 pi,ck(2),t1,t2,mu,dr(2)
  integer i,ii

  if(norb/=2.or.model%natom/=2)stop 'open ladder with more than 2 atom per unit cell is not yet defined @ openladderhk'
  pi=asin(1.d0)*2; one=cmplx(0,1); t1=1; t2=-0.; mu=4*t2+model%mu

  hk=0; do i=1,norb; hk(i,i)=-mu-2*t1*cos(kv(1)*pi); end do

  dr=model%ra(1:2,2)-model%ra(1:2,1)
  hk(1,2)=hk(1,2)-t1*exp(one*kv(2)*pi*dr(2))

  dr=model%ra(1:2,2)-model%ra(1:2,1)+(/1,0/)	 
  hk(1,2)=hk(1,2)-t2*exp(one*sum(kv*dr)*pi)

  dr=model%ra(1:2,2)-model%ra(1:2,1)+(/-1,0/)	 
  hk(1,2)=hk(1,2)-t2*exp(one*sum(kv*dr)*pi)

  hk(2,1)=conjg(hk(1,2))

  return
end subroutine openladderhk
