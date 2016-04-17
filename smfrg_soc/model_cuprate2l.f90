subroutine cuprate2L()
  use workplace
  implicit none

  real*8, dimension (3) :: da,db
  integer iform

  pi=asin(1.d0)*2; twopi=pi*2;  one=cmplx(0,1)

  DataReady=.false.
  QuickSearch=.false.
  AppendBreakPoint=.false.
  AppendKmesh=.true.; kmeshfile='cuprate2Lkmesh.dat'
  useprojectors=.false.;  AppendProjectors=.false.; Projectorfile='cuprate2Lproj.dat'
  useprojectors=useprojectors.and.(.not.QuickSearch)  

  norb=4; natom=2    
  model%norb=norb; model%natom=natom; model%sample='cuprate2L'
  print*,'Sample = ',model%sample; pause 'continue?'

  allocate(model%ra(3,natom))
  model%ra(:,1)=0; model%ra(:,2)=(/0,0,1/)
  model%a=(/1,0,0/); model%b=(/0,1,0/); model%c=(/0,0,1/)
  call dualvector(model%b,model%c,model%a,da); model%ka=da(1:2)*2
  call dualvector(model%c,model%a,model%b,db); model%kb=db(1:2)*2
   
  usegroup=.true.;  group='C4v';  ng=8
  !within C2v, x2, y2, x and y are all 1d irreducible representations
  model%group=group; model%ng=ng

  model%U=0.5; model%Vnn=1; model%VL2L=0.1
  model%Vbonding=1.5; model%Vantibonding=-1.5
  model%mu=0.9

  model%propergauge=.false.
  if(.not.model%propergauge)call checkHsymmetry_()
  call checkTsymmetry()
  call FermiSurf(norb,model); call BandStructure(norb,model)

  nkx=40; nky=40   
  nw=400; wmax=1.e2
  wmin=1.e-4; wir=1.e-3
  diverge=40

  nLform=8; Lcontact=2
  allocate(Lformtable(nLform))
  allocate(Lformbasis(3,nLform))

  !on-site form
  iform=0
  call addform( (/0.d0,0.d0,0.d0/), 'A1g' )

  !in-plane 1st neighbor
  call addform( (/1,0,0/)*1.d0, 'A1g')
  !call addform( (/1,0,0/)*1.d0, 'B1g')   
  call addform( (/1,0,0/)*1.d0, 'E1u')
  call addform( (/0,1,0/)*1.d0, 'E2u')

  !interlayer 1st neighbor
  call addform( (/0,0,1/)*1.d0,  'A1g')
  call addform( (/0,0,-1/)*1.d0, 'A1g')


  !in-plane 2nd neighbor
  !call addform( (/1,1,0/), 'A1g')
  !call addform( (/1,1,0/), 'B2g')
  call addform( (/1,1,0/)*1.d0, 'Y+X')
  call addform( (/-1,1,0/)*1.d0, 'Y-X')

  !inter-plane 1st neighbor
  !call addform( (/0,0,1/)*1.d0, 'A1g')
  !call addform( (/0,0,-1/)*1.d0, 'A1g')

  !inter-plane 2nd neighbor
  !call addform( (/1,0,1/), 'A1g')
  !call addform( (/0,1,1/), 'A1g')
  !call addform( (/1,0,1/), 'E1u')
  !call addform( (/0,1,1/), 'E2u')

  !call addform( (/1,0,-1/), 'A1g')
  !call addform( (/0,1,-1/), 'A1g')
  !call addform( (/1,0,-1/), 'E1u')
  !call addform( (/0,1,-1/), 'E2u')

  if(iform/=nLform)stop 'nLform not correct'
	
  call setMolecules_()
	
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
end subroutine cuprate2L


subroutine cuprate2Lhk(norb,kv,hk,model)
  use standard_derived_types
  implicit none

  integer norb
  real*8 kv(2)
  complex*16 hk(norb,norb)
  type (modelblock) :: model

  real*8 pi,t1,mu,la,tperb
  complex*16 one
  integer i

  pi=asin(1.d0)*2; one=cmplx(0,1)
    
  t1=-0.3; tperb=2; la=1.; mu=4*t1+model%mu

  hk=0

  !1=aup, 2=bup, 3=adn, 4=bdn

  !intralayer normal part
  hk(1,1)=-2*sum(cos(kv*pi))-4*t1*cos(kv(1)*pi)*cos(kv(2)*pi)-mu
  do i=2,4; hk(i,i)=hk(1,1); end do

  !intralayer Rashba
  hk(1,3)=la*( -sin(kv(2)*pi) - one*sin(kv(1)*pi) )
  hk(3,1)=conjg(hk(1,3))
  hk(2,4)=hk(1,3); hk(4,2)=hk(3,1)

  !interlayer hopping
  hk(1,2)=-tperb; hk(2,1)=-tperb
  hk(3,4)=-tperb; hk(4,3)=-tperb

  return
end subroutine cuprate2Lhk

