subroutine cuprate()
  use workplace
  implicit none

  real*8, dimension (3) :: da,db
  integer iform

  bcsflow=.false.
  dataready=.false.
  quicksearch=.false.
  appendbreakpoint=.false.
  appendkmesh=.false.
  appendprojectors=.false.
  useprojectors=.false.
      pi=asin(1.d0)*2; twopi=pi*2;  one=cmplx(0,1)

  kmeshfile='cupratemesh.dat'
  Projectorfile='cuprateproj.dat'
  useprojectors=useprojectors.and.(.not.QuickSearch)  

  norb=2; natom=1    
  model%norb=norb; model%natom=natom; model%sample='cuprt'
  print*,'Sample = ',model%sample; !pause 'continue?'

  allocate(model%ra(3,natom))
  model%ra(:,1)=0
  model%a=(/1,0,0/); model%b=(/0,1,0/); model%c=(/0,0,1/)
  call dualvector(model%b,model%c,model%a,da); model%ka=da(1:2)*2
  call dualvector(model%c,model%a,model%b,db); model%kb=db(1:2)*2
   
  usegroup=.true.;  group='C4v';  ng=8
  !within C2v, x2, y2, x and y are all 1d irreducible representations
  model%group=group; model%ng=ng
  

  !-------------------------------------------------------------------!
  !modelA1:
 ! model%U=8.0; model%Vnn=-2.0; model%mu=-3.0;model%t1=-0.0; model%la=0.01
  !modelA2:
  !model%U=8.5; model%Vnn=-2.4; model%mu=-2.0; model%t1=0.0; model%la=0.01
  !modelB1:
  !model%U=2.1; model%Vnn=0.0; model%mu=-0.05; model%t1=-0.475; model%la=0.01
  !modelB2:
  !model%U=2.2; model%Vnn=0.0; model%mu=-0.05; model%t1=-0.475; model%la=0.01
  !modelB3:
  !model%U=2.5; model%Vnn=0.0; model%mu=-0.10; model%t1=-0.475; model%la=0.02
  !modelB4
  model%U=2.6; model%Vnn=0.0; model%mu=-0.10; model%t1=-0.475; model%la=0.02
  !modelB5:
  !model%U=4.0; model%Vnn=0.0; model%mu=-0.10; model%t1=-0.45; model%la=0.01
  !-------------------------------------------------------------------!
  model%vl2l=0
  model%vbonding=0
  model%jnn=0

  model%propergauge=.false.
  if(.not.model%propergauge)call checkHsymmetry_()
  call checkTsymmetry()
  call FermiSurf(norb,model); call BandStructure(norb,model)

  nkx=40; nky=40   
  !nw=400; wmax=1.e2
  !wmin=1.e-4; wir=1.e-3
  !diverge=40

  nLform=9   !open(10,file='nlform.input'); read(10,*)nLform; close(10)
  Lcontact=2
  allocate(Lformtable(nLform))
  allocate(Lformbasis(3,nLform))

  !on-site form
  iform=0
  call addform( (/0.d0,0.d0,0.d0/), 'A1g' )

  if(group=='C2v')then
    !in-plane 1st neighbor for C2v
    call addform( (/1,0,0/)*1.d0, 'Ax2')
    call addform( (/0,1,0/)*1.d0, 'Ay2')   
    call addform( (/1,0,0/)*1.d0, 'E1u')
    call addform( (/0,1,0/)*1.d0, 'E2u')
  else if(group=='C4v')then
    !in-plane 1st neighbor for C4v
    call addform( (/1,0,0/)*1.d0, 'A1g')
    call addform( (/1,0,0/)*1.d0, 'B1g')   
    call addform( (/1,0,0/)*1.d0, 'E1u')
    call addform( (/0,1,0/)*1.d0, 'E2u')
  end if

  !in-plane 2nd neighbor
  call addform( (/1,1,0/)*1.d0, 'A1g')
  call addform( (/1,1,0/)*1.d0, 'B2g')
  call addform( (/1,1,0/)*1.d0, 'Y+X')
  call addform( (/-1,1,0/)*1.d0, 'Y-X')

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
	
  nw=400
  wmax=100
  wmin=1.d-4
  wir=1.d-5
  diverge=40

  call setMolecules_()

  call cuprate_edgestate(edgemodel)
	
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
end subroutine cuprate

subroutine cupratehk(norb,kv,hk,model)
  use standard_derived_types
  implicit none

  integer norb
  real*8 kv(2)
  complex*16 hk(norb,norb)
  type (modelblock) :: model

  real*8 pi,t1,mu,la
  complex*16 one

  pi=asin(1.d0)*2; one=cmplx(0,1)
    
  t1=model%t1; la=model%la;
  mu=model%mu

  hk=0
  hk(1,1)=-2*sum(cos(kv*pi))-4*t1*cos(kv(1)*pi)*cos(kv(2)*pi)-mu
  hk(2,2)=hk(1,1)

  hk(1,2)=la*( -sin(kv(2)*pi) - one*sin(kv(1)*pi) )
  hk(2,1)=conjg(hk(1,2))

  return
end subroutine cupratehk

subroutine cuprate_edgestate(edgemodel)
  use standard_derived_types
  implicit none
  type (edgemodelblock) :: edgemodel

  integer nb,nq,Lb,norb,nambu,ndim,natom
  real*8 da(3),db(3)

  natom=1; norb=2; nambu=norb*2
  nb=1; Nq=30; Lb=60; ndim=nambu*Lb
  open(37,file='edgeinput.dat')
  write(37,*)Nq
  write(37,*)Lb
  write(37,*)nambu
  write(37,*)ndim
  close(37)

  edgemodel%sample='cuprt'; edgemodel%periodicslab=.false.
  edgemodel%natom=natom; edgemodel%norb=norb; edgemodel%nambu=nambu 
  edgemodel%nb=nb;  edgemodel%nq=nq; edgemodel%Lb=Lb;  edgemodel%ndim=ndim
  print*,'Sample = ',edgemodel%sample; !pause 'continue?'

  allocate(edgemodel%ra(3,natom))
  edgemodel%ra(:,1)=0
  edgemodel%a=(/1,0,0/); edgemodel%b=(/0,1,0/); edgemodel%c=(/0,0,1/)
  call dualvector(edgemodel%b,edgemodel%c,edgemodel%a,da); edgemodel%ka=da(1:2)*2
  call dualvector(edgemodel%c,edgemodel%a,edgemodel%b,db); edgemodel%kb=db(1:2)*2

  return
end subroutine cuprate_edgestate

subroutine cuprategk(norb,vk,gk)
  use standard_derived_types
  use workplace, only : Mformfunc
  implicit none

  integer norb
  real*8 vk(2)
  complex*16 gk(norb,norb)

  integer iS,S,spin,nbspin,idim,stat
  real*8 a,u1,u2
  complex*16 fk,u,formk
  type (Mformconfig), pointer :: Mi
  
  integer error

  a=.3; gk=0
  
  open(50,file='ppmode.dat')
  read(50,*)iS,S
  idim=0
  do while(.true.)
     read(50,*,iostat=error)u1,u2
	 if(error/=0)exit
	 u=cmplx(u1,u2)
	 !read(50,*,iostat=stat)u; if(stat<0)exit
	 idim=idim+1
	 Mi=>Mformfunc(idim)
	 spin=Mi%orb; nbspin=Mi%nborb
	 fk=formk(vk,Mi%Lformfunc)
	 gk(spin,nbspin)=gk(spin,nbspin)+u*fk
  end do
  close(50)

  gk=gk*a
  
  !gk=conjg(transpose(gk))

  return
end subroutine

 
