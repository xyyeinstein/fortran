subroutine lee()
  use workplace
  implicit none

  real*8, dimension (3) :: da,db
  integer iform

  pi=asin(1.d0)*2; twopi=pi*2;  one=cmplx(0,1)

  DataReady=.false.
  QuickSearch=.false.
  AppendBreakPoint=.false.
  AppendKmesh=.true.; kmeshfile='leekmesh.dat'
  useprojectors=.true.;  AppendProjectors=.true.; Projectorfile='leeproj.dat'
  useprojectors=useprojectors.and.(.not.QuickSearch)  

  norb=4; natom=2    
  model%norb=norb; model%natom=natom; model%sample='lee'
  print*,'Sample = ',model%sample; pause 'continue?'

  allocate(model%ra(3,natom))
  model%ra(:,1)=0; model%ra(:,2)=(/0,0,1/)
  model%a=(/1,0,0/); model%b=(/0,1,0/); model%c=(/0,0,1/)
  call dualvector(model%b,model%c,model%a,da); model%ka=da(1:2)*2
  call dualvector(model%c,model%a,model%b,db); model%kb=db(1:2)*2
   
  usegroup=.true.;  group='C4v';  ng=8
  model%group=group; model%ng=ng

  model%U=2; model%Vnn=0; model%VL2L=0; model%mu=0

  model%propergauge=.false.
  if(.not.model%propergauge)call checkHsymmetry_()
  call checkTsymmetry()
  call FermiSurf(norb,model); call BandStructure(norb,model)

  nkx=40; nky=40   
  nw=400; wmax=1.e2
  wmin=1.e-4; wir=1.e-3
  diverge=40

  nLform=6; Lcontact=2
  allocate(Lformtable(nLform))
  allocate(Lformbasis(3,nLform))

  !on-site form
  iform=0
  call addform( (/0.d0,0.d0,0.d0/), 'A1g' )

  !in-plane 1st neighbor
  call addform( (/1,0,0/)*1.d0, 'A1g')
  !call addform( (/0,1,0/)*1.d0, 'B1g')
  call addform( (/1,0,0/)*1.d0, 'E1u')
  call addform( (/0,1,0/)*1.d0, 'E2u')

  !in-plane 2nd neighbor
  !call addform( (/1,1,0/), 'A1g')
  !call addform( (/1,1,0/), 'B2g')
  !call addform( (/1,1,0/), 'E1u')
  !call addform( (/1,1,0/), 'E2u')

  !inter-plane 1st neighbor
  call addform( (/0,0,1/)*1.d0, 'A1g')
  call addform( (/0,0,-1/)*1.d0, 'A1g')

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
end subroutine lee


subroutine leehk(norb,kv,hk,model)
  use standard_derived_types
  implicit none

  integer norb
  real*8 kv(2)
  complex*16 hk(norb,norb)
  type (modelblock) :: model

  real*8 pi,mu,del,la,t,A,ek,skx,sky
  complex*16 one
  integer ia,ja,ispin,jspin,i,isign,idim,jdim

  pi=asin(1.d0)*2; one=cmplx(0,1)
  
  !mu=2.7; del=1; la=1; t=0.9; a=0.1  
  mu=0; del=1.6; la=1; t=1; a=0.  

  hk=0

  ek=-t*sum(cos(kv*pi)) + a * ( 1-cos(kv(1)*pi)*cos(kv(2)*pi) ) -mu -model%mu
 
  do i=1,norb; hk(i,i)=ek; end do

  skx=sin(kv(1)*pi); sky=sin(kv(2)*pi)

  do ja=1,2; do jspin=1,2; jdim=ja+(jspin-1)*2
  do ia=1,2; do ispin=1,2; idim=ia+(ispin-1)*2

     if(ia==ja.and.ispin/=jspin)then
	   isign=1; if(ispin<jspin)isign=-1
	   !hk(idim,jdim)=hk(idim,jdim)+la*(skx+isign*sky*one)
	   hk(idim,jdim)=hk(idim,jdim)+la*(skx*isign*one-sky)
	 end if

	 if(ia/=ja.and.ispin==jspin)hk(idim,jdim)=hk(idim,jdim)+del
  end do; end do; end do; end do
  
  return
end subroutine leehk

