subroutine kagome()
  use workplace
  implicit none

  integer iform,idim
  real*8, dimension (3) :: da,db
  integer i,nband
  real*8 Q; data Q/0.d0/

  one=cmplx(0,1)
  pi=asin(1.d0)*2
  twopi=pi*2

  useX0=.true.; appendX0=.true.
  BCSflow=.false.
  QuickSearch=.false.
  DataReady=.false.
  AppendBreakPoint=.false.
  AppendKmesh=.true.; kmeshfile='kagomekmesh.dat'
  useprojectors=.false.; appendprojectors=.false.; projectorfile='kagomeproj.dat'
  useprojectors=useprojectors.and.(.not.QuickSearch)

    norb=3; natom=3    
	model%norb=norb; model%natom=natom; model%sample='kagme'
    print*,'Sample = ',model%sample 

	allocate(model%ra(3,natom),model%orbits(3))
	model%ra(:,1)=0; model%ra(:,2)=(/-0.25,0.25*sqrt(3.),0./)
	model%ra(:,3)=(/0.25,0.25*sqrt(3.),0./); model%orbits='A1g'
	model%a=(/1,0,0/); model%b=(/0.5,0.5*sqrt(3.),0./); model%c=(/0,0,1/)
    model%ra(:,3)=model%ra(:,3)-model%b
	do i=1,3; model%ra(:,i)=model%ra(:,i)-(/0.5,0.,0./); end do
	call dualvector(model%b,model%c,model%a,da); model%ka=da(1:2)*2
	call dualvector(model%c,model%a,model%b,db); model%kb=db(1:2)*2
   
	usegroup=.true.;  group='C2v';  ng=2
    model%group=group; model%ng=ng

    model%U=2; model%Vnn=1.5; model%mu=-2
	model%propergauge=.false.

	nband=1; model%nband=nband
	allocate(model%npatch(nband),model%nativeband(nband),model%pocket(nband),model%mushift(nband))
	model%nativeband=1; model%pocket='HX'; model%npatch=96; model%mushift=0.
	call meshfermisurface(model)

    model%nest=6; model%qnest=0
	model%qnest(:,1)=(/0.,2./sqrt(3.)/)
	do i=2,model%nest
	   model%qnest(:,i)=model%qnest(:,i-1)
	   call rotate60(model%qnest(:,i))
    end do

    !call suscep(norb,model); 
	!call rpasuscep(norb,model); stop

	Q=0; if(.not.model%propergauge)call checksymmetry(Q);	!call FermiSurf(norb,Q,model); !pause 'continue?'
    call BandStructure(norb,Q,model)

	nkx=40; nky=40  
    nw=400; wmax=1.e2
    wmin=1.e-4; wir=1.e-4
    diverge=40

    nLform=19; Lcontact=3
    allocate(Lformtable(nLform))
	allocate(Lformbasis(3,nLform))
	
	!on-site form
	iform=1; Lformbasis(:,iform)=0; Lformtable(iform)='A1g'

    !1st neighbors
    call addform( (/0.5,0.,0./)*1.d0,  'A1g' )
    call addform( (/0.25,0.25*sqrt(3.),0./)*1.d0,  'A1g' )
	call addform( (/-0.25,0.25*sqrt(3.),0./)*1.d0, 'A1g' ) 
    call addform( (/0.5,0.,0./)*1.d0,  'Odd' )
    call addform( (/0.25,0.25*sqrt(3.),0./)*1.d0,  'Odd' )
	call addform( (/-0.25,0.25*sqrt(3.),0./)*1.d0, 'Odd' ) 

    !2nd neighbors
	call addform( (/0.,0.5*sqrt(3.),0./)*1.d0, 'A1g' )
	call addform( (/0.75,0.25*sqrt(3.),0./)*1.d0, 'A1g' )
    call addform( (/-0.75,0.25*sqrt(3.),0./)*1.d0,  'A1g' )
	call addform( (/0.,0.5*sqrt(3.),0./)*1.d0, 'Odd' )
	call addform( (/0.75,0.25*sqrt(3.),0./)*1.d0, 'Odd' )
    call addform( (/-0.75,0.25*sqrt(3.),0./)*1.d0,  'Odd' )

    !3rd neighbors
	call addform( (/1,0,0/)*1.d0, 'A1g' )
    call addform( (/0.5,0.5*sqrt(3.),0./)*1.d0,  'A1g' )
	call addform( (/-0.5,0.5*sqrt(3.),0./)*1.d0, 'A1g' ) 
	call addform( (/1,0,0/)*1.d0, 'Odd' )
    call addform( (/0.5,0.5*sqrt(3.),0./)*1.d0,  'Odd' )
	call addform( (/-0.5,0.5*sqrt(3.),0./)*1.d0, 'Odd' ) 

    if(iform/=nLform)stop 'nLform not correct'

	call setLformfunc(ng)
	call setOformfunc()
	guess='***';  call setMolecules()
    call hexbounds()
	
	group='C6v';  ng=12
    model%group=group; model%ng=ng

	if(useX0.and.appendX0)then
	  print*,'U,Vnn='; read*,model%U,model%Vnn
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
end subroutine kagome

subroutine kagomehk(norb,kv,hk,model,mf,slvbsn)
  use standard_derived_types
  implicit none

  integer norb
  real*8 kv(2)
  complex*16 hk(norb,norb)
  type (modelblock) :: model
  type (mfblock) :: mf
	logical slvbsn

  real*8 pi,bond(2),dr(2)
  complex*16 one
  integer i,j

  pi=asin(1.d0)*2; one=cmplx(0,1)

  hk=0
  do i=1,3; hk(i,i)=-model%mu; end do
  
  !NN-bonds:

  !1-2 bond
  bond=(/-0.25,0.25*sqrt(3.)/)
  dr=0; if(model%propergauge)dr=model%ra(1:2,2)-model%ra(1:2,1)  
  do i=1,-1,-2; bond=bond*i
     hk(1,2)=hk(1,2)-exp(one*pi*sum((bond-dr)*kv))
  end do
 
  !1-3 bond
  bond=(/0.25,0.25*sqrt(3.)/)
  dr=0; if(model%propergauge)dr=model%ra(1:2,3)-model%ra(1:2,1)  
  do i=1,-1,-2; bond=bond*i
     hk(1,3)=hk(1,3)-exp(one*pi*sum((bond-dr)*kv))
  end do

  !2-3 bond
  bond=(/0.5,0./)
  dr=0; if(model%propergauge)dr=model%ra(1:2,3)-model%ra(1:2,2)  
  do i=1,-1,-2; bond=bond*i
     hk(2,3)=hk(2,3)-exp(one*pi*sum((bond-dr)*kv))
  end do

  do j=1,3; do i=1,j-1
     hk(j,i)=conjg( hk(i,j) )
  end do; end do

  return
end subroutine kagomehk
