subroutine BiS2_1band()
  use workplace
  implicit none

  real*8 Q,cutoff
  integer nband,iform,i
  logical dorpa

  useX0=.false.; AppendX0=.false.
  dorpa=.false.
  BCSflow=.false.
  DataReady=.false.
  AppendBreakPoint=.false.
  SkipInterOrbitPair=.false.
  skipOddOrbPair=.false.
  Appendkmesh=.true.; kmeshfile='BiS2_1bkmesh.dat'
  useprojectors=.false.

    norb=1; natom=1
    model%norb=norb; model%natom=natom; model%sample='BiS2_1b'
    print*,'The sample chosen is ',model%sample; 

	allocate(model%ra(3,natom)); model%ra(:,1)=0
    model%a=(/1,0,0/); model%b=(/0,1,0/); model%c=(/0,0,1/)

    allocate(model%orbits(norb)); model%orbits(1:2)='A1g'
	model%periodiclayers=.false.

    usegroup=.true.;  group='C4v';  ng=8    
	model%group=group; model%ng=ng

    model%U=8; model%JH=0; model%Uab=0; model%Vnn=0.5; model%mu=.925+0.2-0.3
	nband=1; model%nband=nband
	allocate(model%npatch(nband),model%nativeband(nband),model%pocket(nband),model%mushift(nband))
	model%nativeband=(/1/); model%pocket=(/'XY'/); model%npatch=(/240/); model%mushift=0.


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
	if(iform/=nLform)stop 'nform error @ BiS2_1band'


	call setLformfunc(1)   !use identity group (by setting ng=1) to setup fundamental lattice forms
	call setOformfunc()
	guess='***';  call setMolecules()

    nw=400; wmax=1.e2
    wmin=1.e-4; wir=4.e-4
    diverge=40

    if(.not.dorpa)return

	call writeinputs()
	call setupworkplace()
	cutoff=1.e-3; call rpa(cutoff,0,'BS1bX0.dat'); stop

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
end subroutine BiS2_1band

subroutine BiS2hk_1band(norb,vk,hk,model,mf,slvbsn)
  use standard_derived_types
  implicit none
  type (modelblock) :: model
  integer norb
  real*8 vk(2)
  complex*16 one
  complex*16 hk(norb,norb)
  type (mfblock) :: mf
	logical slvbsn
  
  complex*16 h(2,2),U(2,2),V(2,2)
  real*8 mu,kx,ky,pi,soc
  real*8 t,t1,t2,txy,txy3,t21xx,t21yy,t21xy

  pi=2*asin(1.) ; one=cmplx(0,1) 
  t=-0.167 ; txy=-0.107 ; t1=0.880 ; t2=0.094 ; txy3=0.028 ; t21xx=0.069 ; t21yy=0.014 ; t21xy=-0.02 
  kx=vk(1)*pi ; ky=vk(2)*pi
  
  mu=-2.811+model%mu 

  h(1,1)=2*t*(cos(kx)+cos(ky))+2*t1*cos(kx+ky)+2*t2*cos(kx-ky)-mu
  h(2,2)=2*t*(cos(kx)+cos(ky))+2*t1*cos(kx-ky)+2*t2*cos(kx+ky)-mu
  h(1,2)=2*txy*(cos(kx)-cos(ky))+2*txy3*(cos(2*kx)-cos(2*ky))

  h(1,1)=h(1,1)+2*t21xx*(cos(2*kx+ky)+cos(kx+2*ky))+2*t21yy*(cos(2*kx-ky)+cos(kx-2*ky))
  h(2,2)=h(2,2)+2*t21yy*(cos(2*kx+ky)+cos(kx+2*ky))+2*t21xx*(cos(2*kx-ky)+cos(kx-2*ky))
  h(1,2)=h(1,2)+2*t21xy*(cos(2*kx+ky)-cos(kx+2*ky)+cos(2*kx-ky)-cos(kx-2*ky))
  h(2,1)=conjg(h(1,2))

  !rotate to px and py basis
  U=reshape((/1,-1,1,1/),(/2,2/))/sqrt(2.)
  V=conjg(transpose(U))
  h=matmul(V,matmul(h(1:2,1:2),U))

  !rotate to chiral basis
  U(:,1)=(/1,1/); U(:,2)=(/one,-one/); U=U/sqrt(2.); V=conjg(transpose(U))
  h = matmul( U,matmul(h,V) )

  !add atomic soc

  soc=0.7
  h(1,1)=h(1,1)-soc; h(2,2)=h(2,2)+soc
  hk = h(1,1) - abs(h(1,2))**2/h(2,2) 
  
  return
end subroutine BiS2hk_1band
