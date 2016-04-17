subroutine BiS2()
  use workplace
  implicit none

  real*8 Q,cutoff
  integer nband,iform

  !useX0=.true.; AppendX0=.false.
  BCSflow=.false.
  DataReady=.false.
  AppendBreakPoint=.false.
  SkipInterOrbitPair=.false.
  skipOddOrbPair=.false.
  Appendkmesh=.false.; kmeshfile='BiS2kmesh.dat'
  useprojectors=.false.

    norb=2; natom=1
    model%norb=norb; model%natom=natom; model%sample='BiS2'
    print*,'The sample chosen is ',model%sample

	allocate(model%ra(3,natom)); model%ra(:,1)=0
    model%a=(/1,0,0/); model%b=(/0,1,0/); model%c=(/0,0,1/)

    allocate(model%orbits(norb)); model%orbits(1:2)=(/'E1u','E2u'/)
	model%periodiclayers=.false.

    usegroup=.true.;  group='C4v';  ng=8    
	model%group=group; model%ng=ng

    !mu=1.15
    !model%U=1; model%JH=0.5; model%Uab=model%U-2*model%JH; model%Vnn=0.; model%mu=1.15  
	!nband=2; model%nband=nband
	!allocate(model%npatch(nband),model%nativeband(nband),model%pocket(nband),model%mushift(nband))
	!model%nativeband=(/1,2/); model%pocket=(/'GM','XY'/); model%npatch=(/96,48/); model%mushift=0.

    !mu=1.25
    model%U=4.5; model%JH=0.; model%Uab=model%U-2*model%JH; model%Vnn=0.6; model%mu=0.8
	nband=2; model%nband=nband
	allocate(model%npatch(nband),model%nativeband(nband),model%pocket(nband),model%mushift(nband))
	model%nativeband=(/1,2/); model%pocket=(/'XY','XY'/); model%npatch=(/96,96/); model%mushift=0.

    !mu=1
    !model%U=4.; model%JH=0; model%Uab=model%U-2*model%JH; model%Vnn=0.5; model%mu=1.05
	!nband=1; model%nband=nband
	!allocate(model%npatch(nband),model%nativeband(nband),model%pocket(nband),model%mushift(nband))
	!model%nativeband=1;	model%pocket='XY'; model%npatch=120; model%mushift=0.


    !---------------------------------tune interactions interactively----------------------
	if(useX0.and.appendX0)then
	  print*,'Since X0 is appended, you are now free to tune interactions:'
	  print*,'U='; read*,model%U
	  print*,'JH='; read*,model%JH; model%Uab=model%U-2*model%JH
	  print*,'Vnn='; read*,model%Vnn
	end if
    !---------------------------------------------------------------------------------------

	Q=0; call checksymmetry(Q);	
	call FermiSurf(norb,Q,model); call bandstructure(norb,Q,model)
	call meshfermisurface(model)

    nkx=40; nky=40     
    nLform=5; Lcontact=2
	allocate(Lformbasis(3,nLform),Lformtable(nLform))

    iform=0
	call addform( (/0,0,0/)*1.d0,'A1g')
	!call addform( (/1,0,0/)*1.d0,'A1g')
	!call addform( (/1,0,0/)*1.d0,'E1u')
	!call addform( (/0,1,0/)*1.d0,'A1g')
	!call addform( (/0,1,0/)*1.d0,'E2u')
	call addform( (/1,1,0/)*1.d0,'A1g')
	call addform( (/1,-1,0/)*1.d0,'A1g')
	call addform( (/-1,1,0/)*1.d0,'A1g')
	call addform( (/-1,-1,0/)*1.d0,'A1g')
	!call addform( (/2,0,0/)*1.d0,'A1g')
	!call addform( (/0,2,0/)*1.d0,'A1g')
	!call addform( (/2,0,0/)*1.d0,'E1u')
	!call addform( (/0,2,0/)*1.d0,'E2u')
	if(iform/=nLform)stop 'nform error @ BiS2'


	call setLformfunc(1)   !use identity group (by setting ng=1) to setup fundamental lattice forms
	call setOformfunc()
	guess='***';  call setMolecules()

    nw=400; wmax=1.e2
    wmin=1.e-4; wir=1.e-3
    diverge=40


	!call writeinputs()
	!call setupworkplace()
	!cutoff=1.e-2; call rpa(cutoff,1,'BiS2X0.dat'); stop

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
end subroutine BiS2

subroutine BiS2hk(ndim,vk,Q,hk,model,mf,slvbsn)
  use standard_derived_types
  implicit none
  type (modelblock) :: model
  integer ndim
  real*8 vk(2),pi,Q
  complex*16 one
  complex*16 hk(ndim,ndim),U(2,2),V(2,2)
  real*8 mu
  real*8 kx,ky
  real*8 t,t1,t2,txy,txy3,t21xx,t21yy,t21xy,soc
  type (mfblock) :: mf
	logical slvbsn

  pi=2*asin(1.) ; one=cmplx(0,1) 
  t=-0.167 ; txy=-0.107 ; t1=0.880 ; t2=0.094 ; txy3=0.028 ; t21xx=0.069 ; t21yy=0.014 ; t21xy=-0.02 
  kx=vk(1)*pi ; ky=vk(2)*pi
  
  mu=-2.811+model%mu 

  hk(1,1)=2*t*(cos(kx)+cos(ky))+2*t1*cos(kx+ky)+2*t2*cos(kx-ky)-mu
  hk(2,2)=2*t*(cos(kx)+cos(ky))+2*t1*cos(kx-ky)+2*t2*cos(kx+ky)-mu
  hk(1,2)=2*txy*(cos(kx)-cos(ky))+2*txy3*(cos(2*kx)-cos(2*ky))
  hk(2,1)=conjg(hk(1,2))

  hk(1,1)=hk(1,1)+2*t21xx*(cos(2*kx+ky)+cos(kx+2*ky))+2*t21yy*(cos(2*kx-ky)+cos(kx-2*ky))
  hk(2,2)=hk(2,2)+2*t21yy*(cos(2*kx+ky)+cos(kx+2*ky))+2*t21xx*(cos(2*kx-ky)+cos(kx-2*ky))
  hk(1,2)=hk(1,2)+2*t21xy*(cos(2*kx+ky)-cos(kx+2*ky)+cos(2*kx-ky)-cos(kx-2*ky))
  hk(2,1)=conjg(hk(1,2))

  !rotate to px and py basis
  U=reshape((/1,-1,1,1/),(/2,2/))/sqrt(2.)
  V=conjg(transpose(U))
  hk=matmul(V,matmul(hk(1:2,1:2),U))
  
  return
end subroutine BiS2hk
