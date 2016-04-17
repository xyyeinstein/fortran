subroutine LiFeAs2L5d()
  use workplace 
  implicit none

  real*8 Q,x,y,z,tt,factor
  integer ix,iy,iz,a,b,error,nband,nbondxy,nbondz
  logical dorpa
  integer iform
  
  interface
    subroutine CheckHopLiFeAs2L5d(norb,model)
      use standard_derived_types 
      implicit none
      integer norb
      type(modelblock) , target :: model
    end subroutine
  end interface

  DataReady=.false.
  dorpa=.false.
  uniformkmesh=.false.
  Appendkmesh=.false.
  Appendfs=.false.
  useX0=.true.
  appendX0=.false.
  QuickSearch=.false.
  AppendBreakPoint=.false.
  kmeshfile='LFA2L5dkmesh.dat'
  useprojectors=.false.; appendprojectors=.false.; projectorfile='LFA2L5dproj.dat'
  useprojectors=useprojectors.and.(.not.QuickSearch)
  ! this is a test of supercell method, and should be compared to the 5-orbit model

  
    norb=5; natom=1; nambu=norb*2
    model%norb=norb; model%natom=natom; model%sample='LFA2L5d'
    print*,'The sample chosen is ',model%sample 

    allocate(model%ra(3,natom)); model%ra(:,1)=(/0.0,0.0,0.0/)
    model%a=(/1,0,0/); model%b=(/0,1,0/); model%c=(/0,0,1/)

    allocate(model%orbits(norb))
    model%orbits(1:5)=(/'A1g','E1g','E2g','B2g','B1g'/)
    model%periodiclayers=.false.

    usegroup=.true.;  group='S4B';  ng=8     
    model%group=group; model%ng=ng
    model%propergauge=.true.

    nbondxy=6
    nbondz=2
    allocate(model%t(5,5,-nbondxy:nbondxy,-nbondxy:nbondxy,-nbondz:nbondz))
    
    !----------------------------------------------------------
    model%U=1.20; model%JH=0.30; model%Uab=model%U-2*model%JH; 

    open(10,file='NaFeAs2L5d.model')
    read(10,*)model%t
    close(10)

    nband=3; model%nband=nband; model%mu=9.36  !9.48 for LiFeAs
    allocate(model%npatch(nband),model%nativeband(nband))
    allocate(model%pocket(nband),model%mushift(nband))
    model%nativeband=(/2,3,4/)
    model%pocket=(/'M_','GM','XY'/)
    model%npatch=(/32,64,64/)
    model%mushift=0.
    !----------------------------------------------

    Q=0
    call BasisRotateLFA2L5d(norb,model)
    !call FindMaxHopLFA2L5d(norb,model); pause
    call CheckHopLiFeAs2L5d(norb,model)
    call checksymmetry(Q);
    call bandstructure(norb,Q,model)
    call FermiSurf(norb,Q,model)
    call meshfermisurface(model)


    nkx=100
    nky=nkx     
    nLform=9; Lcontact=2
    allocate(Lformtable(nLform)); allocate(Lformbasis(3,nLform))
    iform=0;
    call addform((/0,0,0/)*1.d0,'A1g')
    call addform((/1,0,0/)*1.d0,'A1g')
    call addform((/1,1,0/)*1.d0,'A1g')
    if(nLform==9)then
      call addform( (/-1,0,0/)*1.d0,'A1g')
      call addform( (/0,1,0/)*1.d0,'A1g')
      call addform( (/0,-1,0/)*1.d0,'A1g')
      call addform( (/-1,-1,0/)*1.d0,'A1g')
      call addform( (/1,-1,0/)*1.d0,'A1g')
      call addform( (/-1,1,0/)*1.d0,'A1g')
      call setLformfunc(1)
    else if(nLform==5)then
      call addform((/1,0,0/)*1.d0,'B1g')
      call addform((/1,1,0/)*1.d0,'B2g')
      call setLformfunc(ng)
    end if
    call setOformfunc()

    SkipInterOrbitPair=.false.
    guess='***';  call setMolecules()

    nw=400; wmax=100
    wmin=1.e-4; wir=1.e-5
    wirx0=1.e-5    
    diverge=80

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
end subroutine

subroutine LiFeAs2L5dhk(norb,kv,Q,hk,model,mf,slvbsn)
  use standard_derived_types
  implicit none
  integer norb
  real*8 kv(2),Q
  complex*16 hk(norb,norb)
  type (modelblock), target :: model
  type (mfblock) :: mf
	logical slvbsn

  integer ix,iy,iz,iorb,nbondxy,nbondz
  real*8 r(3),QM(2)
  real*8 pi,p(3)
  complex*16, pointer, dimension (:,:) :: t
  complex*16 one,phase

  one=cmplx(0,1); pi=asin(1.d0)*2
  
  QM=(/pi,pi/)
  p(1:2)=kv*pi; p(3)=Q !-sum(p(1:2))
  !p(1:2)=p(1:2)+QM

  hk=0; 
  nbondxy=6
  nbondz=2
  do iz=-nbondz,nbondz
  do iy=-nbondxy,nbondxy; do ix=-nbondxy,nbondxy
     t=>model%t(:,:,ix,iy,iz)
     r=(/ix,iy,iz/);
	 phase=exp(-one*sum(p*r))
	 hk=hk+t*phase
  end do; end do
  end do

  do iorb=1,norb; hk(iorb,iorb)=hk(iorb,iorb)-model%mu; end do
  
  return
end subroutine

subroutine CheckHopLiFeAs2L5d(norb,model)
  use standard_derived_types 
  implicit none

  integer norb
  type(modelblock) , target :: model
  
  integer norb1atom,atomimage,nbondxy,nbondz
  integer ix,iy,iz,ixg,iyg,izg,ig
  integer iorb,jorb,iatom,jatom,idim,jdim
  integer iorbg,jorbg,iatomg,jatomg,idimg,jdimg
  real*8 r(3),rg(3)
  complex*16, pointer :: t(:,:),tg(:,:)

  integer indexgi(2),indexgj(2)
  real*8 Xgi(2),Xgj(2)

  nbondxy=6
  nbondz=2


  print*,'Check Hopping symmetry ...'
  do iz=-nbondz,nbondz;  do iy=-nbondxy,nbondxy; do ix=-nbondxy,nbondxy
    t=>model%t(:,:,ix,iy,iz)
	r=(/ix,iy,iz/);
	do ig=2,model%ng
	  call groupaction(ig,r,rg,model%group)
	  ixg=nint(rg(1)); iyg=nint(rg(2)); izg=nint(rg(3))
	  tg=>model%t(:,:,ixg,iyg,izg)
      do iorb=1,norb
        call orbitimage(ig,iorb,norb,model%orbits,model%group,indexgi,Xgi)
		iorbg=indexgi(1)
		do jorb=1,norb
          call orbitimage(ig,jorb,norb,model%orbits(1:5),model%group,indexgj,Xgj)
		  jorbg=indexgj(1)
		  if(abs(t(iorb,jorb)-tg(iorbg,jorbg)*Xgi(1)*Xgj(1))>5.e-3)then
			print*,r
			print*,rg
			print*,'  ix ,  iy ,  iz ,  iorb ,jorb'
			write(*,'(5i6)')ix,iy,iz,iorb,jorb
			print*,'ig:',ig
			print*,'ixg  ,iyg  ,izg  ,iorbg,jorbg'
			write(*,'(5i6)')ixg,iyg,izg,iorbg,jorbg
			print*,t(iorb,jorb)
			print*,tg(iorbg,jorbg),Xgi(1),Xgj(1)
			pause
		  end if
		end do
	  end do
	end do
  end do; end do; end do

  print*,'Hoppings satisfy the point group symmetry'
  
  return
end subroutine

subroutine BasisRotateLFA2L5d(norb,model)
  use standard_derived_types
  implicit none

  integer norb
  type(modelblock) :: model
  real*8 U(5,5),V(5,5)
  integer i,ix,iy,iz,nbondxy,nbondz

  nbondxy=6
  nbondz=2
  U=0
  do i=1,5
    U(i,i)=1
  end do
  U(2,2)=1; U(2,3)=1
  U(3,2)=1; U(3,3)=-1
  U(2:3,2:3)=U(2:3,2:3)/sqrt(2.0)
  V=transpose(U)
  do iz=-nbondz,nbondz; do iy=-nbondxy,nbondxy; do ix=-nbondxy,nbondxy
    model%t(:,:,ix,iy,iz)=matmul(U, matmul(model%t(:,:,ix,iy,iz),V))
  end do; end do; end do

  return
end subroutine

subroutine FindMaxHopLFA2L5d(norb,model)
  use standard_derived_types
  implicit none

  integer norb
  type (modelblock) :: model

  integer iorb,jorb,ix,iy,iz,rmax(3),orbits(2)
  real*8 tmax

  tmax=0; rmax=0; orbits=0
  do iorb=1,norb; do jorb=1,norb; 
  do iz=-2,2; do ix=-6,6; do iy=-6,6
    !if(.not.(ix*ix+iy*iy==1))cycle !specifiy on particular bond
	if(ix==0.and.iy==0)cycle !avoid the on-site energy
	!if(.not.(ix==0.and.iy==0.and.iz==0))cycle  !check the on-site energy
	if(.not.(iorb==4.and.jorb==2))cycle
	if(.false.)then
	  if((ix*ix+iy*iy==1).and.iz==1)then
	    if(iorb==4.and.jorb==4)then
	      print*,iorb,jorb,ix,iy,iz
	      print*,model%t(iorb,jorb,ix,iy,iz)
	    end if
	  end if
    else
	  if(abs(model%t(iorb,jorb,ix,iy,iz))>tmax)then
	    tmax=model%t(iorb,jorb,ix,iy,iz)
	    rmax=(/ix,iy,iz/)
	    orbits=(/iorb,jorb/)
	  end if
	end if
  end do; end do; end do
  end do; end do
  
  print*,orbits
  print*,rmax
  print*,tmax
  return
end subroutine
