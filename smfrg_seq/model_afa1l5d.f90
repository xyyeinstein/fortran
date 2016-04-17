subroutine LiFeAs1L5d()
  use workplace 
  implicit none

  real*8 Q,x,y,z,tt
  integer ix,iy,iz,a,b,error,nband,nbondxy,nbondz
  logical dorpa
  integer icase
  
  interface
    subroutine CheckHopLiFeAs1L5d(norb,model)
      use standard_derived_types 
      implicit none
      integer norb
      type(modelblock) , target :: model
    end subroutine
  end interface

  dorpa=.false.
  Appendkmesh=.false.
  useX0=.false.
  appendX0=.false.
  DataReady=.false.
  QuickSearch=.false.
  AppendBreakPoint=.false.
  kmeshfile='lfa1l5dkmesh.dat'
  useprojectors=.false.; appendprojectors=.false.; projectorfile='lfs1l5dproj.dat'
  useprojectors=useprojectors.and.(.not.QuickSearch)
  ! this is a test of supercell method, and should be compared to the 5-orbit model

  
    norb=5; natom=1
    model%norb=norb; model%natom=natom; model%sample='LFA1L5d'
    print*,'The sample chosen is ',model%sample 

	allocate(model%ra(3,natom)); model%ra(:,1)=(/0.0,0.0,0.0/)
    model%a=(/1,0,0/); model%b=(/0,1,0/); model%c=(/0,0,1/)

    allocate(model%orbits(norb)); model%orbits(1:5)=(/'A1g','E2g','E1g','B2g','B1g'/)
	model%periodiclayers=.false.

    usegroup=.true.;  group='D4h';  ng=8     
	model%group=group; model%ng=ng
	model%propergauge=.true.

    model%U=4/2; model%JH=0.54/2; model%Uab=model%U-2*model%JH; 

	nbondxy=5
	nbondz=0
	allocate(model%t(5,5,-nbondxy:nbondxy,-nbondxy:nbondxy,-nbondz:nbondz))
    open(10,file='LiFeAs1L5d.model')
	  read(10,*)model%t
	close(10)

    call BasisRotateLFA1L5d(norb,model)
	call CheckHopLiFeAs1L5d(norb,model)
    
    nband=2; model%nband=nband; model%mu=9.48    !9.48 for LiFeAs
    allocate(model%npatch(nband),model%nativeband(nband),model%pocket(nband),model%mushift(nband))
    model%nativeband=(/3,4/);	model%pocket=(/'GM','XY'/); model%npatch=240; model%mushift=0.
	Q=0
	call checksymmetry(Q);	
	call bandstructure(norb,Q,model)
    call FermiSurf(norb,Q,model)
    call meshfermisurface(model)


	nkx=40; nky=40     
    nLform=3; Lcontact=2
    allocate(Lformtable(nLform)); Lformtable=(/'A1g','A1g','B2g'/)
    allocate(Lformbasis(3,nLform));	Lformbasis(:,1)=0
	Lformbasis(:,2)=(/1,1,0/); Lformbasis(:,3)=(/1,1,0/)

	call setLformfunc(ng)
	call setOformfunc()

    SkipInterOrbitPair=.false.
	guess='A1g';  call setMolecules()

    nw=400; wmax=100
    wmin=1.e-4; wir=1.e-3
    diverge=40

  return
end subroutine

subroutine LiFeAs1L5dhk(norb,kv,Q,hk,model,mf,slvbsn)
  use standard_derived_types
  implicit none
  integer norb
  real*8 kv(2),Q
  complex*16 hk(norb,norb)
  type (modelblock), target :: model
	type (mfblock) :: mf
	logical slvbsn

  integer ix,iy,iz,iorb,nbondxy,nbondz
  real*8 r(3)
  real*8 pi,p(3)
  complex*16, pointer, dimension (:,:) :: t
  complex*16 one,phase

  one=cmplx(0,1); pi=asin(1.d0)*2
  
  p(1:2)=kv*pi; p(3)=Q !-sum(p(1:2))
 
  hk=0; 
  nbondxy=5
  do iy=-nbondxy,nbondxy; do ix=-nbondxy,nbondxy
     t=>model%t(:,:,ix,iy,0)
     r=(/ix,iy,0/);
	 phase=exp(-one*sum(p*r))
	 hk=hk+t*phase
  end do; end do

  do iorb=1,norb; hk(iorb,iorb)=hk(iorb,iorb)-model%mu; end do
  
  return
end subroutine

subroutine CheckHopLiFeAs1L5d(norb,model)
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

  nbondxy=5
  nbondz=0


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
		  if(abs(t(iorb,jorb)-tg(iorbg,jorbg)*Xgi(1)*Xgj(1))>1.e-4)then
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

subroutine BasisRotateLFA1L5d(norb,model)
  use standard_derived_types
  implicit none

  integer norb
  type(modelblock) :: model
  real*8 U(5,5),V(5,5)
  integer i,ix,iy,iz,nbondxy,nbondz

  nbondxy=5
  nbondz=0
  U=0
  do i=1,5
    U(i,i)=1
  end do
  U(2,2)=-1; U(2,3)=1
  U(3,2)=1; U(3,3)=1
  U(2:3,2:3)=U(2:3,2:3)/sqrt(2.0)
  V=transpose(U)
  do iz=-nbondz,nbondz; do iy=-nbondxy,nbondxy; do ix=-nbondxy,nbondxy
    model%t(:,:,ix,iy,iz)=matmul(U, matmul(model%t(:,:,ix,iy,iz),V))
  end do; end do; end do

  return
end subroutine
