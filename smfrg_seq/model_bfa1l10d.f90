subroutine BFA1L10d()
  use workplace 
  implicit none

  real*8 Q,x,y,z,tt
  integer ix,iy,iz,a,b,error,nband,nbondxy,nbondz
  logical dorpa
  integer iform
  
  interface
    subroutine CheckHopBaFeAs1L10d(norb,model)
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
  kmeshfile='lfa1l10dkmesh.dat'
  useprojectors=.false.; appendprojectors=.false.; projectorfile='lfs1l10dproj.dat'
  useprojectors=useprojectors.and.(.not.QuickSearch)
  ! this is a test of supercell method, and should be compared to the 5-orbit model

  
    norb=10; natom=2
    model%norb=norb; model%natom=natom; model%sample='BFA1L10d'
    print*,'The sample chosen is ',model%sample 

	allocate(model%ra(3,natom)); model%ra(:,1)=(/0.0,0.0,0.0/); model%ra(:,2)=(/0.5,0.5,0.0/)
    model%a=(/1,0,0/); model%b=(/0,1,0/); model%c=(/0,0,1/)

    allocate(model%orbits(norb)); model%orbits(1:5)=(/'A1g','E1g','E2g','B1g','B2g'/)
	model%orbits(6:10)=model%orbits(1:5)
	model%periodiclayers=.false.

    usegroup=.true.;  group='S4A';  ng=8   
	model%group=group; model%ng=ng
	model%propergauge=.false.

    model%U=4/2; model%JH=0.54/2; model%Uab=model%U-2*model%JH; 

	nbondxy=8
	nbondz=4
	open(10,file='BaFeAs1L10d.model')
	allocate(model%t(10,10,-nbondxy:nbondxy,-nbondxy:nbondxy,-nbondz:nbondz))
    
    model%zerotz=.false.
	model%periodiclayers=.true.
	model%t=0
	do while (.true.)
      read(10,*,iostat=error)a,b,x,y,z,tt
      ix=nint(2*x); iy=nint(2*y); iz=nint(2*z)
      model%t(a,b,ix,iy,iz)=tt
      if(error/=0)exit
    end do
	close(10)
    call BasisRotateBaFeAs1L10d(norb,model)
	call CheckHopBaFeAs1L10d(norb,model)
    
    
	nband=4; model%nband=nband; model%mu=11.6
    allocate(model%npatch(nband),model%nativeband(nband),model%pocket(nband),model%mushift(nband))
    model%nativeband=(/5,6,7,8/);	model%pocket=(/'G_','G_','M_','M_'/); model%npatch=60; model%mushift=0.
	Q=0
	!call checkHsymmetry(Q); 	
	call bandstructure(norb,Q,model); stop
    call FermiSurf(norb,Q,model)
    call meshfermisurface(model)


	nkx=40; nky=40     
    nLform=3; Lcontact=2
    allocate(Lformtable(nLform)); allocate(Lformbasis(3,nLform))
	iform=0
	call addform( (/0,0,0/)*1.d0,'A1g')
    call addform( (/0.5,0.5,0./)*1.d0,'A1g')
    call addform( (/0.5,0.5,0./)*1.d0,'B2g')
    !call addform( (/-0.5,0.5,0./)*1.d0,'A1g')
    !call addform( (/-0.5,-0.5,0./)*1.d0,'A1g')
    !call addform( (/0.5,-0.5,0./)*1.d0,'A1g')
    call addform( (/1,0,0/)*1.d0,'A1g')
    call addform( (/1,0,0/)*1.d0,'B1g')
    !call addform( (/-1,0,0/)*1.d0,'A1g')
    !call addform( (/0,1,0/)*1.d0,'A1g')
    !call addform( (/0,-1,0/)*1.d0,'A1g')

	call setLformfunc(ng)
	call setOformfunc()

    SkipInterOrbitPair=.false.
	guess='A1g';  call setMolecules()

    nw=400; wmax=100
    wmin=1.e-4; wir=1.e-5
    diverge=40

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

subroutine BFA1L10dhk(norb,kv,Q,hk,model,mf,slvbsn)
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
  
  p(1:2)=kv*pi; p(3)=Q
 
  hk=0; 
  nbondxy=8
  nbondz=3
  do iz=-nbondz,nbondz; if(model%zerotz.and.iz/=0)cycle; do iy=-nbondxy,nbondxy; do ix=-nbondxy,nbondxy
     t=>model%t(:,:,ix,iy,iz)
     r=(/ix,iy,iz/)
	 r=r*0.5
	 phase=exp(-one*sum(p*r))
	 hk=hk+t*phase
  end do; end do; end do

  if(model%propergauge)then
    hk(1:5,6:10)=hk(1:5,6:10)*exp(-one*sum(p*(model%ra(:,1)-model%ra(:,2))))
    hk(6:10,1:5)=hk(6:10,1:5)*exp(one*sum(p*(model%ra(:,1)-model%ra(:,2))))
  end if

  do iorb=1,norb; hk(iorb,iorb)=hk(iorb,iorb)-model%mu; end do
  
  return
end subroutine

subroutine CheckHopBaFeAs1L10d(norb,model)
  use standard_derived_types 
  implicit none

  integer norb
  type(modelblock) , target :: model
  
  integer norb1atom,atomimage,nbondxy,nbondz
  integer ix,iy,iz,ixg,iyg,izg,ig,ixx,iyy
  integer iorb,jorb,iatom,jatom,idim,jdim
  integer iorbg,jorbg,iatomg,jatomg,idimg,jdimg
  real*8 r(3),rg(3),rx,ry
  complex*16, pointer :: t(:,:),tg(:,:)

  integer indexgi(2),indexgj(2)
  real*8 Xgi(2),Xgj(2)

  complex*16 tdf(5,5,-6:6,-6:6,-2:2)

  norb1atom=5; 
  nbondxy=7
  nbondz=2

  do while(.false.)
    print*,'idim,jdim,ix,iy,iz:'
	read(*,*)idim,jdim,ix,iy,iz
	print*,model%t(idim,jdim,ix,iy,iz)
  end do
  print*,'Check Hopping symmetry ...'
  do iz=-nbondz,nbondz;  do iy=-nbondxy,nbondxy; do ix=-nbondxy,nbondxy
    t=>model%t(:,:,ix,iy,iz)
	r=(/ix,iy,iz/); r(1:2)=r(1:2)/2
	do ig=2,model%ng
	  call groupaction(ig,r,rg,model%group)
	  !if(ig<=4)rg(3)=-rg(3)
      ixg=nint(rg(1)*2); iyg=nint(rg(2)*2); izg=nint(rg(3))
	  tg=>model%t(:,:,ixg,iyg,izg)
      do iatom=1,model%natom
	    iatomg=atomimage(ig,iatom,model)
		do iorb=1,norb1atom
          idim=iorb+(iatom-1)*norb1atom
		  call orbitimage(ig,iorb,norb1atom,model%orbits,model%group,indexgi,Xgi)
		  iorbg=indexgi(1)
		  idimg=iorbg+(iatomg-1)*norb1atom
	      do jatom=1,model%natom
	        jatomg=atomimage(ig,jatom,model)
		    do jorb=1,norb1atom
              jdim=jorb+(jatom-1)*norb1atom
		      call orbitimage(ig,jorb,norb1atom,model%orbits(1:5),model%group,indexgj,Xgj)
		      jorbg=indexgj(1)
		      jdimg=jorbg+(jatomg-1)*norb1atom

			  if(abs(t(idim,jdim)-tg(idimg,jdimg)*Xgi(1)*Xgj(1))>1.e-4)then
			    print*,r
				print*,rg
			    print*,'  ix ,  iy ,  iz ,  iorb ,iatom ,idim,jorb, jatom, jdim '
				write(*,'(9i6)')ix,iy,iz,iorb,iatom,idim,jorb,jatom,jdim
				print*,'ig:',ig
				print*,'ixg  ,iyg  ,izg  ,iorbg,iatomg,idimg,jorbg,jatomg,jdimg'
				write(*,'(9i6)')ixg,iyg,izg,iorbg,iatomg,idimg,jorbg,jatomg,jdimg
				print*,t(idim,jdim)
				print*,tg(idimg,jdimg),Xgi(1),Xgj(1)
				pause
			  end if
			end do
		  end do
		end do
	  end do
	end do
  end do; end do; end do

  print*,'hoppings satisfy group symmetry'

  do iz=-nbondz,nbondz;  do iy=-nbondxy,nbondxy; do ix=-nbondxy,nbondxy
    t=>model%t(:,:,ix,iy,iz); tg=>model%t(:,:,ix,iy,-iz)
	do iatom=1,model%natom; do iorb=1,norb1atom; idim=iorb+(iatom-1)*norb1atom
	do jatom=1,model%natom; do jorb=1,norb1atom; jdim=jorb+(jatom-1)*norb1atom
	  iatomg=mod(iatom,2)+1; idimg=iorb+(iatomg-1)*norb1atom
	  jatomg=mod(jatom,2)+1; jdimg=jorb+(jatomg-1)*norb1atom
	  if(abs(t(idim,jdim)-tg(idimg,jdimg))>1.e-6)then
	    print*,'  ix , iy , iz'
		write(*,'(3i5)')ix,iy,iz
	    print*,'iatom,iorb,jatom,jorb'
		write(*,'(4i5)')iatom,iorb,jatom,jorb
		print*,t(idim,jdim)
        print*,'iatomg,iorb,jatomg,jorb'
		write(*,'(4i5)')iatomg,iorb,jatomg,jorb
		print*,tg(idimg,jdimg)
	    pause
	  end if
	end do; end do
	end do; end do
  end do; end do; end do

  print*,'Hoppings can be downfold'
  
  tdf=0
  iatom=1
  do iorb=1,norb1atom; idim=iorb
  do iz=-2,2; do iy=-7,7; do ix=-7,7
    if(mod(ix+10,2)==0.and.mod(iy+10,2)==0)then
	  jatom=1
    else if(mod(ix+10,2)/=0.and.mod(iy+10,2)/=0)then
	  jatom=2
	else
	  cycle
	end if
	rx=(ix+iy)
	rx=rx/2
	ry=(ix-iy)
    ry=ry/2
	ixx=nint(rx)
	iyy=nint(ry)
	if(abs(ixx)>6.or.abs(iyy)>6)cycle
	do jorb=1,norb1atom; jdim=jorb+(jatom-1)*norb1atom
	  tdf(iorb,jorb,ixx,iyy,iz)=model%t(idim,jdim,ix,iy,iz)
	  if(jatom==2.and.(jorb==2.or.jorb==3))tdf(iorb,jorb,ixx,iyy,iz)=tdf(iorb,jorb,ixx,iyy,iz)
	end do
  end do; end do; end do
  end do

  open(10,file='BaFeAs1L5d.model')
    write(10,*)tdf
  close(10)



  return
end subroutine

subroutine BasisRotateBaFeAs1L10d(norb,model)
  use standard_derived_types
  implicit none

  integer norb
  type(modelblock) :: model
  complex*16 U(10,10),V(10,10)
  integer i,ix,iy,iz,nbondxy,nbondz
  complex*16 one

  one=cmplx(0.0,1.0)
  nbondxy=7
  nbondz=2
  U=0
  do i=1,10
    U(i,i)=1
  end do
  !U(2,2)=one; U(3,3)=one; U(7,7)=one; U(8,8)=one
  !U(7,7)=1; U(7,8)=-1; U(8,7)=1; U(8,8)=1; U(7:8,7:8)=U(7:8,7:8)/sqrt(2.d0); U(2:3,2:3)=U(7:8,7:8)
  U(7,7)=-1; U(8,8)=-1

  V=transpose(U)
  do iz=-nbondz,nbondz; do iy=-nbondxy,nbondxy; do ix=-nbondxy,nbondxy
    model%t(:,:,ix,iy,iz)=matmul(U, matmul(model%t(:,:,ix,iy,iz),V))
  end do; end do; end do

  return
end subroutine


  
	

