subroutine Hujp2()
  use workplace
  implicit none

  real*8 Q,cutoff
  integer iform,nband,icase
  logical dorpa,appendrpax0

  DataReady=.false.
  dorpa=.false.
  appendrpax0=.false.
  useX0=.true.
  uniformqmesh=.true.
  appendX0=.false.
  Appendkmesh=.false.
  appendfs=.false.
  if(.true.)then 
    appendx0=.true.
    appendkmesh=.true.
    appendfs=.true.
  end if
  QuickSearch=.false.
  AppendBreakPoint=.false.
  kmeshfile='hujp2kmesh.dat'
  useprojectors=.false.; appendprojectors=.false.; projectorfile='hujp2proj.dat'
  useprojectors=useprojectors.and.(.not.QuickSearch)

  norb=8; natom=4; nambu=16
  model%norb=norb; model%natom=natom; model%sample='hujp2'
  print*,'The sample chosen is ',model%sample 

  allocate(model%ra(3,natom))
  model%ra(:,1)=(/0.0,0.0,0.0/)
  model%ra(:,2)=(/0.5,0.5,0.0/)
  model%ra(:,3)=(/0.5,0.5,1.0/)
  model%ra(:,4)=(/0.0,0.0,1.0/)

  model%a=(/1,0,0/); model%b=(/0,1,0/); model%c=(/0,0,1/)
  allocate(model%orbits(norb)); model%orbits=(/'E1u','E2u','E1g','E2g','E1u','E2u','E1g','E2g'/)
  model%periodiclayers=.false.
    
  usegroup=.true.;  group='S4_';  ng=4
  model%group=group; model%ng=ng
    
  open(10,file='frg.input')
  read(10,*)model%U,model%Jh,model%mu,model%t1,model%split
  close(10)
  !model%U=2.4
  !model%Jh=0.4
  if(dorpa)then
    print*,'U='
	read(*,*)model%U
	print*,'JH='
	read(*,*)model%JH
  end if
  model%Uab=model%U-2*model%JH
  !model%mu=0

  model%propergauge=.false.

  icase=1
  select case(icase)
    case(1)
      !for KFe2Se2
      nband=4; model%nband=nband; !model%mu=-0.30
      allocate(model%npatch(nband),model%nativeband(nband),model%pocket(nband),model%mushift(nband))
      model%nativeband=(/5,6,7,8/)
	  model%pocket=(/'M_','M_','M_','M_'/)
	  model%npatch=240
	  model%mushift=0.
    case(2)
      !for BaFe2As2
      nband=8; model%nband=nband; model%mu=-.30
      allocate(model%npatch(nband),model%nativeband(nband),model%pocket(nband),model%mushift(nband))
      model%nativeband=(/1,2,3,4,5,6,7,8/)
	  model%pocket=(/'G_','G_','G_','G_','M_','M_','M_','M_'/)
	  model%npatch=60
	  model%mushift=0.
	case(3)
      !for KFe2As2
      nband=4; model%nband=nband; model%mu=-0.76
      allocate(model%npatch(nband),model%nativeband(nband),model%pocket(nband),model%mushift(nband))
      model%nativeband=(/1,2,3,4/)
	  model%pocket=(/'G_','G_','G_','G_'/)
	  model%npatch=60
	  model%mushift=0.
  end select


  Q=0; call checksymmetry(Q); 
  call fermisurf(norb,Q,model)
  call bandstructure(norb,Q,model) 
  call meshfermisurface(model)


  nkx=320
  nky=nkx    
  
  if(dorpa)then
    nLform=1
    Lcontact=0
  else
    nLform=7
    Lcontact=2
  end if
  allocate(Lformbasis(3,nLform),Lformtable(nLform))

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
  call addform( (/0,0,1/)*1.d0,'A1g')
  call addform( (/0,0,-1/)*1.d0,'A1g')

  call setLformfunc(ng)  !set fundamental lattice forms instead of irreducible harmonics

  skipOddOrbPair=.true.
  call setOformfunc()

  guess='***';  call setMolecules()

  nw=400
  wmax=1.e2
  wmin=1.e-4
  wir=1.e-5
  wirx0=1.e-5
  diverge=100

  if(.not.dorpa)return

  call writeinputs()
  call setupworkplace()
  cutoff=1.e-3; call rpa(cutoff,appendrpax0,'hujpX0.dat'); stop

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
  
subroutine Hujp2hk(norb,kv,hk,model,mf,slvbsn)
  use standard_derived_types
  implicit none
  integer norb
  real*8 kv(2)
  complex*16 hk(norb,norb),one
  type (modelblock) :: model
  type (mfblock) :: mf
	logical slvbsn
 
  real*8 t1s,t1d,t2s,t2d,t3s,t3d
  real*8 t1x,t1y,t2,tt2,t3x,t3y,tc,t4a,t4b
  real*8 f1s,f1d,f2s,f2d,f3s,f3d
  
  real*8 vk(2)

  complex*16 t(norb,norb,-2:2,-2:2),V(4,4)
  integer idim,jdim,iorb,jorb,iatom,jatom,ix,iy,norb1layer
  real*8 rx,ry,split

  t4a=-0.0; t4b=-0.0
  !general case
  t1s=model%t1; t1d=-0.0; t2s=0.25; t2d=0.6; t3s=0.05; t3d=-0.00; tc=0.05
  t1d=-0.03; t3d=-0.05
  split=model%split
  
  !FeAs
  !t1s=0.3; t1d=0.1; t2s=0.3; t2d=1.0; t3s=0.; t3d=0.0; tc=0.05
    
  !FeSe
  !t1s=0.45; t1d=0.05; t2s=0.3; t2d=1.0; t3s=0; t3d=0; tc=0.05
  
  t1x=t1s+t1d; t1y=t1s-t1d; t2=t2s+t2d; tt2=t2s-t2d; t3x=t3s+t3d; t3y=t3s-t3d
  
  t=0;
  !t(1,1,0,0)=0.1
  !t(2,2,0,0)=0.1
  !t(3,3,0,0)=-0.1
  !t(4,4,0,0)=-0.1
  
  norb1layer=4
  do ix=-2,2; do iy=-2,2
     do idim=1,norb1layer; iatom=(idim-1)/2+1; iorb=idim-2*(iatom-1)
     do jdim=1,norb1layer; jatom=(jdim-1)/2+1; jorb=jdim-2*(jatom-1)
        rx=ix+model%ra(1,jatom)-model%ra(1,iatom)
        ry=iy+model%ra(2,jatom)-model%ra(2,iatom)
        !nn
        if(abs(rx**2+ry**2-0.5)<1.e-5)then
          if(iatom==jatom)then
            print*,'atom position wrong!'
            cycle
          end if
		  if(iorb==jorb)t(idim,jdim,ix,iy)=tc
		  if(iorb/=jorb)then 
		    if(rx*ry>0.and.iorb==iatom)t(idim,jdim,ix,iy)=t1y
            if(rx*ry<0.and.iorb==iatom)t(idim,jdim,ix,iy)=t1x
            if(rx*ry>0.and.iorb/=iatom)t(idim,jdim,ix,iy)=-t1x
            if(rx*ry<0.and.iorb/=iatom)t(idim,jdim,ix,iy)=-t1y
		  end if
        end if
        !nnn
        if(abs(rx**2+ry**2-1)<1.e-5)then
          if(iatom/=jatom)then
            print*,'atom position wrong!'
            cycle
          end if
          if(iorb/=jorb)cycle
          if(abs(ry)<1.e-5.and.iorb==iatom.and.iatom==1)t(idim,jdim,ix,iy)=t2
          if(abs(ry)<1.e-5.and.iorb==iatom.and.iatom==2)t(idim,jdim,ix,iy)=tt2
          if(abs(rx)<1.e-5.and.iorb==iatom.and.iatom==1)t(idim,jdim,ix,iy)=tt2
          if(abs(rx)<1.e-5.and.iorb==iatom.and.iatom==2)t(idim,jdim,ix,iy)=t2

          if(abs(ry)<1.e-5.and.iorb/=iatom.and.iatom==1)t(idim,jdim,ix,iy)=tt2
          if(abs(ry)<1.e-5.and.iorb/=iatom.and.iatom==2)t(idim,jdim,ix,iy)=t2
          if(abs(rx)<1.e-5.and.iorb/=iatom.and.iatom==1)t(idim,jdim,ix,iy)=t2
          if(abs(rx)<1.e-5.and.iorb/=iatom.and.iatom==2)t(idim,jdim,ix,iy)=tt2
        end if
        !3nn
        if(abs(rx**2+ry**2-2)<1.e-5)then
          if(iatom/=jatom)then
            print*,'atom position wrong!'
            cycle
          end if
          if(iorb/=jorb)cycle
          if(rx*ry>0.and.iatom==iorb)t(idim,jdim,ix,iy)=t3y
          if(rx*ry<0.and.iatom==iorb)t(idim,jdim,ix,iy)=t3x
          if(rx*ry>0.and.iatom/=iorb)t(idim,jdim,ix,iy)=t3x
          if(rx*ry<0.and.iatom/=iorb)t(idim,jdim,ix,iy)=t3y
        end if
        !4nn
		if(abs(rx**2+ry**2-2.5)<1.e-5)then
		  if(iatom==jatom)then
		    print*,'atom position wrong!'
			cycle
		  end if
		  if(iorb/=jorb)cycle
		  if(iorb==1)then
		    if(abs(rx)>abs(ry))t(idim,jdim,ix,iy)=t4a
			if(abs(rx)<abs(ry))t(idim,jdim,ix,iy)=t4b
		  end if
		  if(iorb==2)then
		    if(abs(rx)>abs(ry))t(idim,jdim,ix,iy)=t4b
			if(abs(rx)<abs(ry))t(idim,jdim,ix,iy)=t4a
		  end if
		end if
			
		    
     end do; end do
  end do
  end do
  
  V=0
  V(1,1)=1; V(2,2)=1; V(3,3)=1; V(4,4)=-1

  do ix=-2,2; do iy=-2,2
    !t(1:4,1:4,ix,iy)=matmul(V,matmul(t(1:4,1:4,ix,iy),V))
  end do; end do

  !do iorb=1,norb; do jorb=1,norb
  !  do ix=-2,2; do iy=-2,2
  !    if(t(iorb,jorb,ix,iy)/=t(jorb,iorb,-ix,-iy))then
  !      print*,'t is not symmetric'
  !      stop
  !    end if
  !  end do; end do
  !end do; end do
  
 
  vk=kv*3.1415926
  
  hk=0; one=cmplx(0,1)
    
  do ix=-2,2; do iy=-2,2
     hk=hk+t(:,:,ix,iy)*exp(one*(vk(1)*ix+vk(2)*iy))
  end do; end do
  
  hk(1:2,3:4)=hk(1:2,3:4)*exp(one*(vk(1)*0.5+vk(2)*0.5))
  hk(3:4,1:2)=hk(3:4,1:2)*exp(-one*(vk(1)*0.5+vk(2)*0.5))

  do iorb=1,norb; hk(iorb,iorb)=hk(iorb,iorb)-model%mu; end do

  hk(5:8,5:8)=hk(1:4,1:4)

  hk(1,7)=split
  hk(2,8)=split
  hk(3,5)=split
  hk(4,6)=split
  hk(5,3)=split
  hk(6,4)=split
  hk(7,1)=split
  hk(8,2)=split

  do iorb=5,8
    hk(iorb,iorb)=hk(iorb,iorb) !+1.e-3
  end do

  return
end subroutine
