subroutine Hujp()
  use workplace
  implicit none

  real*8 Q,cutoff
  integer iform,nband
  logical dorpa,appendrpax0

  DataReady=.false.
  dorpa=.false.
  doMF=.false.
  appendrpax0=.false.
  useX0=.false.
  appendX0=.false.
  Appendkmesh=.false.
  appendfs=.false.
  uniformkmesh=.false.
  uniformqmesh=.false.
  QuickSearch=.false.
  AppendBreakPoint=.false.
  kmeshfile='hujpkmesh.dat'
  useprojectors=.false.; appendprojectors=.false.; projectorfile='hujpproj.dat'
  useprojectors=useprojectors.and.(.not.QuickSearch)

  norb=4; natom=2; nambu=8; model%n4logqmesh=2
  model%norb=norb; model%natom=natom; model%sample='hujp'
  print*,'The sample chosen is ',model%sample 

  allocate(model%ra(3,natom)); model%ra(:,1)=(/0.,0.,0./); model%ra(:,2)=(/0.5,0.5,0./)
  model%a=(/1,0,0/); model%b=(/0,1,0/); model%c=(/0,0,1/)
  allocate(model%orbits(norb)); model%orbits=(/'E1u','E2u','E1u','E2u'/)
  model%periodiclayers=.false.
    
  usegroup=.true.;  group='S4_';  ng=4
  model%group=group; model%ng=ng
  
  !open(10,file='frg.input')
  !read(10,*)model%U,model%JH,model%mu,model%t1
  !close(10)
  model%U=6
  model%Jh=1.5
  model%Uab=model%U-2*model%JH
  model%mu=-0.30
  nkx=100
  nky=100
  !call filling2mu(.false.)

  model%propergauge=.false.

  if(.true.)then
    !for KFe2Se2
    nband=2; model%nband=nband; 
    allocate(model%npatch(nband),model%nativeband(nband),model%pocket(nband),model%mushift(nband))
    model%nativeband=(/3,4/); model%pocket=(/'M_','M_'/); model%npatch=240; model%mushift=0.0
  else
    !for BaFe2As2
    nband=4; model%nband=nband; 
    allocate(model%npatch(nband),model%nativeband(nband),model%pocket(nband),model%mushift(nband))
    model%nativeband=(/1,2,3,4/);	model%pocket=(/'G_','G_','M_','M_'/); model%npatch=64; model%mushift=0.
  end if

  Q=0; call checksymmetry(Q); 
  call meshfermisurface(model)
  call fermisurf(norb,Q,model)
  call bandstructure(norb,Q,model) 
  
     
  
  if(dorpa)then
    nLform=1
    Lcontact=0
  else
    nLform=5
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
			subroutine addmfbond(atom,nbatom,basis)
	  		integer atom,nbatom
	  		real*8 basis(2)
	  		if(iform>MF%nbond)return
	  		iform=iform+1
	  		MF%chi(iform)%atom=atom
	  		MF%chi(iform)%nbatom=nbatom
	  		MF%chi(iform)%bond=basis
	  		MF%delta(iform)%atom=atom
	  		MF%delta(iform)%nbatom=nbatom
	  		MF%delta(iform)%bond=basis
    	end subroutine
end subroutine
  
subroutine Hujphk(norb,kv,hk,model,mf,slvbsn)
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

  complex*16 t(norb,norb,-2:2,-2:2),t1,t1g,U(norb,norb),V(norb,norb)
  integer idim,jdim,iorb,jorb,iatom,jatom,ix,iy,iatomg,jatomg,norb1atom,idimg,jdimg,ixg,iyg
  real*8 rx,ry,splita,splitb,splitc

  t4a=-0.0; t4b=-0.0
  !general case
  t1s=model%t1; t1d=-0.03; t2s=0.25; t2d=0.6; t3s=0.05; t3d=-0.05; tc=0.05
  
  splita=0.0; splitb=0.0; splitc=0.

  t1x=t1s+t1d; t1y=t1s-t1d; t2=t2s+t2d; tt2=t2s-t2d; t3x=t3s+t3d; t3y=t3s-t3d
  t=0;
    
  do ix=-2,2; do iy=-2,2
    do idim=1,norb; iatom=(idim-1)/2+1; iorb=idim-2*(iatom-1)
    do jdim=1,norb; jatom=(jdim-1)/2+1; jorb=jdim-2*(jatom-1)
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
      	if(abs(rx)<1.e-5.and.iorb==iatom.and.iatom==1)t(idim,jdim,ix,iy)=tt2
      	if(abs(ry)<1.e-5.and.iorb/=iatom.and.iatom==1)t(idim,jdim,ix,iy)=tt2
      	if(abs(rx)<1.e-5.and.iorb/=iatom.and.iatom==1)t(idim,jdim,ix,iy)=t2
          
      	if(abs(ry)<1.e-5.and.iorb==iatom.and.iatom==2)t(idim,jdim,ix,iy)=tt2
      	if(abs(rx)<1.e-5.and.iorb==iatom.and.iatom==2)t(idim,jdim,ix,iy)=t2
      	if(abs(ry)<1.e-5.and.iorb/=iatom.and.iatom==2)t(idim,jdim,ix,iy)=t2
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
  
  vk=kv*3.1415926
  
  hk=0; one=cmplx(0,1)
    
  do ix=-2,2; do iy=-2,2
     hk=hk+t(:,:,ix,iy)*exp(one*(vk(1)*ix+vk(2)*iy))
  end do; end do
  
  if(model%propergauge.eqv..false.)then
    hk(1:2,3:4)=hk(1:2,3:4)*exp(one*(vk(1)*0.5+vk(2)*0.5))
    hk(3:4,1:2)=hk(3:4,1:2)*exp(-one*(vk(1)*0.5+vk(2)*0.5))
  end if

  do iorb=1,norb; hk(iorb,iorb)=hk(iorb,iorb)-model%mu; end do

  splita=splita*cos(vk(1))*cos(vk(2))
  hk(1,3)=hk(1,3)+splita; hk(3,1)=hk(3,1)+splita
  hk(2,4)=hk(2,4)+splita; hk(4,2)=hk(4,2)+splita

  splitb=splitb*sum(cos(vk))
  hk(1,1)=hk(1,1)+splitb; hk(2,2)=hk(2,2)+splitb
  hk(3,3)=hk(3,3)-splitb; hk(4,4)=hk(4,4)-splitb

  splitc=splitc*sin(0.5*vk(1))*sin(0.5*vk(2))
  hk(1,4)=hk(1,4)+splitc; hk(4,1)=hk(4,1)+splitc
  hk(2,3)=hk(2,3)+splitc; hk(3,2)=hk(3,2)+splitc
  
  !print*,sum(abs(imag(hk)))
  return
end subroutine
