subroutine KFeSe2L20d()
  use workplace
  implicit none

  real*8 Q
  integer error

  DataReady=.false.
  AppendBreakPoint=.false.
  Appendkmesh=.false.; kmeshfile='KFS20dkmesh.dat'

    norb=20; natom=4    
    model%norb=norb; model%natom=natom; model%sample='KFS20d'
    print*,'The sample chosen is ',model%sample 

	allocate(model%ra(3,natom))
	model%ra(:,1)=(/0.5,0.,-0.25/); model%ra(:,2)=(/0.,0.5,-0.25/)
	model%ra(:,3)=(/0.,0.5,0.25/);  model%ra(:,4)=(/0.5,0.,0.25/)
	model%a=(/1,0,0/); model%b=(/0,1,0/); model%c=(/0,0,1/)
	model%periodiclayers=.true.

    allocate(model%orbits(norb))
	model%orbits(1:5)=(/'A1g','E1g','E2g','B1g','B2g'/)
	model%orbits(6:10)=model%orbits(1:5)  
	model%orbits(11:20)=model%orbits(1:10)
  
    usegroup=.true.;  group='D4h';  ng=8
    model%group=group; model%ng=ng

	model%U=4/2; model%JH=0.54/2; model%Uab=model%U-2*model%JH; model%mu=9.27
			
	model%zerotz=.false.
	call gethop20d(model)
	Q=0; call checksymmetry(Q);	call FermiSurf(norb,Q,model) 
	
    nkx=40; nky=40      
    nLform=3; Lcontact=2
    allocate(Lformtable(nLform)); Lformtable=(/'A1g','B1g','A1g'/)
    allocate(Lformbasis(3,nLform));	Lformbasis(:,1)=0  
	Lformbasis(:,2)=(/1,0,0/); Lformbasis(:,3)=(/1,0,0/)
	call setLformfunc(ng)

	call setOformfunc()

    SkipInterOrbitPair=.true.
	guess='A1g';  call setMolecules()

    nw=400; wmax=100
    wmin=1.e-4; wir=1.e-3
    diverge=40

  return
end subroutine KFeSe2L20d

subroutine gethop20d(model)
  use standard_derived_types
  implicit none
  type (modelblock) :: model

  real*8 t(10,10,-4:4,-4:4,-1:1)
  integer a,b,ix,iy,iz,error
  real*8 x,y,z,tt
  real*8 ra(3,4),bond(3)
  integer Fe(4),iFe,jFe
  integer idim,iidim
  integer jdim,jjdim
  integer iatom,jatom
  integer iorb,jorb
  integer findatom

  open(10,file='KFeSe10d.model')
  t=0
  do while (.true.)
     read(10,*,iostat=error)a,b,x,y,z,tt
     ix=nint(x*2); iy=nint(y*2); iz=nint(z*2)
     t(a,b,ix,iy,iz)=tt
     if(error/=0)exit
  end do; close(10)

  ra=model%ra
  Fe=(/1,2,1,2/)

  allocate(model%t(20,20,-4:4,-4:4,-1:1))
  model%t=0
  do jatom=1,4; jFe=Fe(jatom)
  do iz=-1,1; do iy=-4,4; do ix=-4,4; bond=(/ix,iy,iz/)/2.
  iatom=findatom(ra(:,jatom)+bond); if(iatom<=0)cycle; iFe=Fe(iatom)
  do jorb=1,5; jdim=jorb+(jatom-1)*5; jjdim=jorb+(jFe-1)*5
  do iorb=1,5; idim=iorb+(iatom-1)*5; iidim=iorb+(iFe-1)*5
     model%t(idim,jdim,ix,iy,iz)=t(iidim,jjdim,ix,iy,iz)
  end do; end do
  end do; end do; end do
  end do
  
  return
end subroutine gethop20d

subroutine KFeSe2L20dhk(norb,kv,hk,model,mf,slvbsn)
  use standard_derived_types
  implicit none
  integer norb
  real*8 kv(2)
  complex*16 hk(norb,norb)
  type (modelblock), target :: model
  type (mfblock) :: mf
	logical slvbsn

  real*8 p(3),bond(3),pi
  complex*16 one
  integer ix,iy,iz,iorb
  complex*16, pointer, dimension (:,:) :: t

  pi=asin(1.d0)*2; one=cmplx(0,1)
  p=0; p(1:2)=kv*pi/2

  hk=0
  do iorb=1,20; hk(iorb,iorb)=-model%mu; end do
  do iz=-1,1; do iy=-4,4; do ix=-4,4
     bond=(/ix,iy,iz/)
     t=>model%t(:,:,ix,iy,iz)
	 hk=hk+t*exp( -one*sum(p*bond) )
  end do; end do; end do

  return
end subroutine KFeSe2L20dhk
