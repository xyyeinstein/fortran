
subroutine EigenBCS_R(Rvmax)
  use workplace
  implicit none

  integer Rvmax

  complex*16, dimension (natom,natom,-Rvmax:Rvmax,-Rvmax:Rvmax) :: PR,CR,DR
  complex*16, dimension (natom,natom,natom,natom,-Rvmax*2:Rvmax*2,-Rvmax*2:Rvmax*2) :: Xpp
  complex*16, allocatable, dimension (:,:) :: A,W(:)

  integer nA,idim,jdim,aatom,batom,catom,a1,b2,ix,iy,jx,jy,imax

  if(ndim/=natom)stop 'EigenBCS_R is designed for onsite forms only'

  call NetPCD_R(PR,CR,DR,Rvmax)

  call getXpp(Xpp,Rvmax*2)
	
  nA=natom*natom*(2*Rvmax+1)**2
  allocate(A(nA,nA),W(nA))
  A=0

  !test the t-U-J case: with strong U, and with varying J1 and J2, there is a transition from d-wave singlet to p-wave triplet. 
  PR=0; CR=0; DR=0; PR(:,:,0,0)=-12
  do ix=-1,1; do iy=-1,1
     if(abs(ix)+abs(iy)==1)CR(:,:,ix,iy)=1
  	 if(abs(ix)+abs(iy)==2)CR(:,:,ix,iy)=-0.5
  end do; end do

  jdim=0; do b2=1,natom; do a1=1,natom; do jy=-Rvmax,Rvmax; do jx=-Rvmax,Rvmax; jdim=jdim+1
  idim=0; do batom=1,natom; do aatom=1,natom; do iy=-Rvmax,Rvmax; do ix=-Rvmax,Rvmax; idim=idim+1

    if(aatom==batom.and.ix==0.and.iy==0)then
	   do catom=1,natom; A(idim,jdim)=A(idim,jdim) + PR(aatom,catom,ix,iy)*Xpp(catom,catom,a1,b2,-jx,-jy); end do
	end if

	A(idim,jdim)=A(idim,jdim)+CR(aatom,batom,ix,iy)*Xpp(batom,aatom,a1,b2,-jx-ix,-jy-iy)
	A(idim,jdim)=A(idim,jdim)+DR(aatom,batom,ix,iy)*Xpp(aatom,batom,a1,b2,-jx+ix,-jy+iy)

  end do; end do; end do; end do
  end do; end do; end do; end do

  call ZGEIGEN(nA,A,W)

  imax=1; wmax=real(W(1))
  do idim=2,nA; if(real(W(idim))>wmax)then; wmax=real(W(idim)); imax=idim; end if; end do

  open(10,file='BCSeigen.dat')

  do jdim=1,nA; if(real(W(jdim))<=0)cycle
     write(10,*)'lambda=',W(jdim)
	 idim=0; do batom=1,natom; do aatom=1,natom; do iy=-Rvmax,Rvmax; do ix=-Rvmax,Rvmax; idim=idim+1
	   if(abs(A(idim,jdim))<1.e-3)cycle
	   write(10,100)aatom,batom,ix,iy,A(idim,jdim)
	 end do; end do; end do; end do
  end do
  close(10)

100 format(1x, 4i6,2f20.10)
  return
end subroutine EigenBCS_R
  



subroutine NetPCD_R(PR,CR,DR,Rvmax)
  use workplace
  implicit none

  integer Rvmax
  complex*16, dimension (ndim,ndim,nq) :: dP,dC,dD
  complex*16, dimension (ndim,ndim,-Rvmax:Rvmax,-Rvmax:Rvmax) :: PR,CR,DR
  integer i1,i2,iq,idim,jdim
  complex*16 factor,factorc
  real*8, dimension (3) :: ri,rj,S

  !subtracting overlapping contributions in C and D from P, so that P+C+D contains no overcounting
  dP=0; dC=C; dD=D; call contact(dP,dC,dD); P=P-dP

  !Fourier transform of P, C and D with respect to the given Rmesh
  PR=0; PR(:,:,0,0)=P(:,:,nq); CR=0; DR=0

  do i2=-Rvmax,Rvmax; do i1=-Rvmax,Rvmax; S=i1*model%a+i2*model%b
	 do iq=1,nq; factor=exp(-one*pi*sum(qv(1:2,iq)*S(1:2)))*qv(3,iq)
     do jdim=1,ndim; rj=model%ra(:,Mformfunc(jdim)%atom )
	 do idim=1,ndim; ri=model%ra(:,Mformfunc(idim)%atom )
	    factorc = exp( -one*pi*sum( qv(1:2,iq)*(rj(1:2)-ri(1:2)) )  )  
	    CR(idim,jdim,i1,i2)=CR(idim,jdim,i1,i2)+factor*C(idim,jdim,iq)*factorc
	    DR(idim,jdim,i1,i2)=DR(idim,jdim,i1,i2)+factor*D(idim,jdim,iq)*factorc
     end do; end do; end do
  end do; end do

  !print*,real(CR(:,:,0,0))
  !print*,real(CR(:,:,0,1))
  !print*,real(CR(:,:,1,1))
  !pause 'CR ok?'
  return
end subroutine NetPCD_R



subroutine getXpp(Xpp,Rvmax)
  use workplace
  implicit none

  integer Rvmax
  complex*16, dimension (natom,natom,natom,natom,-Rvmax:Rvmax,-Rvmax:Rvmax) :: Xpp

  complex*16 chi(norb,norb,norb,norb)
  real*8 kab(2),evk(norb),evkq(norb)
  complex*16, dimension (norb,norb) :: Ak,Akq

  real*8 sk,factor,x1,x2,error,Q
  real*8 fermi

  integer flv1,flv2,flv3,flv4,idim,jdim,iidim,jjdim,m,n,i1,i2,i3,x,y,z,ik

  model%propergauge=.true.

  Xpp=0
  call flexiblemesh(1.e-4)

  !print*,'nk=',nk; pause 'continue?'

  do ik=1,nk
     kab=kv(1:2,ik); sk=kv(3,ik)

	 call gethk(norb,kab,Q,Ak,model,MF,.false.)
	 call ZHEIGEN(norb,Ak,evk)
	 
	 call gethk(norb,-kab,Q,Akq,model,MF,.false.)
	 call ZHEIGEN(norb,Akq,evkq)

     chi=0
	 do m=1,norb; do n=1,norb
	    x1=evk(m); x2=-evkq(n); if(abs(x1-x2)<1.e-8)then; x1=min(x1,x2); x2=x1+1.e-8; end if
		factor=( fermi(x1,1.e-3) - fermi(x2,1.e-3) ) / (x2-x1)
		do flv4=1,norb; do flv3=1,norb;	do flv2=1,norb; do flv1=1,norb
		   chi(flv1,flv2,flv3,flv4)=chi(flv1,flv2,flv3,flv4)+Ak(flv1,m)*conjg(Ak(flv3,m))*Akq(flv2,n)*conjg(Akq(flv4,n))*factor
		end do; end do; end do; end do
	 end do; end do
	 chi=chi*sk

	 do y=-Rvmax,Rvmax; do x=-Rvmax,Rvmax
	    Xpp(:,:,:,:,x,y)=Xpp(:,:,:,:,x,y)+chi*exp( one* ( sum( kab*(x*model%a(1:2)+y*model%b(1:2))) )  )
     end do; end do

  end do


  open(10,file='Xpp.dat')
  write(10,*)Xpp; close(10)

  return
end subroutine getXpp

