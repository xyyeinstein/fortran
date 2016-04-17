subroutine meshqvec()
  use workplace, only : nq,qv,qmesh,model,uniformqmesh
  use standard_derived_types
  implicit none

  integer ng,ig,ip,iq,nest
  real*8 wmax,wmin,qmin,q(2),S,cut

  interface
    function sinesq(p)
	  real*8 sinesq,p(2)
	end function sinesq
	function hexsinesq(p)
	  real*8 hexsinesq,p(2)
	end function hexsinesq
  end interface

  logical hexagonal
  
	if(uniformqmesh)then
	  call meshqvec_uniform()
		return
	end if

  if(model%ladder)then; call meshqvec_ladder(); return; end if

  hexagonal=( abs( sum(model%a*model%b) )> 1.e-5 )

  ng=8; wmax=4; wmin=1.e-3
  if(hexagonal)then
    call HeXBZmesh(ng,wmax,wmin,qmesh,hexsinesq)
  else  
    call BZmesh(ng,wmax,wmin,qmesh,sinesq)
  end if

  nest=model%nest
  nq=qmesh%npoints+1+nest
  allocate(qv(3,nq)); qv=0

  iq=0
  do ig=1,ng; do ip=1,qmesh%g(ig)%np
     if(qmesh%g(ig)%info(ip)/=0)cycle
	 iq=iq+1
     qv(1:2,iq)=qmesh%g(ig)%p(:,ip)
	 qv(3,iq)=qmesh%g(ig)%area
  end do; end do
  if(abs( sum(qv(3,:))-1 )>1.e-5)stop 'jacobi error @ meshqvec'

  qv(:,nq-nest)=0    !zone center

  qv(1:2,nq-nest+1:nq)=model%qnest(:,1:nest)
  qv(3,nq-nest+1:nq)=0

  print*,'# qmesh points = ',nq

  open(10,file='qv.dat')
  write(10,100)qv; close(10)
100 format(1x,3f20.10)

  return
end subroutine meshqvec
  
subroutine meshqvec_uniform()
  use workplace, only : nq,qv,qmesh,model
  use standard_derived_types
  implicit none

  integer nl,ix,iy,iq
	real*8 qx,qy

  nl=40
	nq=nl*nl+1
	allocate(qv(3,nq))
	qv=0
	iq=0
	do ix=1,nl
	do iy=1,nl
	  qx=-1.d0+(2.d0*ix-1.d0)/nl
		qy=-1.d0+(2.d0*iy-1.d0)/nl
		iq=iq+1
		qv(1:2,iq)=(/qx,qy/)
	end do
	end do
	qv(3,:)=1.d0/nq
	qv(:,nq)=0.d0

  print*,'# qmesh points = ',nq

  open(10,file='qv.dat')
  write(10,100)qv; close(10)
100 format(1x,3f20.10)
 
  return
end subroutine
	



function sinesq(p)
  implicit none
  real*8 sinesq,p(2)
  integer i,j
  real*8 pi; data pi/3.1415926/

  real*8 sinesq_
  
  sinesq=sinesq_(p); return

  sinesq=0  
  do i=-1,1; do j=-1,1
     sinesq=sinesq + sum( sin(  ( p-(/i,j/) )*pi  )**2   ) 
   end do; end do
   return
end function sinesq

function sinesq_(p)
  use workplace,only : model
  implicit none
  real*8 sinesq_,p(2)
  integer i,j,n
  real*8 pi; data pi/3.1415926/
  n=model%n4logqmesh
  sinesq_=0  
  do i=-n,n; do j=-n,n
     sinesq_=sinesq_ + sum( sin( n* ( p-(/i,j/)*1./n )*pi  )**2   ) 
   end do; end do
   return
end function sinesq_

function hexsinesq(p)
  implicit none
  real*8 hexsinesq,p(2)
  integer i,j
  real*8 pi; data pi/3.1415926/
  real*8 a(2,6)
  real*8 hexsinesq_

  !hexsinesq=hexsinesq_(p); return
  
  a(:,1)=(/1,0/)                                   !for important spots around (0,2/sqrt(3)) and its symmetric points
  !a(:,1)=(/0,1/)*sqrt(3.)/2 !.5,0.5*sqrt(3.)/)    !for important spots around (4/3,0) and its symmetric points

  do i=1,5; a(:,i+1)=a(:,i); call rotate60(a(:,i+1)); end do

  hexsinesq=0
  do i=1,6
     hexsinesq=hexsinesq + sin(sum(p*a(:,i))*pi)**2   
  end do
  return
end function hexsinesq

function hexsinesq_(p)
  implicit none
  real*8 hexsinesq_,p(2)
  integer i,j
  real*8 pi; data pi/3.1415926/
  
  real*8 a(2),b(2)
  a=(/1,0/)                                   !for important spots around (0,2/sqrt(3)) and its symmetric points
  b=(/0,1/)*sqrt(3.)/2                        !for important spots around (4/3,0) and its symmetric points

  hexsinesq_=0  
  do i=1,6; if(i>1)then; call rotate60(a); call rotate60(b); end if
     hexsinesq_=hexsinesq_ + ( sin(sum(p*a)*pi) * sin(sum(p*b)*pi)   )**2   
  end do
  return
end function hexsinesq_


subroutine meshqvec_ladder()
  use standard_derived_types
  use workplace, only : nq,qv,model
  implicit none

  integer ng,ig,ip,nest
  real*8 wmax,wmin
  type (meshconfig), dimension (2*model%nleg) :: line

  integer iq,iqy,nqy,icount
  real*8 qy,qL(2),qU(2),vq(2),dq

  nqy=model%nleg; if(nqy==0)stop 'invalid # of legs in the ladder model @ meshqvec_ladder'

  ng=15; wmax=6; wmin=1.e-5; dq=2./nqy
  
  icount=0
  do iqy=0,nqy-1; qy=iqy*dq; if(qy>1)qy=qy-2
     qL=(/0.d0,qy/); qU=(/1.d0,qy/); icount=icount+1; call linemesh(ng,wmax,wmin,qL,qU,line(icount),sinexsq)
     qL=(/-1.d0,qy/); qU=(/0.d0,qy/); icount=icount+1; call linemesh(ng,wmax,wmin,qL,qU,line(icount),sinexsq)
  end do; if(icount/=2*nqy)stop '# of linemeshes /= 2*nleg @ meshqvec_ladder'

  nest=model%nest; nq=sum(line(:)%npoints)+1+nest;  allocate(qv(3,nq)); qv=0

  iq=0; icount=0
  do iqy=0,nqy-1; qy=iqy*dq; if(qy>1)qy=qy-2
     qL=(/0.d0,qy/); qU=(/1.d0,qy/); icount=icount+1
     do ig=1,ng; do ip=1,line(icount)%g(ig)%np
       if(line(icount)%g(ig)%info(ip)/=0)cycle
	   iq=iq+1; qv(1:2,iq)=line(icount)%g(ig)%p(1,ip)*(qU-qL)+qL
	   qv(3,iq)=line(icount)%g(ig)%dp
     end do; end do

     qL=(/-1.d0,qy/); qU=(/0.d0,qy/); icount=icount+1
     do ig=1,ng; do ip=1,line(icount)%g(ig)%np
       if(line(icount)%g(ig)%info(ip)/=0)cycle
	   iq=iq+1; qv(1:2,iq)=line(icount)%g(ig)%p(1,ip)*(qU-qL)+qL
	   qv(3,iq)=line(icount)%g(ig)%dp
    end do; end do
  end do
  qv(3,:)=qv(3,:)/sum(qv(3,:))

  qv(1:2,nq-nest+1:nq)=model%qnest(:,1:nest)
  
  open(10,file='qv.dat')
  write(10,100)qv(:,1:nq); close(10)
100 format(1x,3f20.10)

  return
  contains
      function sinexsq(p)
        implicit none
        real*8 sinexsq,p,qx
		qx=qL(1)+(qU(1)-qL(1))*p
        sinexsq=abs(sin(qx*3.1415926))
        return
      end function sinexsq
end subroutine meshqvec_ladder
