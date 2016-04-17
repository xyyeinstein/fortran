subroutine meshqvec()
  use workplace, only : nq,qv,qmesh
  use standard_derived_types
  implicit none

  integer ng,ig,ip,iq
  real*8 wmax,wmin,qmin,q(2),S,cut

  interface
    function sinesq(p)
	  real*8 sinesq,p(2)
	end function sinesq
  end interface

  ng=8; wmax=4; wmin=1.e-3
  call BZmesh(ng,wmax,wmin,qmesh,sinesq)

  nq=qmesh%npoints+1
  allocate(qv(3,nq)); qv=0

  iq=0
  do ig=1,ng; do ip=1,qmesh%g(ig)%np
     if(qmesh%g(ig)%info(ip)/=0)cycle
	 iq=iq+1
     qv(1:2,iq)=qmesh%g(ig)%p(:,ip)
	 qv(3,iq)=qmesh%g(ig)%dp**2
  end do; end do
  if(sum(qv(3,:))/=1)stop 'jacobi error @ meshqvec'

  q=1.e-5; call search(q,iq,qmesh)
  q=qv(1:2,iq); qmin=sum(q*q)+1.e-6; S=0; cut=0
  do iq=1,nq-1; if(sum(qv(1:2,iq)**2)>qmin)cycle
     qv(1:2,iq)=qv(1:2,iq)*(1+2*cut)
	 S=S+qv(3,iq)*cut
	 qv(3,iq)=qv(3,iq)*(1-cut)
  end do
  qv(1:2,nq)=0; qv(3,nq)=S

  open(10,file='qv.dat')
  write(10,100)qv; close(10)
100 format(1x,3f20.10)

  return
end subroutine meshqvec
   

function sinesq(p)
  implicit none
  real*8 sinesq,p(2)
  integer i,j
  real*8 pi; data pi/3.1415926/
  
  sinesq=0  
  do i=-1,1; do j=-1,1
     sinesq=sinesq + sum( sin(  ( p-(/i,j/) )*pi  )**2   ) 
   end do; end do
   return
end function sinesq
