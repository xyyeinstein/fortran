

subroutine findlimits(L,U,K,x)
  integer L,U,K
  real*8 x,y

  y=x*K; L=max(int(y-0.5),-K); U=min(int(y+0.5),K)
  if(y>=L+1)L=L+1; if(y<=U-1)U=U-1
  do while(U-L>1)
     M=(L+U)/2
	 if(y>M)then
	   L=M
	 else if(y<M)then
	   U=M
	 else
	   L=M; U=M; return
	 end if
   end do
   return
end subroutine findlimits  

function spline(n,zn,tn,z)
  implicit none
  real*8 spline
  integer n
  complex*16 zn(n),z
  real*8 tn(n)

  integer i,j
  complex*16 sum,factor

  sum=0
  do i=1,n; factor=1
     do j=1,n; if(j==i)cycle
        factor=factor*(z-zn(j))/(zn(i)-zn(j))
     end do
     sum=sum+tn(i)*factor
  end do

  spline=real(sum)

  return
end function spline

subroutine ascendingorder(n,x,info)
  implicit none
  integer n
  real*8 x(n)
  integer info(n)

  real*8 ox(n)
  integer sn(n)
  integer i,L
  integer LowerBound

  ox=x
  do i=1,n; info(i)=i; end do
  sn(1)=1

  do i=2,n
     L=LowerBound(x(i),i-1,ox(1:i-1))+1
	 if(L==i)then
	   ox(i)=x(i); sn(i)=info(i)
	 else
	   ox(L+1:i)=ox(L:i-1); ox(L)=x(i)
	   sn(L+1:i)=sn(L:i-1); sn(L)=info(i)
	 end if
  end do

  info=sn
  x=ox

  return
end subroutine ascendingorder

function LowerBound(x,n,xn)
  integer LowerBound,n
  real*8 x,xn(n)

  integer L,U,M

  if(x<=xn(1))then
    L=0
  else if(x>=xn(n))then
    L=n
  else
    L=1; U=n
    do while (U-L>1)
       M=(L+U)/2
	   if(x<xn(M))then
	     U=M
	   else if(x>xn(M))then
	     L=M
	   else
	     L=M; U=M; exit
	   end if
    end do
  end if
  LowerBound=L
  return
end function LowerBound


