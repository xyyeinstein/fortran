function searchtable(kv,nkx,nky,table)
  implicit none
  real*8 searchtable
  integer nkx,nky
  real*8 kv(2),table(-nkx:nkx,-nky:nky)
  
  real*8 kvv(2),t(4),d(4),v(2,4),tL,tU,x,y
  integer i,xL,xU,yL,yU
  complex*16 zx,z(4)
  real*8 spline


  kvv=mod(kv,2.d0); where(kvv>1)kvv=kvv-2; where(kvv<-1)kvv=kvv+2
  
  call findlimits(xL,xU,nkx,kvv(1))
  call findlimits(yL,yU,nky,kvv(2))

  if(xL==xU.and.yL==yU)then
     searchtable=table(xL,yL)
  else if(xL==xU.and.yL/=yU)then
    searchtable = table(xL,yL) + ( table(xL,yU)-table(xL,yL) )*(kvv(2)*nky-yL)/(yU-yL)
  else if(xL/=xU.and.yL==yU)then
    searchtable = table(xL,yL) + ( table(xU,yL)-table(xL,yL) )*(kvv(1)*nkx-xL)/(xU-xL)
  else
    !the last possibility is xL/=xU, yL/=yU, which will be handled by complex spline

    !simple linear interpolation
	x=kvv(1)*nkx; y=kvv(2)*nky
    tL=table(xL,yL); tU=table(xU,yL)
	t(1)=tL+(tU-tL)*(x-xL)/(xU-xL)
    tL=table(xL,yU); tU=table(xU,yU)
	t(2)=tL+(tU-tL)*(x-xL)/(xU-xL)
	tL=t(1); tU=t(2)
	searchtable=tL+(tU-tL)*(y-yL)/(yU-yL)

    !spline interpolation
    !zx=cmplx(kvv(1)*nkx,kvv(2)*nky)
    !z(1)=cmplx(xL,yL); z(2)=cmplx(xL,yU); z(3)=cmplx(xU,yL); z(4)=cmplx(xU,yU)
    !t(1)=table(xL,yL); t(2)=table(xL,yU); t(3)=table(xU,yL); t(4)=table(xU,yU)
    !searchtable=spline(4,z,t,zx)
  end if
  return
end function searchtable

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


