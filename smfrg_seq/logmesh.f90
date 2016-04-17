subroutine search(px,ipx,mesh,uniform)
  use standard_derived_types
  implicit none
  type (meshconfig) :: mesh 
  real*8 px(2)
  integer ipx
	logical uniform

  integer ng,ig,ip,ipmax,icount
  
  integer i,ii
  real*8 dist,disti,distance

  if(uniform)then
	  call search_uniform(px,ipx)
		return
	end if

  if(sum(abs(px))<1.e-10)then; ipx=mesh%npoints+1; return; end if

  ng=mesh%ng;  ip=1
  do ig=1,ng
     ii=ip; dist=distance(px-mesh%g(ig)%p(:,ii))
     do i=ip+1,ip+3
        disti=distance(px-mesh%g(ig)%p(:,i))
        if(disti<dist)then; ii=i; dist=disti; end if
     end do
     ip=mesh%g(ig)%info(ii);	 if(ip==0)exit
  end do

  icount=0
  do i=1,ig; ipmax=mesh%g(i)%np; if(i==ig)ipmax=ii
  do ip=1,ipmax; if(mesh%g(i)%info(ip)/=0)cycle
     icount=icount+1
  end do; end do

  ipx=icount
  return
end subroutine search

subroutine search_uniform(px,ipx)
  use workplace, only: qv,nq
  implicit none

	integer ipx
	real*8 px(2)

	integer iq
	real*8 dista,distb,distance

  ipx=1
  dista=100
	do iq=1,nq
	  distb=distance(px-qv(1:2,iq))
    if(distb<dista)then
		  ipx=iq
			dista=distb
		end if
	end do
	return
end subroutine

function distance(p)
  implicit none
  real*8 distance,p(2),pp(2),pi,twopi
  integer i
  real*8 ka(2),kb(2)

  pi=1; twopi=2   !pi set as unit here
  pp=mod(p,twopi)   
  where(pp>pi)pp=pp-twopi; where(pp<-pi)pp=pp+twopi

  distance=sum(abs(pp))
  return
end function distance

subroutine BZmesh(ng,wmax,wmin,mesh,epfunc)
  use standard_derived_types
  implicit none
  integer ng,nmesh
  real*8 wmax,wmin
  type (meshconfig) :: mesh
  interface
    function epfunc(p)
	  real*8 epfunc,p(2)
	end function epfunc
  end interface

  real*8 pc(2),dp
  real*8 w,b,pi,twopi
  integer i,ip,ii,ig,np,icount,ix,iy,nx,nsearch
  real*8 ep
  logical expand

  b=exp(log(wmax/wmin)/ng)

  mesh%ng=ng; allocate(mesh%g(ng))
  mesh%g(:)%np=0   !initialize mesh

  pi=1; twopi=2
  !notice: pi is set to unity in this subroutine

  nx=2; np=nx*nx; dp=pi/2    
  mesh%g(1)%np=np; mesh%g(1)%w=wmax; mesh%g(1)%dp=dp; mesh%g(1)%area=dp*dp

  allocate(mesh%g(1)%info(np),mesh%g(1)%p(2,np)); mesh%g(1)%info=0

  icount=0
  do iy=1,nx; do ix=1,nx; icount=icount+1
     pc=(/ix-nx/2-0.5,iy-nx/2-0.5/)*twopi/nx
     mesh%g(1)%p(:,icount)=pc
  end do; end do
  
  ig=1; w=wmax; nsearch=10
  do while(w>=wmin)
     icount=0; dp=mesh%g(ig)%dp/2; np=mesh%g(ig)%np
     do i=1,np; pc=mesh%g(ig)%p(:,i); expand=.false.
        do ix=-nsearch,nsearch; do iy=-nsearch,nsearch
           ep= epfunc(pc+(/ix,iy/)*dp*2./nsearch) 
           if(abs(ep)<=w)expand=.true.; if(expand)exit
        end do; if(expand)exit; end do	     
        if(.not.expand)cycle
        mesh%g(ig)%info(i)=icount+1
        icount=icount+4
     end do

     if(icount==0)then; w=w/b; cycle; end if
  
     if(ig+1>ng)then; mesh%g(ig)%info=0; exit; end if

     mesh%g(ig+1)%np=icount; mesh%g(ig+1)%dp=dp; mesh%g(ig+1)%area=dp*dp; mesh%g(ig+1)%w=w
     allocate(mesh%g(ig+1)%info(icount),mesh%g(ig+1)%p(2,icount))
     mesh%g(ig+1)%info=0

     icount=0
     do i=1,np; if(mesh%g(ig)%info(i)==0)cycle; pc=mesh%g(ig)%p(:,i)
        do ix=-1,1,2; do iy=-1,1,2; icount=icount+1
           mesh%g(ig+1)%p(:,icount)=pc+(/ix,iy/)*dp
        end do; end do	  
     end do

     w=w/b; ig=ig+1
  end do

  icount=0
  do ig=1,ng; do ip=1,mesh%g(ig)%np; if(mesh%g(ig)%info(ip)/=0)cycle
     icount=icount+1
  end do; end do
  mesh%npoints=icount

  return
end subroutine BZmesh



function area(r)
  implicit none
  real*8 area,r(2,3)
  real*8, dimension (2) :: a,b

  a=r(:,2)-r(:,1); b=r(:,3)-r(:,1)
  area=a(1)*b(2)-a(2)*b(1)
  area=abs(area)

  return
end function area


subroutine hexsearch(px,ipx,mesh)
  use standard_derived_types
  implicit none
  type (meshconfig) :: mesh 
  real*8 px(2)
  integer ipx

  integer ng,ig,ip,ipmax,icount
  
  integer i,ii,n
  real*8 dist,disti,hexdistance

  if(sum(abs(px))<1.e-10)then; ipx=mesh%npoints+1; return; end if

  ng=mesh%ng;  ip=1
  do ig=1,ng
     ii=ip; dist=hexdistance(px-mesh%g(ig)%p(:,ii))
	 if(ig==1)then; n=5; else; n=3; end if
     do i=ip+1,ip+n
        disti=hexdistance(px-mesh%g(ig)%p(:,i))
        if(disti<dist)then; ii=i; dist=disti; end if
     end do
     ip=mesh%g(ig)%info(ii);	 if(ip==0)exit
  end do

  icount=0
  do i=1,ig; ipmax=mesh%g(i)%np; if(i==ig)ipmax=ii
  do ip=1,ipmax; if(mesh%g(i)%info(ip)/=0)cycle
     icount=icount+1
  end do; end do

  ipx=icount
  return
end subroutine hexsearch


function hexdistance(p)
  implicit none
  real*8 hexdistance,p(2),pp(2)
  integer i
  real*8 a(3),b(3)

  a=(/1,0,0/); b=(/0.5,0.5*sqrt(3.),0./)
  call RBZkv(p,a,b,pp)
  hexdistance=sum(pp*pp)
  return
end function hexdistance


subroutine HexBZmesh(ng,wmax,wmin,mesh,epfunc)
  use standard_derived_types
  implicit none
  integer ng,nmesh
  real*8 wmax,wmin
  type (meshconfig) :: mesh
  interface
    function epfunc(p)
	  real*8 epfunc,p(2)
	end function epfunc
  end interface

  real*8 pc(2),ps(2,200)
  real*8 w,b
  integer i,j,ip,ig,np,icount,is,nps
  real*8 ep
  logical expand,mannual

  real*8 a,vertex(2,3),vertices(2,3,4),Stot,area
  real*8 point(2)

  b=exp(log(wmax/wmin)/ng)

  mesh%ng=ng; allocate(mesh%g(ng))
  mesh%g(:)%np=0   !initialize mesh

  np=6; mesh%g(1)%np=np; mesh%g(1)%w=wmax

  allocate(mesh%g(1)%info(np),mesh%g(1)%p(2,np),mesh%g(1)%tri(2,3,np)); mesh%g(1)%info=0

  vertex(:,1)=0; vertex(:,2)=(/4./3,0./); vertex(:,3)=(/2./3,2./sqrt(3.)/)
  mesh%g(1)%tri(:,:,1)=vertex; Stot=area(vertex)*6;  mesh%g(1)%area=1./6

  do ip=2,np; do i=1,3; call rotate60(vertex(:,i)); end do
     mesh%g(1)%tri(:,:,ip)=vertex
  end do

  do ip=1,np; pc=0
     do i=1,3; pc=pc+mesh%g(1)%tri(:,i,ip)/3; end do
	 mesh%g(1)%p(:,ip)=pc
  end do

  ig=1; w=wmax; mannual=.false.
  do while(w>=wmin)
     icount=0; np=mesh%g(ig)%np
     do ip=1,np; vertex=mesh%g(ig)%tri(:,:,ip); expand=.false.
        if(mannual)then; call searchpoints_mannual(vertex,16,ps); nps=16
		else;  call searchpoints_auto(vertex,14,nps,ps); end if
		do is=1,nps; ep= epfunc(ps(:,is)) 
           if(abs(ep)<w)expand=.true.; if(expand)exit
        end do
        if(.not.expand)cycle
        mesh%g(ig)%info(ip)=icount+1
        icount=icount+4
     end do

     if(icount==0)then; w=w/b; cycle; end if
  
     if(ig+1>ng)then; mesh%g(ig)%info=0; exit; end if

     mesh%g(ig+1)%np=icount; mesh%g(ig+1)%w=w
     allocate(mesh%g(ig+1)%info(icount),mesh%g(ig+1)%p(2,icount),mesh%g(ig+1)%tri(2,3,icount))
     mesh%g(ig+1)%info=0

     icount=0
     do ip=1,np; if(mesh%g(ig)%info(ip)==0)cycle
	    vertex=mesh%g(ig)%tri(:,:,ip); call newvertices(vertex,vertices)
        do i=1,4; icount=icount+1
           mesh%g(ig+1)%tri(:,:,icount)=vertices(:,:,i)	   
		   pc=0
		   do j=1,3; pc=pc+vertices(:,j,i)/3; end do
		   mesh%g(ig+1)%p(:,icount)=pc
        end do
     end do
	 mesh%g(ig+1)%area=area(mesh%g(ig+1)%tri(:,:,1))/Stot

     w=w/b; ig=ig+1
  end do

  icount=0
  do ig=1,ng; do ip=1,mesh%g(ig)%np; if(mesh%g(ig)%info(ip)/=0)cycle
     icount=icount+1
  end do; end do
  mesh%npoints=icount

  return
end subroutine HexBZmesh


subroutine searchpoints_mannual(vertex,n,p)
  implicit none
  real*8 vertex(2,3)
  integer n
  real*8 p(2,*)

  real*8 v1(2),p1(2,15)
  integer i,j,icount
 
  if(n/=16)stop 'n/=16 @ searchpoints'

  icount=1; p1(:,1)=vertex(:,1)
  do i=1,4; v1=vertex(:,1)+(vertex(:,2)-vertex(:,1))*i/4.
     do j=0,i; icount=icount+1
        p1(:,icount)=(vertex(:,3)-vertex(:,2))*j/4.+v1
	 end do
  end do

  p(:,1)=(p1(:,1)+p1(:,2)+p1(:,3))/3
  p(:,2)=(p1(:,2)+p1(:,4)+p1(:,5))/3
  p(:,3)=(p1(:,5)+p1(:,2)+p1(:,3))/3
  p(:,4)=(p1(:,3)+p1(:,4)+p1(:,6))/3
  p(:,5)=(p1(:,4)+p1(:,7)+p1(:,8))/3
  p(:,6)=(p1(:,8)+p1(:,4)+p1(:,5))/3
  p(:,7)=(p1(:,5)+p1(:,8)+p1(:,9))/3
  p(:,8)=(p1(:,9)+p1(:,5)+p1(:,6))/3
  p(:,9)=(p1(:,6)+p1(:,9)+p1(:,10))/3
  p(:,10)=(p1(:,7)+p1(:,11)+p1(:,12))/3
  p(:,11)=(p1(:,12)+p1(:,7)+p1(:,8))/3
  p(:,12)=(p1(:,8)+p1(:,12)+p1(:,13))/3
  p(:,13)=(p1(:,13)+p1(:,8)+p1(:,9))/3
  p(:,14)=(p1(:,9)+p1(:,13)+p1(:,14))/3
  p(:,15)=(p1(:,14)+p1(:,9)+p1(:,10))/3
  p(:,16)=(p1(:,10)+p1(:,14)+p1(:,15))/3

  return
end subroutine searchpoints_mannual
  

subroutine searchpoints_auto(vertex,n,np,p)
  implicit none
  real*8 vertex(2,3)
  integer n,np
  real*8 p(2,*)

  real*8 a(2),b(2),x(2),y(2),w(2),c(2),wx,wy,sqrt3
  integer ia,ib
 
  a=vertex(:,3)-vertex(:,2); b=vertex(:,1)-vertex(:,2); c=vertex(:,3)
  x=a/sqrt(sum(a*a)); y= vertex(:,1)-(vertex(:,2)+vertex(:,3))/2; y=y/sqrt(sum(y*y))
  if(abs(sum(x*y))>1.e-5)stop 'unit vectors not orthorgonal @ searchpoints_'

  np=0; sqrt3=sqrt(3.)
  do ib=1,n; do ia=1,n
     w = vertex(:,2) + (ia-0.25)*a/n + (ib-0.25)*b/n
	 wx=sum( (w-c)*x )
	 wy=sum( (w-c)*y ); if(wy<0)stop 'wy can not be negative here'
	 if( wy>=-wx*sqrt3)cycle
	 np=np+1
	 p(:,np)=w
  end do; end do    

  do ib=1,n; do ia=1,n
     w = vertex(:,2) + (ia-0.75)*a/n + (ib-0.75)*b/n
	 wx=sum( (w-c)*x )
	 wy=sum( (w-c)*y ); if(wy<0)stop 'wy can not be negative here'
	 if( wy>=-wx*sqrt3)cycle
	 np=np+1
	 p(:,np)=w
  end do; end do    

  return
end subroutine searchpoints_auto


subroutine newvertices(vertex,vertices)
  implicit none
  real*8 vertex(2,3),vertices(2,3,4)
  
  !purpose: break a triangle into 4 small triangles

  real*8 v1(2),v2(2),p(2,6)


  v1=vertex(:,2)-vertex(:,1); v2=vertex(:,3)-vertex(:,2)
  v1=v1/2; v2=v2/2

  p(:,1)=vertex(:,1)
  p(:,2)=p(:,1)+v1
  p(:,3)=p(:,2)+v2
  p(:,4)=p(:,2)+v1
  p(:,5)=p(:,4)+v2
  p(:,6)=p(:,5)+v2

  vertices(:,:,1)=p(:,1:3)
  
  vertices(:,1,2)=p(:,2)
  vertices(:,2:3,2)=p(:,4:5)

  vertices(:,1,3)=p(:,3)
  vertices(:,2:3,3)=p(:,5:6)

  vertices(:,1,4)=p(:,5)
  vertices(:,2:3,4)=p(:,2:3)

  return
end subroutine newvertices



subroutine linemesh(ng,wmax,wmin,kstart,kend,line,epfunc)
  use standard_derived_types
  implicit none

  integer ng,nl
  real*8 wmax,wmin
	real*8, dimension (2) :: kstart,kend
  type (meshconfig) :: line
    interface
    function epfunc(p)
	    real*8 epfunc,p
    end function epfunc
  end interface

  real*8 pc,dp,ltot
  real*8 w,b,ep
  integer i,ip,ii,ig,np,icount,ix,nx,nsearch
  logical expand
    
  b=exp(log(wmax/wmin)/ng)

  line%ng=ng; allocate(line%g(ng)); line%g(:)%np=0

  nx=2; np=nx; dp=0.25
  line%g(1)%np=np; line%g(1)%w=wmax; line%g(1)%dp=dp

  allocate(line%g(1)%info(np),line%g(1)%p(1,np)); line%g(1)%info=0

  icount=0
  do ix=1,nx; icount=icount+1
    pc=(ix-nx/2+0.5)/nx
    line%g(1)%p(1,icount)=pc
  end do
  
  ig=1; w=wmax; nsearch=10
  do while(w>=wmin)
    icount=0; dp=line%g(ig)%dp/2; np=line%g(ig)%np
    do i=1,np; pc=line%g(ig)%p(1,i); expand=.false.
      do ix=-nsearch,nsearch
        ep= epfunc(pc+ix*dp*2./nsearch) 
        if(abs(ep)<=w)expand=.true.; if(expand)exit
      end do     
      if(.not.expand)cycle
      line%g(ig)%info(i)=icount+1
      icount=icount+2
    end do

    if(icount==0)then; w=w/b; cycle; end if
  
    if(ig+1>ng)then; line%g(ig)%info=0; exit; end if

    line%g(ig+1)%np=icount; line%g(ig+1)%dp=dp
    line%g(ig+1)%w=w
    allocate(line%g(ig+1)%info(icount),line%g(ig+1)%p(1,icount))
    line%g(ig+1)%info=0
    icount=0
    do i=1,np; if(line%g(ig)%info(i)==0)cycle; pc=line%g(ig)%p(1,i)
      do ix=-1,1,2; icount=icount+1
        line%g(ig+1)%p(1,icount)=pc+ix*dp
      end do 
    end do

    w=w/b; ig=ig+1
  end do
      
  icount=0; ltot=0
  do ig=1,ng; do ip=1,line%g(ig)%np; if(line%g(ig)%info(ip)/=0)cycle
    icount=icount+1; ltot=ltot+line%g(ig)%dp*2
  end do; end do
  line%npoints=icount; if(abs(ltot-1)>1.e-5)stop 'jaccobi error @ linemesh'
  
  return
end subroutine linemesh
