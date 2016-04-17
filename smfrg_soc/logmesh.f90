subroutine search(px,ipx,mesh)
  use standard_derived_types
  implicit none
  type (meshconfig) :: mesh 
  real*8 px(2)
  integer ipx

  integer ng,ig,ip,ipmax,icount
  
  integer i,ii
  real*8 dist,disti,distance

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

function distance(p)
  implicit none
  real*8 distance,p(2),pp(2),pi,twopi
  integer i
  pi=1; twopi=2   !pi set as unity here
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
  mesh%g(1)%np=np; mesh%g(1)%w=wmax; mesh%g(1)%dp=dp

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

     mesh%g(ig+1)%np=icount; mesh%g(ig+1)%dp=dp; mesh%g(ig+1)%w=w
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

    real*8 pc,dp
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
      
    icount=0
	do ig=1,ng; do ip=1,line%g(ig)%np; if(line%g(ig)%info(ip)/=0)cycle
	   icount=icount+1
	end do; end do
    line%npoints=icount
    return
  end subroutine linemesh
