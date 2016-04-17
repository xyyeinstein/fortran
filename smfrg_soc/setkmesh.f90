subroutine meshkvec()
  use workplace, only : norb,nk,kv,kmesh,nkx,nky,ek,model,appendkmesh
  implicit none

  integer ng,ig,ip
  real*8 wmax,wmin

  integer ix,iy
  real*8 vk(2)
  complex*16 hk(norb,norb)
  real*8 gapmesh
  integer ik,iorb

  allocate(ek(norb,-nkx:nkx,-nky:nky))
  !allocate(uk(norb,norb,-nkx:nkx,-nky:nky))

  do iy=-nky,nky; do ix=-nkx,nkx; vk=(/ix,iy/)*1./(/nkx,nky/)
     call gethk(norb,vk,hk,model)
     call ZHEIGEN(norb,hk,ek(:,ix,iy))
  end do; end do
  print*,'eigentable done.'
  call outputdos()
  print*,'dos done.'

  if(appendkmesh)then
    call openkmesh(); ng=kmesh%ng
  else
    ng=8 !open(31,file='ng.input'); read(31,*)ng; close(31)
    wmax=15; wmin=1.e-5
    call BZmesh(ng,wmax,wmin,kmesh,ekmin)
    call savekmesh()
  end if

  nk=kmesh%npoints; print*,'# kmesh points =',nk
  allocate(kv(3,nk))

  nk=0
  do ig=1,ng; do ip=1,kmesh%g(ig)%np
     if(kmesh%g(ig)%info(ip)/=0)cycle
	 nk=nk+1
     kv(1:2,nk)=kmesh%g(ig)%p(:,ip)
	 kv(3,nk)=kmesh%g(ig)%dp**2
  end do; end do

  if(sum(kv(3,:))/=1)stop 'jacobi error @ meshkvec'

  open(10,file='kv.dat')
  write(10,100)kv; close(10)
100 format(1x,3f20.10)

  open(10,file='fs.dat')
  do ip=1,kmesh%g(ng)%np; write(10,100)kmesh%g(ng)%p(:,ip); end do
  close(10)

  gapmesh=1.e10
  do ik=1,nk; vk=kv(1:2,ik)
     call gethk(norb,vk,hk,model)
	 call ZHEIGEN(norb,hk,ek(:,1,1))
	 do iorb=1,norb; gapmesh=min(gapmesh,abs(ek(iorb,1,1))); end do
  end do
  print*,'mesh gap =',gapmesh   

  return
  contains
    function ekmin(vk)
	  real*8 ekmin,vk(2)
	  real*8 eval(norb),searchtable
	  complex*16 hk(norb,norb)
	  integer iorb
	  
	  call gethk(norb,vk,hk,model)
	  call ZHEIGEN(norb,hk,eval)
      !do iorb=1,norb; eval(iorb)=searchtable(vk,nkx,nky,ek(iorb,:,:)); end do

	  ekmin=1.d10
	  do iorb=1,norb
	     ekmin=min( ekmin, abs(eval(iorb)) )
	  end do
	  
	  return
   end function ekmin
end subroutine meshkvec



 subroutine FlexibleMesh(La)
    use workplace, only : nk,kv,kmesh
    implicit none
    real*8 La

    integer ig,ip,ngmin

    ngmin=0
    do ig=1,kmesh%ng; if(kmesh%g(ig)%w>La/100)then; ngmin=ngmin+1; end if; end do
    ngmin=max(ngmin,2)

    nk=0
    do ig=1,ngmin; do ip=1,kmesh%g(ig)%np
       if(kmesh%g(ig)%info(ip)/=0.and.ig/=ngmin)cycle
	   nk=nk+1; kv(1:2,nk)=Kmesh%g(ig)%p(:,ip); kv(3,nk)=Kmesh%g(ig)%dp**2
    end do; end do

    if(sum(kv(3,1:nk))/=1)stop 'jacobi error'

    return
 end subroutine FlexibleMesh


  subroutine Savekmesh()
    use workplace, only : kmesh,kmeshfile
    integer ig,np,ng
    open(17,file=kmeshfile)
	ng=kmesh%ng
	write(17,*)ng
	do ig=1,ng
	   np=kmesh%g(ig)%np
	   write(17,*)np
	   write(17,*)kmesh%g(ig)%w
	   if(np==0)cycle
	   write(17,*)kmesh%g(ig)%dp
	   write(17,*)kmesh%g(ig)%info
	   write(17,*)kmesh%g(ig)%p
	end do
    write(17,*)kmesh%npoints
	close(17)
    return
  end subroutine SaveKmesh
   
  subroutine OpenKmesh()
    use workplace, only : kmesh,kmeshfile
    integer ig,ng,np
	open(17,file=kmeshfile)
	read(17,*)ng; kmesh%ng=ng
	allocate(kmesh%g(ng))
	do ig=1,ng
	   read(17,*)np; kmesh%g(ig)%np=np
	   read(17,*)kmesh%g(ig)%w
	   read(17,*)kmesh%g(ig)%dp
	   if(np==0)cycle
	   allocate(kmesh%g(ig)%info(np))
	   allocate(kmesh%g(ig)%p(2,np))
	   read(17,*)kmesh%g(ig)%info
	   read(17,*)kmesh%g(ig)%p
	end do
	read(17,*)kmesh%npoints
	close(17)
    return
  end subroutine OpenKmesh


subroutine symmetrizeKv()
  use workplace
  implicit none

  integer ik,ikg,ig,ngmax,icount
  logical reducible
  real*8 pv(3),pg(3),weight
  real*8, allocatable, dimension (:,:) :: kvv
  
  allocate(kvv(3,2*nk))

  icount=0; pv=0
  do ik=1,nk; pv(1:2)=kv(1:2,ik); if(reducible(pv(1:2),group))cycle
     ngmax=8; if(abs(pv(2))<1.e-5.or.abs(pv(1)-pv(2))<1.e-5)ngmax=4
	 weight=kv(3,ik)
	 do ig=1,ngmax; call groupaction(ig,pv,pg,'C4v')
	    icount=icount+1
		kvv(1:2,icount)=pg(1:2)
		kvv(3,icount)=weight
     end do
  end do

  nk=icount
  if(abs(sum(kvv(3,1:nk))-1)>1.e-5)pause 'jaccobi error @ symmetrizekv'

  if(nk>size(kv(1,:)))then
    deallocate(kv)
	allocate(kv(3,nk))
  end if
  kv(:,1:nk)=kvv(:,1:nk)
  deallocate(kvv)
  return
end subroutine symmetrizekv
