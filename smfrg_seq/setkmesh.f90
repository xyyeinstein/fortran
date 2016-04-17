subroutine meshkvec()
  use workplace, only : norb,nk,kv,kmesh,model,appendkmesh,usegroup,wir,useX0,uniformkmesh,MF
  implicit none

  integer ng,ig,ip
  real*8 wmax,wmin

  integer ik,iorb
  real*8 vk(2),area,gapmesh,Q,eval(norb)
  complex*16 hk(norb,norb)
  logical hexagonal

  call eigentable()

  if(model%ladder)then; call meshkvec_ladder(); return; end if

  if(uniformkmesh)then; call meshkvec_uniform(); return; end if
  Q=0
  hexagonal=.not.( sum(model%a*model%b)==0 )

  

  if(appendkmesh)then
    call openkmesh(); ng=kmesh%ng
  else
    ng=8; wmax=4; wmin=1.e-5
    !ng=10; wmax=6; wmin=1.e-3
    if(hexagonal)then
	  call HexBZmesh(ng,wmax,wmin,kmesh,ekmin)
    else
	  call BZmesh(ng,wmax,wmin,kmesh,ekmin)
    end if
	call savekmesh()
  end if

  nk=kmesh%npoints; print*,'number of k-mesh points=',nk
  allocate(kv(3,nk))

  nk=0
  do ig=1,ng; do ip=1,kmesh%g(ig)%np
     if(kmesh%g(ig)%info(ip)/=0)cycle
	 nk=nk+1
     kv(1:2,nk)=kmesh%g(ig)%p(:,ip)
	 kv(3,nk)=kmesh%g(ig)%area
  end do; end do

  if(abs( sum(kv(3,1:nk))-1 )>1.e-5)then
    print*,sum(kv(3,1:nk))
	stop 'jaccobi error @ meshkvec'
  end if

  if(hexagonal)then
    call symmetrizeKv()
    print*,'number of k-mesh points after symmetrization =',nk
  end if

  open(10,file='kv.dat')
  write(10,100)kv(:,1:nk); close(10)
100 format(1x,3f20.10)

  open(10,file='fs_.dat')
  do ip=1,kmesh%g(ng)%np; write(10,100)kmesh%g(ng)%p(:,ip); end do
  !do ip=1,kmesh%g(ng-1)%np; write(10,100)kmesh%g(ng-1)%p(:,ip); end do
  close(10)

  gapmesh=1.e10
  do ik=1,nk; vk=kv(1:2,ik)
     call gethk(norb,vk,Q,hk,model,MF,.false.)
	 call ZHEIGEN(norb,hk,eval)
	 do iorb=1,norb; gapmesh=min(gapmesh,abs(eval(iorb))); end do
  end do

  print*,'meshgap =',gapmesh
  if(gapmesh>wir/5)then
    print*,'meshgap too large versus the infrared limit wir:'
	print*,'wir =',wir
	if(useX0)print*,'X0 is to be used, but may be unreliable because of the bad kmesh.'
  end if

  return
  contains
    function ekmin(vk)
	  real*8 ekmin,vk(2)
	  real*8 eval(norb),searchtable
	  complex*16 hk(norb,norb)
	  integer iorb
	  
	  call gethk(norb,vk,Q,hk,model,MF,.false.)
	  call ZHEIGEN(norb,hk,eval)

	  ekmin=1.d10
	  do iorb=1,norb
	     ekmin=min( ekmin, abs(eval(iorb)) )
	  end do
	  
	  return
   end function ekmin
end subroutine meshkvec

subroutine meshkvec_uniform()
  use workplace
  implicit none

  integer ikx,iky,ik
  real*8 kx,ky

  nk=nkx*nky
  allocate(kv(3,nk))
  
  ik=0
  do ikx=1,nkx; do iky=1,nky
    kx=-1.d0+2.d0*(ikx-1)/nkx+1.d0/nkx
	ky=-1.d0+2.d0*(iky-1)/nkx+1.d0/nkx
	ik=ik+1
    kv(1:2,ik)=(/kx,ky/)
	kv(3,ik)=1.d0/nk
  end do; end do

  return
end subroutine

subroutine FlexibleMesh(La)
    use workplace, only : nk,kv,kmesh,model
    implicit none
    real*8 La

    integer ig,ip,ngmin

    ngmin=0
    do ig=1,kmesh%ng; if(kmesh%g(ig)%w>La/100)then; ngmin=ngmin+1; end if; end do
    ngmin=max(ngmin,2)

    nk=0
    do ig=1,ngmin; do ip=1,kmesh%g(ig)%np
       if(kmesh%g(ig)%info(ip)/=0.and.ig/=ngmin)cycle
	   nk=nk+1; kv(1:2,nk)=Kmesh%g(ig)%p(:,ip)
	   kv(3,nk)=kmesh%g(ig)%area
    end do; end do

    if(abs( sum(kv(3,1:nk))-1 )>1.e-5)stop 'jaccobi error @ flexiblemesh'
    if( sum(model%a*model%b) /=0 )call symmetrizeKv()
    return
end subroutine FlexibleMesh

subroutine symmetrizeKv()
  use workplace
  implicit none

  integer ik,ikg,ig,ngmax,icount
  logical reducible
  real*8 pv(3),pg(3),weight
  real*8, allocatable, dimension (:,:) :: kvv
  
  allocate(kvv(3,2*nk))

  icount=0; pv=0
  do ik=1,nk; pv(1:2)=kv(1:2,ik); if(reducible(pv(1:2),'C6v'))cycle
     ngmax=12; if(abs(pv(1)-pv(2)*sqrt(3.))<1.e-4)ngmax=6
	 weight=kv(3,ik)
	 do ig=1,ngmax; call groupaction(ig,pv,pg,'C6v')
	    icount=icount+1; if(icount>2*nk)pause 'icount exceeds limit'
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

subroutine Savekmesh()
    use workplace, only : kmesh,kmeshfile,model
    integer ig,np,ng

    if(model%ladder)return

    open(17,file=kmeshfile)
	ng=kmesh%ng
	write(17,*)ng
	do ig=1,ng
	   np=kmesh%g(ig)%np
	   write(17,*)np
	   write(17,*)kmesh%g(ig)%w
	   write(17,*)kmesh%g(ig)%area
	   if(np==0)cycle
	   write(17,*)kmesh%g(ig)%info
	   write(17,*)kmesh%g(ig)%p
	end do
    write(17,*)kmesh%npoints
	close(17)
    return
end subroutine SaveKmesh
   
subroutine OpenKmesh()
    use workplace, only : kmesh,kmeshfile,model
    integer ig,ng,np

    if(model%ladder)then; call meshkvec_ladder(); return; end if    

	open(17,file=kmeshfile)
	read(17,*)ng; kmesh%ng=ng
	allocate(kmesh%g(ng))
	do ig=1,ng
	   read(17,*)np; kmesh%g(ig)%np=np
	   read(17,*)kmesh%g(ig)%w
	   read(17,*)kmesh%g(ig)%area
	   !kmesh%g(ig)%area=kmesh%g(ig)%area**2 !!!!!!!!!!!!!!!!!!!!!!!!!!
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


subroutine meshkvec_ladder()
  use standard_derived_types
  use workplace, only : norb,nk,kv,model,wir,useX0,MF
  implicit none

  integer ng,ig,ip
  real*8 wmax,wmin
  type (meshconfig), dimension (2*model%nleg) :: line

  integer ik,iky,iorb,nky,icount
  real*8 ky,kL(2),kU(2),vk(2),dk,gapmesh,Q,eval(norb)
  complex*16 hk(norb,norb)

  Q=0; nky=model%nleg; if(nky==0)stop 'invalid # of legs in the ladder model @ meshkvec_ladder'

  ng=16; wmax=6; wmin=1.e-5; dk=2./nky
  
  icount=0
  do iky=0,nky-1; ky=iky*dk; if(ky>1)ky=ky-2
     kL=(/0.d0,ky/); kU=(/1.d0,ky/); icount=icount+1; call linemesh(ng,wmax,wmin,kL,kU,line(icount),epfunc)
     kL=(/-1.d0,ky/); kU=(/0.d0,ky/); icount=icount+1; call linemesh(ng,wmax,wmin,kL,kU,line(icount),epfunc)
  end do; if(icount/=2*nky)stop '# of linemeshes /= 2*nleg @ meshkvec_ladder'

  nk=sum(line(:)%npoints);  allocate(kv(3,nk))

  open(10,file='fs.dat')
  ik=0; icount=0
  do iky=0,nky-1; ky=iky*dk; if(ky>1)ky=ky-2
     kL=(/0.d0,ky/); kU=(/1.d0,ky/); icount=icount+1
     do ig=1,ng; do ip=1,line(icount)%g(ig)%np
       if(line(icount)%g(ig)%info(ip)/=0)cycle
	   ik=ik+1; kv(1:2,ik)=line(icount)%g(ig)%p(1,ip)*(kU-kL)+kL
	   kv(3,ik)=line(icount)%g(ig)%dp
	   if(ig==ng)write(10,100)kv(:,ik)
     end do; end do

     kL=(/-1.d0,ky/); kU=(/0.d0,ky/); icount=icount+1
     do ig=1,ng; do ip=1,line(icount)%g(ig)%np
       if(line(icount)%g(ig)%info(ip)/=0)cycle
	   ik=ik+1; kv(1:2,ik)=line(icount)%g(ig)%p(1,ip)*(kU-kL)+kL
	   kv(3,ik)=line(icount)%g(ig)%dp
	   if(ig==ng)write(10,100)kv(:,ik)
    end do; end do
  end do; close(10)
  kv(3,:)=kv(3,:)/sum(kv(3,:))
  
  open(10,file='kv.dat')
  write(10,100)kv(:,1:nk); close(10)
100 format(1x,3f20.10)

  gapmesh=1.e10
  do ik=1,nk; vk=kv(1:2,ik)
     call gethk(norb,vk,Q,hk,model,MF,.false.)
	 call ZHEIGEN(norb,hk,eval)
	 do iorb=1,norb; gapmesh=min(gapmesh,abs(eval(iorb))); end do
  end do
  if(gapmesh>wir/5)then
    print*,'meshgap too large versus the infrared limit wir:'
	print*,'meshgap =',gapmesh
	print*,'wir =',wir
	if(useX0)print*,'X0 is to be used, but may be unreliable because of the bad kmesh.'
  end if

  return
  contains
      function epfunc(p)
        implicit none
        real*8 epfunc,p,Q,emin
		real*8, dimension (2) :: kv
		complex*16 hk(norb,norb)
		real*8 ek(norb)
		integer ib
		kv=kL+(kU-kL)*p; emin=1.e10; Q=0
		call gethk(norb,kv,Q,hk,model,MF,.false.); call ZHEIGEN(norb,hk,ek)
        do ib=1,norb; emin=min(emin,abs(ek(ib))); end do
        epfunc=emin
        return
      end function epfunc
end subroutine meshkvec_ladder

subroutine pnest()
  use workplace
  
  integer iorb,jorb
  real*8 k(2),eta,nest,evk(norb),evkq(norb),dos
  complex*16 hk(norb,norb),hkq(norb,norb)
  integer stat1

  eta=1.e-0
  
  open(10,file='kv.dat')
    do while(.true.)
        read(10,*,iostat=stat1)k,dos
        call gethk(norb,k,0.d0,hk,model,mf,.false.)
        call ZHEIGEN(norb,hk,evk)
        k=k+(/1,1/)
        call gethk(norb,k,0.d0,hkq,model,mf,.false.)
        call ZHEIGEN(norb,hkq,evkq)
        do iorb=1,norb
            do jorb=1,norb
                nest=nest+1/(1+(evk(iorb)-evkq(jorb))**2/eta**2)
            end do
        end do
		if(stat1/=0)exit
    end do
    close(10)
    
    nest=nest/norb/norb/nk
    print*,nest

  return
end subroutine

