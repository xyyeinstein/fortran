
subroutine setgroupindex()
  use workplace
  implicit none

  integer idim,ig,iq,iqg,ik,ikg
  real*8 pv(3),pg(3),qq(2),k(2),kg(2)
  integer info(nq)
  logical active,reducible,hexagonal
  integer findq

  hexagonal=( abs(sum(model%a*model%b))>1.e-5)

  allocate(indexgroup(8,ng,ndim),Xg(8,ng,ndim))

  do idim=1,ndim; do ig=1,ng
	 call Mformimage(ig,idim,indexgroup(:,ig,idim),Xg(:,ig,idim))
	 !print*,idim,ig
	 !print*,indexgroup(:,ig,idim)
	 !pause
  end do; end do


  allocate(indexqg(ng,nq))
  allocate(indexkg(ng,nk))
  indexqg=0
  indexkg=0

  pv=0; info=-1
  do iq=1,nq
    if(abs(qv(3,iq))<1.e-10)then
	  info(iq)=iq
	  indexqg(:,iq)=iq
	  cycle
	end if
	pv(1:2)=qv(1:2,iq)
    active=.not.reducible(pv(1:2),group)
    if(active)info(iq)=iq
    do ig=1,ng
	  call groupaction(ig,pv,pg,group)
	  iqg=findq(pg(1:2),hexagonal)
	  call rbzkv(pg(1:2)-qv(1:2,iqg),model%a,model%b,qq)
	  if(sum(abs(qq))>1.e-5)then
	    print*,pg(1:2),qv(1:2,iqg)
        print*,qq(1:2)
		print*, 'an image of qv does not exist and is redirected to qv itself @ setgroupindex.'
		read*
		indexqg(ig,iq)=iq; cycle
	  end if
	  if(active)info(iqg)=iqg
	  indexqg(ig,iq)=iqg
    end do
  end do

  do iq=1,nq; if(info(iq)<=0)pause 'qv not covered by group @ setgroupindex'; end do
  
  return
  do ik=1,nk
    k=kv(1:2,ik)
	if(reducible(k,group))cycle
	print*,ik
	do ig=1,ng
	  call groupaction(ig,k,kg,group)
	  call searchkv(ikg,kg)
	  indexkg(ig,ik)=ikg
	end do
  end do

  return
end subroutine setgroupindex

function findq(q,hexagonal)
  use workplace
  implicit none
  integer findq
  real*8 q(2)
  logical hexagonal

  integer iq
  real*8 qq(2),dist,error

  if(model%ladder)then
    error=1.e10
    do iq=1,nq
      call rbzkv(q-qv(1:2,iq),model%a,model%b,qq)
	  dist=sum( abs(qq) )
	  if(dist<error)then
	    error=dist
		findq=iq
	  end if
	end do
  else 
	if(hexagonal)then
      call hexsearch(q,findq,qmesh)
    else
      call search(q,findq,qmesh,uniformqmesh)
    end if
  end if
  return
end function findq


subroutine MformImage(ig,iMform,indexMg,XMg)
  use workplace
  implicit none

  integer ig,iMform,indexMg(8)
  real*8 XMg(8)

  type (Mformconfig), pointer :: M
  type (Oformconfig), pointer :: O
  type (Lformconfig), pointer :: L
  integer iO,iL,jM,atom,nbatom,icount,morb
  integer atomimage
  integer indexOg(4),indexLg(2)
  real*8 XOg(4),XLg(2),value

  morb=norb/natom; if(model%sample=='CuO2')morb=norb

  M=>Mformfunc(iMform)
  call Oformimage(ig,M%Oform,nOform,Oformfunc,morb,model%orbits(1:morb),group,indexOg,XOg)
  call Lformimage(ig,M%Lform,nLform,Lformfunc,group,indexLg,XLg)
  atom=atomimage(ig,M%atom,model);  nbatom=atomimage(ig,M%nbatom,model)	 

  icount=0; indexMg=0; XMg=0
  do jM=1,ndim; M=>Mformfunc(jM)
     if(M%atom/=atom.or.M%nbatom/=nbatom)cycle
	 do iO=1,4; if(indexOg(iO)<=0)cycle; if(M%Oform/=indexOg(iO))cycle
	 do iL=1,2; if(indexLg(iL)<=0)cycle; if(M%Lform/=indexLg(iL))cycle
	    value=XOg(iO)*XLg(iL); if(abs(value)<1.e-6)cycle
		icount=icount+1; if(icount>8)stop 'too many Mform images @ Mformimage'
		indexMg(icount)=jM; XMg(icount)=value
	 end do; end do
  end do

  if(icount==0)stop 'MformImages not found @ MformImage'

  return
end subroutine MformImage





