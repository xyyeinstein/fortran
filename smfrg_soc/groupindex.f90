subroutine setgroupindex()
  use workplace
  implicit none

  integer idim,ig,iq,iqg
  real*8 pv(3),pg(3),Xg
  character*3 formtable(2)
  integer ip(2),jdim
  integer partner
  integer info(nq)
  logical active,reducible
  real*8 formvalue

  type (Mformconfig), pointer :: Mf,M1,M2
  integer spin,nbspin
  character*2 Sform
  complex*16 Xpp,Xph,SpinXg(2)
  integer SpinIndexg(2)
  character*1 chr(2)
  integer iidim,jjdim,findMform
  character*3 formg

  chr=(/'+','-'/); group=model%group

  allocate(indexgroup(ng,ndim))
  allocate(ppXg(ng,ndim),phXg(ng,ndim))

  do idim=1,ndim; do ig=1,ng
     Mf=>Mformfunc(idim)
	 call spinimage(ig,group,SpinIndexg,SpinXg)
     spin=Mf%spin; nbspin=Mf%nbspin
	 Sform=chr(SpinIndexg(spin))//chr(SpinIndexg(nbspin))
	 Xpp=conjg(SpinXg(spin)*SpinXg(nbspin))
	 Xph=conjg(SpinXg(spin))*SpinXg(nbspin)
	 Xpp=conjg(Xpp); Xph=conjg(Xph)   !??

	 call orbitimage(ig,group,Mf%form,formg,Xg)
	 iidim=findMform(Mf%atom,Mf%nbatom,Mf%Lformfunc%r(:,1),formg,Sform)

     indexgroup(ig,idim)=iidim
	 ppXg(ig,idim)=Xpp*Xg
     phXg(ig,idim)=Xph*Xg
  end do; end do

  allocate(indexqg(ng,nq))
  pv=0; info=-1
  do iq=1,nq; pv(1:2)=qv(1:2,iq); active=.not.reducible(pv(1:2),group)
     do ig=1,ng; call groupaction(ig,pv,pg,group)
        call search(pg(1:2),iqg,qmesh); indexqg(ig,iq)=iqg
		if(active)info(iqg)=iqg
	    if(sum( abs(pg(1:2)-qv(1:2,iqg) ) )>1.e-5)stop 'qv not symmetric @ setgroupindex'
     end do
  end do

  do iq=1,nq
     if(info(iq)<0)then
	   print*,iq,info(iq)
	   print*,qv(1:2,iq)
	   pause 'this qv is not covered by irreducible zone + group @ setgroupindex'
     end if
  end do
  
  return
end subroutine setgroupindex


 function findMform(atom,nbatom,rv,form,Sform)
   use workplace
   implicit none
   integer atom,nbatom,findMform
   real*8 rv(3)
   character*3 form
   character*2 Sform
   type (Mformconfig), pointer :: Mi,Mf

   real*8 rnorm
   integer icount,i

   rnorm=sum(rv*rv)
   icount=0
   do i=1,ndim
      Mf=>Mformfunc(i)
	  if(Mf%atom/=atom)cycle
	  if(Mf%nbatom/=nbatom)cycle
	  if(Mf%form/=form)cycle
      if(Mf%Sform/=Sform)cycle
      if( abs( sum(Mf%Lformfunc%r(:,1)**2)-rnorm )>1.e-4 )cycle
      icount=icount+1
	  findMform=i
   end do

   if(icount==0)stop 'Mform not found @ findMform'
   if(icount>1)stop 'Too many Mform found @ findMform'
   if(findMform<1.or.findMform>ndim)stop 'Mform found but invalid @ findMform'

   return
 end function findMform
   
