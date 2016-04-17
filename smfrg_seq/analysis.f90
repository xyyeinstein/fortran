subroutine FermiSurf(norb,Q,model)
  use standard_derived_types
  implicit none
  integer norb
  real*8 Q
  type (modelblock) :: model
	type (mfblock) :: mf

  integer ix,iy,iorb,iband,nl
  real*8 pi,kv(2),Ak_orb
  complex*16 hk(norb,norb),one
  real*8 ek(norb),eta,dos,scale
  integer stat1

  pi=asin(1.d0)*2; one=cmplx(0,1); eta=0.015
  scale=1; if( abs( sum(model%a*model%b) )>1.e-5 )scale=1.4
 
  open(10,file='Ak.dat'); open(11,file='Ak_orb.dat')
	nl=40
  do iy=-nl,nl;  do ix=-nl,nl; kv=(/ix,iy/)*scale/nl
    call gethk(norb,kv,Q,hk,model,mf,.false.); call ZHEIGEN(norb,hk,ek)
  do iband=1,norb; dos = - aimag( 1/ (one*eta-ek(iband)) )/pi; write(10,*)dos; end do
  do iorb=1,norb; 
    Ak_orb=0;
    do iband=1,norb
      dos=- aimag( 1/ (one*eta-ek(iband)) )/pi
      Ak_orb=Ak_orb+dos*hk(iorb,iband)*conjg(hk(iorb,iband))
    end do
    write(11,*)Ak_orb
  end do
  end do; end do
  close(10); close(11)
  print*,'band-wise Ak at fermi energy output to file'
  open(11,file='fs.dat')
  open(10,file='ak_fs.dat')
  do while(.true.)
    read(11,100,iostat=stat1)kv,iband,dos
    call gethk(norb,kv,Q,hk,model,mf,.false.);
    call ZHEIGEN(norb,hk,ek)
    write(10,'(2e16.8)')hk(:,iband)
    if(stat1/=0)exit
  end do
  close(10)
  close(11)
100 format(1x,2f20.6,i6,f20.6)
  return
end subroutine FermiSurf


subroutine analysis(La,isave,iXq)
  use workplace
  implicit none
  real*8 La
  logical isave,iXq

  integer i,info(3),iq(3)
  real*8 q(2),S(3),Vmax,Smax
  real*8 pv(3),pg(3)
  integer ig,iqg
  real*8 findmax

  complex*16, dimension (ndim,ndim,nq) :: Peff,Ceff
  logical reducible

  Peff=0; Ceff=0
  if(model%ephcoupled)call ElectronPhonon(La,Peff,Ceff)
  Peff=Peff+P; Ceff=Ceff+C

  Vmax=findmax(ndim,nq,Peff,info)
  Vmax=max(Vmax,findmax(ndim,nq,-Ceff,info)); iq(2)=abs(info(3))
  Vmax=max(Vmax,findmax(ndim,nq,2*D-Ceff,info)); iq(3)=abs(info(3))

  Smax=-1.e10
  !do i=1,nq; q=qv(1:2,i); if(usegroup.and.reducible(q,group))cycle
  !   call channel('P',q,ndim,Peff(:,:,i),S(1),.false.,.false.)
  !	 if(S(1)>Smax)then; Smax=S(1); iq(1)=i; end if
  !end do
  iq(1)=nq-model%nest; !print*,qv(1:2,iq(1)); pause 'ok?'
  !call channel('P',qv(1:2,iq(1)),ndim,Peff(:,:,iq(1)),S(1),isave,.false.)
  call pairchannel(ndim,S(1),isave)   

  Smax=-1.e10; iq(2)=0
  do i=1,nq; q=qv(1:2,i); if(usegroup.and.reducible(q,group))cycle
     call channel('C',q,ndim,-Ceff(:,:,i),S(2),.false.,.false.)
     if(S(2)>Smax)then; Smax=S(2); iq(2)=i; end if
  end do
  call channel('C',qv(1:2,iq(2)),ndim,-Ceff(:,:,iq(2)),S(2),isave,iXq)

  Smax=-1.e10; iq(3)=0
  do i=1,nq; q=qv(1:2,i); if(usegroup.and.reducible(q,group))cycle
     call channel('D',q,ndim,2*D(:,:,i)-Ceff(:,:,i),S(3),.false.,.false.)
     if(S(3)>Smax)then; Smax=S(3); iq(3)=i; end if
  end do
  call channel('D',qv(1:2,iq(3)),ndim,2*D(:,:,iq(3))-Ceff(:,:,iq(3)),S(3),isave,iXq)

  q=(1.d0,1.d0)
  call search(q,iq(3),qmesh,.false.)
  call channel('P',qv(1:2,iq(3)),ndim,P(:,:,iq(3)),S(3),isave,iXq)

  output%Vmax(1)=Vmax; output%Vmax(2:4)=S
  do i=1,3; output%q(:,i)=qv(1:2,iq(i)); end do

  return
end subroutine analysis


subroutine channel(name,q,nA,A,Smax,isave,iXq)
    use workplace, only : Mformconfig,Lformconfig,Oformconfig,Mformfunc,Lformfunc,Oformfunc,model,output
    implicit none
	character*1 name
	real*8 q(2)
	integer nA,iS,iSS,info(nA)
	real*8 Smax,Si(nA)
	complex*16, dimension (nA,nA) :: A
	logical isave,iXq

	complex*16, dimension (nA,nA) :: B,U,V
	real*8 Sig(nA)
	integer idim,jdim,iatom,iorb,iform,mode(nA),npairmode

	character*3 Mform,Lform,Oform,form(nA)
	type (Mformconfig), pointer :: Mi
	type (Lformconfig), pointer :: Li
	type (Oformconfig), pointer :: Oi
	real*8, pointer, dimension (:,:) :: ra
	real*8 pi,Umax
	complex*16 phase,one

  B=A; call ZSVD( nA,nA,B,U,Sig,V )
	do iS=1,nA
	   Si(iS)=Sig(iS)*real( sum( U(:,iS)*V(iS,:) ) )
	end do
	Si=-Si;	call ascendingorder(nA,Si,info); Si=-Si
	Smax=Si(1)

	if(name=='PP')then
	  npairmode=min(10,nA)
	  do iS=1,npairmode; iSS=info(iS)
	    Umax=abs(U(1,iSS)); mode(iS)=1; form(iS)=mformfunc(1)%form
	    do idim=2,nA
	      if(abs(U(idim,iSS))>Umax)then
		      Umax=abs(U(idim,iSS)); mode(iS)=idim; form(iS)=mformfunc(idim)%form
		    end if
	    end do
	  end do
    do is=1,10
  	  output%LeadingPair(1:npairmode)%S=Si(1:npairmode)
	    output%LeadingPair(1:npairmode)%iS=mode(1:npairmode)
	    output%LeadingPair(1:npairmode)%form=form(1:npairmode)
    end do
  end if
    if(.not.isave)return

    !gauge transform to proper gauge
	pi=asin(1.d0); one=cmplx(0,1)
	do idim=1,nA
	   Mi=>Mformfunc(idim); ra=>model%ra; phase=exp(one*pi*sum(q*ra(1:2,Mi%atom)))
	   U(idim,:)=U(idim,:)*phase; V(:,idim)=V(:,idim)*conjg(phase)
	end do

	if(name=='P')then;        open(10,file='pairMode.txt'); open(11,file='pairmode.dat')
	  else if(name=='C')then; open(10,file='sdwMode.txt')
	  else if(name=='D')then; open(10,file='cdwMode.txt')
	end if

    write(10,*)'  ordering momentum/pi='; write(10,*)q
	write(10,*)' '
	do iS=1,nA; if(abs(Sig(iS))<1.e-10)cycle; iSS=info(iS)
	   write(10,*)'   iS, S =',iS,Si(iS)
       write(10,*)'   Mform    Lform     rbasis       Oform           Uvec               Vvec'
	   do idim=1,nA; if(abs(U(idim,iSS))+abs(V(idim,iSS))<1.e-4)cycle
	      Mi=>Mformfunc(idim); Mform=Mi%form
		  Li=>Lformfunc(Mi%Lform); Oi=>Oformfunc(Mi%Oform)
	      Lform=Li%form; Oform=Oi%form
	      if(Li%ndim==2)then
		    write(10,100)idim,Mi%atom,Mform,Lform,Li%r(:,1),Oform,Oi%parity,Oi%pair(:,1),U(idim,iSS),V(iSS,idim)
	      else if(Li%ndim==3)then
		    write(10,200)idim,Mi%atom,Mform,Lform,Li%r(:,1),Oform,Oi%parity,Oi%pair(:,1),U(idim,iSS),V(iSS,idim)
	      end if
	   end do; write(10,*)' ' 
	end do; close(10)
100 format(1x, i4, i3, ' ', a3, ' L-', a3, 2f5.1, ' O-', a3,a1, ' ', 2i2, 4f8.4)
200 format(1x, i4, i3, ' ', a3, ' L-', a3, 3f5.1, ' O-', a3,a1, ' ', 2i2, 4f8.4)

    if(name=='P')then
       open(11,file='pairmode.dat')
       do iS=1,nA; iSS=info(iS)
	      write(11,*)iS,Si(iS)
	      write(11,300)U(:,iSS)
	   end do
	   call GapForm_Band(U(:,info(1)),1)
       close(11)
	   if(iXq)call bubbleSusceptibility()
	else 
	   if(name=='C')then; open(11,file='sdwmode.dat')
	   else; open(11,file='cdwmode.dat'); end if
       write(11,*)q
       do iS=1,nA; iSS=info(iS)
          write(11,*)iS,Si(iS)
	      write(11,300)U(:,iSS)
	   end do
	   close(11)
	   if(name=='C')call RealSpaceTexture('sdw',U(:,info(1)),q,1)
	   if(name=='D')call RealSpaceTexture('cdw',U(:,info(1)),q,1)      	    
	end if
300 format(1x,2f20.6)
    return
end subroutine channel


subroutine pairchannel(nA,Smax,isave)
  use workplace, only : P,nq,qv,Mformconfig,Lformconfig,Oformconfig,Mformfunc,Lformfunc,Oformfunc,model,output
  implicit none

  integer nA
  real*8 Smax
  logical isave
  complex*16 A(nA,nA)

  integer imode,is,iss,iq,loc(2)
  real*8 q(2),smin,si(nA),sig(na,nq)
	complex*16, dimension (nA,nA) :: U,V
	integer idim,jdim,iatom,iorb,iform

	character*3 Mform,Lform,Oform
	type (Mformconfig), pointer :: Mi
	type (Lformconfig), pointer :: Li
	type (Oformconfig), pointer :: Oi
	real*8, pointer, dimension (:,:) :: ra
	real*8 pi,explain
	complex*16 phase,one

  do iq=1,nq
    A=P(:,:,iq)
    call ZSVD( nA,nA,A,U,Sig(:,iq),V )
    do iS=1,nA
	    Sig(iS,iq)=Sig(iS,iq)*real( sum( U(:,iS)*V(iS,:) ) )
	  end do
	end do
  smax=maxval(sig)
  smin=minval(Sig)
  do imode=1,10
    output%LeadingPair(imode)%S=maxval(Sig)
    loc=maxloc(sig)
    output%LeadingPair(imode)%is=loc(1)
    output%leadingpair(imode)%iq=loc(2)
    output%leadingpair(imode)%q=qv(1:2,loc(2))
    Sig(loc(1),loc(2))=smin
  end do
	
  if(.not.isave)return
  
  open(10,file='leadingpairMode.txt'); open(11,file='leadingpairmode.dat')
    
  do imode=1,10
    q=output%leadingpair(imode)%q
    is=output%leadingpair(imode)%is
    iq=output%leadingpair(imode)%iq
    A=P(:,:,iq)
    call ZSVD(nA,nA,A,U,Si,V)
    do iss=1,nA
      si(iss)=si(iss)*real(sum(U(:,iss)*V(iss,:)))
    end do
    !gauge transform to proper gauge. 
    !this is necessary for nonzero-q Cooper pairing, 
    !but is an identity operation for q=0
	  pi=asin(1.d0); one=cmplx(0,1)
	  do idim=1,nA
	    Mi=>Mformfunc(idim); ra=>model%ra; phase=exp(one*pi*sum(q*ra(1:2,Mi%atom)))
	    U(idim,:)=U(idim,:)*phase; V(:,idim)=V(:,idim)*conjg(phase)
	  end do
    write(10,*)'   iS, S =',imode,si(is)
    write(10,*)'  ordering momentum/pi='
    write(10,*)q
  	write(10,*)' '
	  write(10,*)'   Mform    Lform     rbasis       Oform           Uvec               Vvec'
	  do idim=1,nA; if(abs(U(idim,is))+abs(V(idim,is))<1.e-4)cycle
	    Mi=>Mformfunc(idim); Mform=Mi%form
		  Li=>Lformfunc(Mi%Lform); Oi=>Oformfunc(Mi%Oform)
	    Lform=Li%form; Oform=Oi%form
	    if(Li%ndim==2)then
		    write(10,100)idim,Mi%atom,Mform,Lform,Li%r(:,1),Oform,Oi%parity,Oi%pair(:,1),U(idim,is),V(is,idim)
	    else if(Li%ndim==3)then
		    write(10,200)idim,Mi%atom,Mform,Lform,Li%r(:,1),Oform,Oi%parity,Oi%pair(:,1),U(idim,is),V(is,idim)
	    end if
	  end do; 
    write(10,*)' ' 
    if(imode==1)then
      write(11,*)imode,si(is)
      write(11,'(2f20.6)')U(:,is)
      call GapForm_Band(U(:,is),1)
    end if
	end do; 
  close(10)
  close(11)
100 format(1x, i4, i3, ' ', a3, ' L-', a3, 2f5.1, ' O-', a3,a1, ' ', 2i2, 4f9.4)
200 format(1x, i4, i3, ' ', a3, ' L-', a3, 3f5.1, ' O-', a3,a1, ' ', 2i2, 4f8.4)

  return
end subroutine pairchannel



subroutine RealSpaceTexture(id,U,q,info)
  use workplace
  implicit none
  character*3 id
  complex*16 U(ndim)
  real*8 q(2)
  integer info

  type (Mformconfig), pointer :: Mi
  integer i1,idim,ix,iy
  real*8 S,u1,u2,ra(2),Rt(2)
  complex*16 phase

  if(info==0)then
    open(10,file=id//'mode.dat')
    read(10,*)i1,S
	read(10,*)q
    do idim=1,ndim
       read(10,*)u1,u2
       U(idim)=cmplx(u1,u2)
    end do; close(10)
  end if

  open(10,file=id//'text.dat')
  do idim=1,ndim; Mi=>Mformfunc(idim); ra=model%ra(1:2,Mi%atom)
  do iy=-10,10; do ix=-10,10; rt=ix*model%a(1:2)+iy*model%b(1:2)
     phase=exp(one*pi*sum(q*rt))
	 write(10,100)rt+ra,U(idim)*phase
  end do; end do; end do
  close(10)
100 format(4f15.6)
  return
end subroutine RealSpaceTexture


subroutine GapForm_Band(U,info)
   use workplace
   implicit none

   complex*16 U(ndim)
   integer info

   integer i1,iS,idim,i
   real*8 u1,u2,S,Smax,factor

   integer ikx,iky,ik,iband,error
   real*8 vk(2),Q
   complex*16 pair(norb,norb),hk(norb,norb),ak(2*norb,2*norb)
   real*8 eval(norb),bdgeval(2*norb),bdggap  
   real*8 eta,dos(norb)

   one=cmplx(0,1); pi=asin(1.d0)*2; eta=0.1

   if(info==0)then
     open(20,file='pairmode.dat')
 	 read(20,*)i1,S
	 do idim=1,ndim
	    read(20,*)u1,u2
	    U(idim)=cmplx(u1,u2)
	 end do; close(20)
   end if
   
   open(20,file='fs.dat'); open(21,file='gapform.dat')

   factor=0.1
   do while(.true.); 
     read(20,*,iostat=error)vk
     if(error/=0)exit
     call gethk(norb,vk,Q,hk,model,MF,.false.) !normal state hamiltonian in the orbit and momentum space
     call getgk(vk,pair,U)        !gap function for annihilation operators in the orbit and momentum space

      !bdg gap on fermi surface
	  ak(1:norb,1:norb)=hk; ak(norb+1:2*norb,norb+1:2*norb)=-conjg(hk)
	  ak(1:norb,norb+1:2*norb)=conjg(transpose(pair))*factor; ak(norb+1:2*norb,1:norb)=pair*factor
	  call ZHEIGEN(norb*2,ak,bdgeval)
	  bdggap=1.e10
	  do iband=1,2*norb; bdggap=min(bdggap,abs(bdgeval(iband))); end do

	  call ZHEIGEN(norb,hk,eval)
	  dos=-aimag(1/(one*eta-eval))/pi
      
      !gap function for annihilation operators in the band space
      pair = matmul( conjg(transpose(hk)), matmul(pair,hk) )
  	  iband=1; do i1=2,norb; if(abs(eval(i1))<abs(eval(i1-1)))iband=i1; end do
	  write(21,100)vk,pair(iband,iband),bdggap,eval(iband)
	  !end do
    !end do
	!end do; end do
  end do; close(20); close(21)

100 format(1x,6f12.6)

  print*,'Band-basis gap, BdG gap and meshgap on fermi surfaces output to file'

  return
end subroutine GapForm_Band

subroutine bubbleSusceptibility()
  use workplace
  implicit none

  real*8 charge(norb*2)
  real*8 eta,Agap
  complex*16 U(ndim)
  real*8, allocatable, dimension (:,:) :: q,v(:)
  complex*16, allocatable, dimension (:,:) :: Xq
  complex*16, allocatable, dimension (:,:,:,:) :: Xsdw
  complex*16, allocatable, dimension (:,:,:) :: Vsdw
  complex*16 Vcdw(norb,norb)

  integer mq,mv
  character*20 explain
  real*8 u1,u2
  integer idim,iidim,iq,iv
  integer iatom,iorb
  integer mode
  logical logkmesh,hexagonal
  real*8, dimension (2) :: G,M,X
  
  complex*16 work(nambu,nambu),Ua(nambu,nambu),Va(nambu,nambu)
  real*8 sig(nambu),si(nambu)
  integer info,ia

  !purpose: calculate the dynamic susceptibility for a given vertex charge,
  !         applicable for both spin and charge susceptibilites, provided
  !         that the vertex charge is defined properly.


  !define the vertex charge:
  !-----------optical mode of spin suscept for 2-atom unit cell----------
  mode=1; norb1atom=norb/natom; charge=1
  do iatom=1,natom; do iorb=1,norb1atom; idim=iorb+(iatom-1)*norb1atom
	 if(mod(iatom,2)==0)charge(idim)=mode
  end do; end do
  !------------end of vertex 'charge'--------------


  !-----------  pairing form ---------------------
  Agap=0.2       !gap amplitude multiplying the form factor
  open(10,file='pairmode.dat')
  read(10,*)explain
  do idim=1,ndim
     read(10,*)u1,u2
     U(idim)=cmplx(u1,u2)*Agap
  end do; close(10)
  !------------end of pairing form-----------------


  !---------- define the q-cut for X(q,v)-----------------
  mq=31; allocate(q(2,mq))
  G=0; M=(/1,1/); X=(/1,0/)
  hexagonal=( abs(sum(model%a*model%b))>1.e-5 )
  if(hexagonal)then
    M=(/2./3,2./sqrt(3.)/); X=(/1.,1./sqrt(3.)/)
  end if
  if(natom==2)X=(/2,0/)

  do iq=1,10; q(:,iq)=G+0.1*(iq-1)*(M-G); end do
  do iq=11,20; q(:,iq)=M+0.1*(iq-11)*(X-M); end do
  do iq=21,31; q(:,iq)=X+0.1*(iq-21)*(G-X); end do
  !---------- end of q-cut ------------------------------

  !----------define the energies v for X(q,v) -------
  mv=51; allocate(v(mv))
  do iv=1,mv; v(iv)=(iv-1)*1./(mv-1); end do
  !mv=1; allocate(v(mv)); v=0  
  !-----------end of v-------------------------------

  logkmesh=.false.      !use logarithmic k-mesh in bubble calculation
  nambu=norb*2          !Nambu dimension in supercond state
  allocate(Xsdw(nambu,nambu,mv,mq),Xq(mv,mq))
  allocate(Vsdw(nambu,nambu,mq))
  !call getXsdw(mq,q,mv,v,U,charge,Xsdw,Xq,logkmesh)
  
  open(10,file='XsdwBdG.dat')
    read(10,*)Xsdw
  close(10)

  do iq=1,mq
    do iv=1,mv
	  work=Xsdw(:,:,iv,iq)
	  call ZSVD(nambu,nambu,work,Ua,sig,Va)
	  do ia=1,nambu
	   Si(ia)=Sig(ia)*real( sum( Ua(:,ia)*Va(ia,:) ) )
	  end do
	  Si=-Si;	call ascendingorder(nambu,Si,info); Si=-Si
	  Xq(iv,iq)=cmplx(si(1),0)
	end do
  end do
    
  open(10,file='Xq.dat')
  do iq=1,mq; do iv=1,mv
     write(10,100)v(iv),q(:,iq),Xq(iv,iq)
  end do; end do
  close(10)
100 format(1x,5f15.6)

  open(11,file='XsdwBdG.dat')
    write(11,*)Xsdw
  close(11)

  !RPA spin resonance
  Vsdw=0
  do iq=1,mq
    call p2pinteraction(q(:,iq),Vsdw(1:norb,1:norb,iq),Vcdw)
    if(nambu==2*norb)Vsdw(norb+1:nambu,norb+1:nambu,iq)=Vsdw(1:norb,1:norb,iq)
  end do
  Vsdw=Vsdw
  
  print*,'Vsdw OK'

  call RPA_SC(mv,mq,Vsdw)
  
  return
end subroutine bubbleSusceptibility

subroutine RPA_SC(mv,mq,Vsdw)
  use workplace
  implicit none
  
  integer mv,mq
  complex*16 Xsdw(nambu,nambu,mv,mq),Vsdw(nambu,nambu,mq)

  integer iq,iv,ia,ib
  complex*16 Xq


  Xsdw=Xsdw+1.e-8

  open(10,file='Xqrpa.dat')
  do iq=1,mq
    do iv=1,mv
	  print*,iq,iv
	  call ZINVERT(nambu,Xsdw(:,:,iv,iq))
      
	  Xsdw(:,:,iv,iq)=Xsdw(:,:,iv,iq)+Vsdw(:,:,iq)
	  call ZINVERT(nambu,Xsdw(:,:,iv,iq))
	  Xq=0
	  do ia=1,nambu
	  do ib=1,nambu
	    Xq = Xq + Xsdw(ia,ib,iv,iq)
	  end do
	  end do
	  write(10,'(2e16.8)')Xq
	end do
  end do
  close(10)
   
  return
end subroutine

subroutine getXsdw(mq,q,mv,v,U,charge,Xsdw,Xq,logkmesh)
  use workplace
  implicit none

  integer mq,mv
  real*8 q(2,mq),v(mv),charge(nambu)
  complex*16 U(ndim)
  complex*16 :: Xsdw(nambu,nambu,mv,mq),Xq(mv,mq)
  logical logkmesh

  integer mx,my,ix,iy,icount,mk
  real*8, allocatable, dimension (:,:) :: k
  real*8 eta,Qc
  complex*16, dimension (nambu,nambu) :: Akq,Ak
  real*8, dimension (nambu) :: evkq,evk
  complex*16, dimension (norb,norb) :: hk,gk

  integer iq,ik,ia,ib,iv
  complex*16, dimension (mv) :: Linhard
  real*8 coherence,ea,eb
  real*8 fermi
  real*8 dos(2*mv),w(2*mv)


  
  Qc=0

  if(.not.logkmesh)then
    mx=30
	my=mx
	mk=2*mx*2*my
	eta=0.5/min(mx,my)
    allocate(k(3,mk))
    icount=0
    do iy=-my,my-1; do ix=-mx,mx-1; icount=icount+1
      k(1:2,icount)=0.5*ix*model%ka/mx+0.5*iy*model%kb/my
    end do; end do
    k(3,:)=1./mk
  else
    mk=nk; allocate(k(3,nk)); k=kv; eta=0.01
  end if
  
  Xsdw=0; Xq=0
  do iq=1,mq
    do ik=1,mk
      !BdG eigenstates at k
      call gethk(norb,k(1:2,ik),Qc,hk,model,MF,.false.)
      Ak(1:norb,1:norb)=hk
      if(nambu==2*norb)then
        call getgk(k(1:2,ik),gk,U)
        Ak(norb+1:nambu,norb+1:nambu)=-hk   !assuming h(-k)^*=h(k)
        Ak(1:norb,norb+1:nambu)=transpose(conjg(gk))
        Ak(norb+1:nambu,1:norb)=gk
      end if
      call ZHEIGEN(nambu,Ak,evk)

      !BdG eigenstates at k+q
      call gethk(norb,k(1:2,ik)+q(1:2,iq),Qc,hk,model,MF,.false.)
      Akq(1:norb,1:norb)=hk
      if(nambu==2*norb)then
        call getgk(k(1:2,ik)+q(:,iq),gk,U)
        Akq(norb+1:nambu,norb+1:nambu)=-hk   
        Akq(1:norb,norb+1:nambu)=transpose(conjg(gk))
        Akq(norb+1:nambu,1:norb)=gk
      end if
      call ZHEIGEN(nambu,Akq,evkq)
     
	  !- Tr[G(k+q)*C*G(k)*C], where C is the vertex
      do ia=1,nambu; do ib=1,nambu
        ea=evkq(ia); eb=evk(ib)
		if(abs(ea-eb)<1.e-8)eb=eb+1.e-8
        Linhard=( fermi(ea,1.e-3)-fermi(eb,1.e-3) ) / ( v+one*eta - (ea-eb) )
		coherence=abs( sum(conjg(Akq(:,ia))*charge*Ak(:,ib)) )**2      !coherence factor
        Xsdw(ia,ib,:,iq) = Xsdw(ia,ib,:,iq) + coherence*Linhard*k(3,ik)
     end do; end do
    
    end do
  end do
 

  return
end subroutine

  

subroutine getXq(mq,q,mv,v,U,charge,Xq,logkmesh)
  use workplace
  implicit none

  integer mq,mv
  real*8 q(2,mq),v(mv),charge(nambu)
  complex*16 U(ndim)
  complex*16 :: Xq(mv,mq)
  logical logkmesh

  integer mx,my,ix,iy,icount,mk
  real*8, allocatable, dimension (:,:) :: k
  real*8 eta,Qc
  complex*16, dimension (nambu,nambu) :: Akq,Ak
  real*8, dimension (nambu) :: evkq,evk
  complex*16, dimension (norb,norb) :: hk,gk

  integer iq,ik,ia,ib,iv
  complex*16, dimension (mv) :: Linhard
  real*8 coherence,ea,eb
  real*8 fermi
  real*8 dos(2*mv),w(2*mv)

  Qc=0

  do iv=1,mv; w(iv)=-v(mv-iv+1); end do; w(mv+1:2*mv)=v
  
  if(.not.logkmesh)then
    mx=100; my=100; mk=2*mx*2*my; eta=0.5/min(mx,my)
    allocate(k(3,mk))
    icount=0
    do iy=-my,my-1; do ix=-mx,mx-1; icount=icount+1
      k(1:2,icount)=0.5*ix*model%ka/mx+0.5*iy*model%kb/my
    end do; end do
    k(3,:)=1./mk
  else
    mk=nk; allocate(k(3,nk)); k=kv; eta=0.01
  end if

  Xq=0; dos=0
  do iq=1,mq
    do ik=1,mk
      !BdG eigenstates at k
      call gethk(norb,k(1:2,ik),Qc,hk,model,MF,.false.)
      Ak(1:norb,1:norb)=hk
      if(nambu==2*norb)then
        call getgk(k(1:2,ik),gk,U)
        Ak(norb+1:nambu,norb+1:nambu)=-hk   !assuming h(-k)^*=h(k)
        Ak(1:norb,norb+1:nambu)=transpose(conjg(gk))
        Ak(norb+1:nambu,1:norb)=gk
      end if
      call ZHEIGEN(nambu,Ak,evk)

      if(iq==1)then
        do ib=1,nambu
          dos=dos-k(3,ik)*sum( abs(Ak(1:norb,ib))**2 ) * aimag( 1./(w+one*eta-evk(ib)) ) /pi
        end do
      end if  

      !BdG eigenstates at k+q
      call gethk(norb,k(1:2,ik)+q(1:2,iq),Qc,hk,model,MF,.false.)
      Akq(1:norb,1:norb)=hk
      if(nambu==2*norb)then
        call getgk(k(:,ik)+q(:,iq),gk,U)
        Akq(norb+1:nambu,norb+1:nambu)=-hk   
        Akq(1:norb,norb+1:nambu)=transpose(conjg(gk))
        Akq(norb+1:nambu,1:norb)=gk
      end if
      call ZHEIGEN(nambu,Akq,evkq)
     
      !- Tr[G(k+q)*C*G(k)*C], where C is the vertex
      do ia=1,nambu; do ib=1,nambu
        coherence=abs( sum(conjg(Akq(:,ia))*charge*Ak(:,ib)) )**2      !coherence factor
        ea=evkq(ia); eb=evk(ib)
		if(abs(ea-eb)<1.e-8)eb=eb+1.e-8
        Linhard=( fermi(ea,1.e-3)-fermi(eb,1.e-3) ) / ( v+one*eta - (ea-eb) )
        Xq(:,iq) = Xq(:,iq) + coherence*Linhard*k(3,ik)
      end do; end do

    end do
  end do

  if(nambu==2*norb)then
    open(10,file='dosSC.dat')
  else 
    open(10,file='dosn.dat')
  end if
  do iv=1,2*mv; write(10,*)w(iv),dos(iv); end do; close(10)
  
  return
end subroutine getXq

subroutine getgk(k,pair,U)
  use workplace
  implicit none
  real*8 k(2)
  complex*16 pair(norb,norb)
  complex*16 U(ndim)

  type (Mformconfig), pointer :: Mi
  type (Lformconfig), pointer :: Li
  type (Oformconfig), pointer :: Oi
  integer iOform,iLform,npi,ip,iorb1,iorb2,idim
  complex*16 fki,formk
  real*8 fOi

   pair=0  
   do idim=1,ndim; Mi=>Mformfunc(idim)
       iOform=Mi%Oform; iLform=Mi%Lform
       Oi=>Oformfunc(iOform); Li=>Lformfunc(iLform)
       npi=Oi%npair; do ip=1,npi; fOi=Oi%f(ip)
       iorb1=Oi%pair(1,ip); iorb2=Oi%pair(2,ip)
       if(natom>1)call modifyindex(iorb1,Mi%atom,iorb2,Mi%nbatom)	   
       fki=formk(k,iLform)
    
	   pair(iorb1,iorb2) = pair(iorb1,iorb2) + fki*fOi*U(idim)

    end do; end do

    !gap function for annihilation operators in the orbit and momentum space
    pair = conjg( transpose(pair) )      

	return
end subroutine getgk


function findmax(ndim,nq,V,info)
  implicit none 
  integer ndim,nq,info(3)
  complex*16 V(ndim,ndim,nq)
  real*8 findmax

  integer iq,idim,jdim
 
  findmax=0; info=1
  do iq=1,nq; do jdim=1,ndim; do idim=1,ndim
     if(abs(V(idim,jdim,iq))>findmax)then
       findmax=abs(V(idim,jdim,iq)); info=(/idim,jdim,iq/)
     end if
  end do; end do; end do

  if(abs(findmax)>4.e2)info=-info
  return
end function findmax

function fermi(ek,temp)
  implicit none

  real*8 fermi,ek,temp
  
	fermi=1/(1+exp(-abs(ek)/temp))
  if(ek>0)fermi=1-fermi

  return
end function fermi

subroutine BandStructure(norb,Q,model)
  use standard_derived_types
  implicit none
  integer norb
  real*8 Q
  type (modelblock) :: model
	type (mfblock) :: mf

  integer i
  real*8 pi,kv(2),dk(2)
  complex*16 hk(norb,norb),one
  real*8 ek(norb)
  real*8, dimension (2) :: G,M,X
  logical hexagonal

  pi=asin(1.d0)*2; Q=0
 
  open(10,file='band.dat')
  
  G=(/0,0/); M=(/1.,1./); X=(/1.,0./)
  hexagonal=( abs(sum(model%a*model%b))>1.e-5 )
  if(hexagonal)then
    M=(/2./3,2./sqrt(3.)/); X=(/1.,1./sqrt(3.)/)
  end if

  !cut-1: G-M
  dk=M/100
  do i=0,99; kv=dk*i
     call gethk(norb,kv,Q,hk,model,mf,.false.)
	 call ZHEIGEN(norb,hk,ek)
     write(10,100)ek
  end do

  !cut-2: M-X
  dk=(X-M)/100
  do i=0,99; kv=M+dk*i
     call gethk(norb,kv,Q,hk,model,mf,.false.)
	 call ZHEIGEN(norb,hk,ek)
     write(10,100)ek
  end do

  !cut-3: X-G
  dk=(G-X)/100
  do i=0,100; kv=X+dk*i
     call gethk(norb,kv,Q,hk,model,mf,.false.)
	 call ZHEIGEN(norb,hk,ek)
     write(10,100)ek
  end do

  print*,'band structure output to file'
  close(10)
100 format(1x,f15.6)
  return
end subroutine BandStructure

subroutine RPAXq()
  use workplace
  implicit none

  !input arguments
  integer mq      !number of momenta
  integer mv      !number of frequencies
  real*8 As       !scaling factor for spin interaction (to avoid SDW divergence)
  real*8 Ac       !scaling factor for charge interaction (to avoid CDW divergence)
  real*8 Agap     !scaling factor for pairing form factor
  logical logkmesh           !decide whether logarithmic or linear k-mesh is to be used in the bubble calculation  


  !internal arguments
  integer mx,my,mk,ix,iy,icount,iq,iv,idim
  real*8 eta      !smearing factor
  real*8 Qc       !conserved momentum aside from the in-plane momenta, used only in multi-layer models with periodic condition 

  real*8, allocatable, dimension (:,:) :: k             !k-mesh
  real*8, allocatable, dimension (:,:) :: q             !q-mesh
  real*8, allocatable, dimension (:) :: v               !v-mesh
  complex*16 pairform(ndim)                             !pairing form factor

  complex*16, target, allocatable, dimension (:,:,:) :: Xsq,Xcq      !bare spin- and charge-suseptibilities
  complex*16, dimension (norb,norb) :: RXsq,RXcq                     !work arrays for RPA spin- and charge-suseptibilities
  complex*16, dimension (norb,norb) :: Vsdw,Vcdw                     !work arrays for p2p interactions in SDW and CDW channels

  real*8, dimension (2) :: vtxs,vtxc  !spin and charge vertex

  real*8 u1,u2
  character*20 explain
  logical hexagonal
  real*8, dimension (2) :: G,M,X

  integer is,info,dmq
  real*8 sig(norb),sigr(norb)
  complex*16 UL(norb,norb),UR(norb,norb),work(norb,norb)
  logical appendx0file

  As= 0.1104
  Ac=5.e-2
  Agap=0.05
  logkmesh=.true.
  appendX0file=.true.

  !prepare pairing form
  open(10,file='pairmode.dat')
  read(10,*)explain
  do idim=1,ndim
     read(10,*)u1,u2
     pairform(idim)=cmplx(u1,u2)*Agap
  end do; close(10)

  !---------- define the q-cut for X(q,v)-----------------
  if(.false.)then
    dmq=10
	mq=3*dmq+1
	allocate(q(2,mq))
	G=(/0.0,0.0/)
	M=(/1.0,1.0/)
	X=(/1.0,0.0/)
	if(natom==2.or.natom==4)X=(/2.0,0.0/)
    do iq=1,dmq
	  q(:,iq)=G+(iq-1)*(M-G)/dmq
	end do
	do iq=dmq+1,2*dmq
	  q(:,iq)=M+(iq-dmq-1)*(X-M)/dmq
	end do
	do iq=2*dmq+1,mq
	  q(:,iq)=X+(iq-2*dmq-1)*(G-X)/dmq
	end do
  else
    !mq=nq; allocate(q(2,mq)); q=qv(1:2,:)
    mq=1; allocate(q(2,mq)); q(:,1)=(/1.d0,1.d0/)
  end if
  !---------- end of q-cut ------------------------------

  !define v-mesh and w-mesh
  if(.true.)then
    mv=101; allocate(v(mv)); do iv=1,mv; v(iv)=(iv-1)*0.1/(mv-1); end do
  else
    mv=1; allocate(v(mv)); v=0
  end if
  
  allocate(Xsq(norb,norb,mv),Xcq(norb,norb,mv))

  !prepare k-mesh:
  if(.not.logkmesh)then
    mx=100
	my=mx
	mk=2*mx*2*my
	eta=1./min(mx,my)
	allocate(k(3,mk))
	icount=0
	do iy=-my,my-1; do ix=-mx,mx-1; icount=icount+1
	   k(1:2,icount)=0.5*ix*model%ka/mx+0.5*iy*model%kb/my
	end do; end do
	k(3,:)=1./mk
  else
    mk=nk; allocate(k(3,nk)); k=kv; eta=0.01
  end if

  open(12,file='qw.dat')    
  do iq=1,mq; do iv=1,mv
    write(12,'(i4,1e16.8)')iq,v(iv)
  end do; end do
  close(12); 

  open(12,file='X0file.dat')
  open(10,file='Xs.dat')
  open(11,file='Xs_.dat')
  do iq=1,mq; print*,'iq=',iq
     
	 !prepare orb-orb Vsdw and Vcdw 
     call p2pinteraction(q(:,iq),Vsdw,Vcdw)

	 !scale Vsdw and Vcdw to avoid divergence or unphysical results
	 Vsdw=Vsdw*As; Vcdw=Vcdw*Ac
     
	 !print*,Vsdw; pause; cycle	 
	 
	 if(appendx0file)then
	   read(12,'(2e16.8)')Xsq,Xcq
	 else
	   !calculate bare suceptibilities
       call Orb2OrbX0q(q(:,iq),eta,mv,v,mk,k,pairform,Xsq,Xcq)
	   write(12,'(2e16.8)')Xsq,Xcq
     end if

	 !RPA operations:
	 do iv=1,mv      
	    RXsq=Xsq(:,:,iv)*2; call ZINVERT(norb,RXsq)           !notice the factor of 2 after Xsq leads to <S+S->, which is treated by RPA 
		RXsq=RXsq+Vsdw; call ZINVERT(norb,RXsq); RXsq=RXsq/2  !transform to RPA renormalized <sz*sz>, or <sx*sx>=<sy*sy> 
		RXcq=Xcq(:,:,iv);   call ZINVERT(norb,RXcq)
        RXcq=RXcq+Vcdw;     call ZINVERT(norb,RXcq)
        
        write(10,100)v(iv),q(:,iq),trace(norb,Xsq(:,:,iv)),trace(norb,RXsq)
		
		work=Xsq(:,:,iv)
		call ZHEIGEN(norb,work,sig)
		work=RXsq
        call ZHEIGEN(norb,work,sigr)
		write(11,'(2e16.8)')sig(norb),sigr(norb)
     end do
	 
  end do; close(10); close(11); close(12)
  
100 format(1x,7f11.5)
      
  return
  contains
    function trace(n,A)
	  integer n
	  complex*16 trace,A(n,n)

	  integer i,j
	  trace=0
	  do i=1,n
	  do j=1,n
	    trace=trace+A(i,j)
	  end do
	  end do
	end function
  
end subroutine RPAXq


subroutine Orb2OrbX0q(q,eta,mv,v,mk,k,pairform,Xsq,Xcq)
  use workplace
  implicit none

  integer mv,mk
  real*8 q(2),k(3,mk),v(mv),eta
  complex*16, dimension (norb,norb,mv) :: Xsq,Xcq

  integer ik,nambu1orb
  integer ia,iorb,in,idim
  integer ib,jorb,jn,jdim
  complex*16 pairform(ndim)
  real*8 vtxs(2),vtxc(2),ea,eb
  complex*16, dimension (norb,norb) :: hk,gk
  complex*16, dimension (nambu,nambu) :: Ak,Akq
  real*8, dimension (nambu) :: evk,evkq
  complex*16 cohs,cohc
  complex*16 Linhard(mv),z
  real*8 Qc,fermi
  
  nambu1orb=nambu/norb; Qc=0

  vtxs=1; vtxc=1
  if(nambu1orb==2)vtxc(2)=-1

  Xsq=0; Xcq=0
  do ik=1,mk
    !BdG eigenstates at k
     call gethk(norb,k(1:2,ik),Qc,hk,model,MF,.false.)
     if(nambu==norb)then
	   Ak=hk
     else
	   call getgk(k(1:2,ik),gk,pairform)
       Ak(1:norb,1:norb)=hk
	   Ak(norb+1:nambu,norb+1:nambu)=-hk   !assuming h(-k)^*=h(k)
       Ak(1:norb,norb+1:nambu)=transpose(conjg(gk))
	   Ak(norb+1:nambu,1:norb)=gk
	 end if
	 call ZHEIGEN(nambu,Ak,evk)

     !BdG eigenstates at k+q
	 call gethk(norb,k(1:2,ik)+q,Qc,hk,model,MF,.false.)
     if(nambu==norb)then
	   Akq=hk
	 else
	   call getgk(k(1:2,ik)+q,gk,pairform)
       Akq(1:norb,1:norb)=hk
	   Akq(norb+1:nambu,norb+1:nambu)=-hk   
       Akq(1:norb,norb+1:nambu)=transpose(conjg(gk))
	   Akq(norb+1:nambu,1:norb)=gk
	 end if
	 call ZHEIGEN(nambu,Akq,evkq)
     
	 !- G(k+q)*Vtx*G(k)*Vtx, where Vtx=Vtxs or Vtxc
     do ia=1,nambu; do ib=1,nambu             !sum over band indices
	    ea=evkq(ia)
		eb=evk(ib)
		if(abs(ea-eb)<1.e-6)then
		  eb=ea
		  ea=ea+1.e-3
		end if
		z=cmplx(0,1.e-3)
	    !Linhard=aimag( log( (z-ea)/(z-eb) ) ) /(ea-eb) /pi
		Linhard=( fermi(ea,1.e-3)-fermi(eb,1.e-3) )*k(3,ik) / ( v+one*eta - (ea-eb) )
	    do jorb=1,norb; do jn=1,nambu1orb; jdim=jorb+(jn-1)*norb
		   do iorb=1,norb; do in=1,nambu1orb; idim=iorb+(in-1)*norb
              cohs=conjg(Akq(idim,ib))*vtxs(in)*Ak(idim,ia) * conjg(Ak(jdim,ia))*vtxs(jn)*Akq(jdim,ib)     
              cohc=conjg(Akq(idim,ib))*vtxc(in)*Ak(idim,ia) * conjg(Ak(jdim,ia))*vtxc(jn)*Akq(jdim,ib)     
		      Xsq(iorb,jorb,:) = Xsq(iorb,jorb,:) + cohs*Linhard/4          !notice the factor of 1/4 according to definition of <sz*sz>
		      Xcq(iorb,jorb,:) = Xcq(iorb,jorb,:) + cohc*Linhard            !the factor is 1 for charge susceptibility, by definition of <n*n>
	   end do; end do; end do; end do
	 end do; end do  
  end do  !end of sum over k

  return
end subroutine Orb2OrbX0q


subroutine p2pinteraction(q,Vsdw,Vcdw)
  use workplace
  implicit none
  real*8 q(2)
  complex*16, dimension (norb,norb) :: Vsdw,Vcdw                     

  integer iq
  integer iorb1,iorb2,npi,ip,iOform,idim
  integer jorb1,jorb2,npj,jp,jOform,jdim
  type (Mformconfig), pointer :: Mi,Mj
  type (Oformconfig), pointer :: Oi,Oj
  real*8 fOi,fOj,factor,w
  logical hexagonal
  integer findq

  hexagonal=( abs( sum(model%a*model%b) )>1.e-5 )
  iq=findq(q,hexagonal)   !if(hexagonal)then; call hexsearch(q,iq,qmesh); else; call search(q,iq,qmesh); end if
  
  Vsdw=0; Vcdw=0
  do jdim=1,ndim; Mj=>Mformfunc(jdim); if(Mj%Lform/=1)cycle
  jOform=Mj%Oform; Oj=>Oformfunc(jOform)
  npj=Oj%npair; do jp=1,npj; fOj=Oj%f(jp)
  jorb1=Oj%pair(1,jp); jorb2=Oj%pair(2,jp); if(jorb1/=jorb2)cycle
  if(natom>1)call modifyindex(jorb1,Mj%atom,jorb2,Mj%nbatom)

     do idim=1,ndim; Mi=>Mformfunc(idim); if(Mi%Lform/=1)cycle
	 iOform=Mi%Oform; Oi=>Oformfunc(iOform)
	 npi=Oi%npair; do ip=1,npi; fOi=Oi%f(ip)
	 iorb1=Oi%pair(1,ip); iorb2=Oi%pair(2,ip); if(iorb1/=iorb2)cycle
	 if(natom>1)call modifyindex(iorb1,Mi%atom,iorb2,Mi%nbatom)

        factor=fOi*fOj
	    Vsdw(iorb1,jorb1)=Vsdw(iorb1,jorb1)+C(idim,jdim,iq)*factor
	    Vcdw(iorb1,jorb1)=Vcdw(iorb1,jorb1)+( C(idim,jdim,iq)-2*D(idim,jdim,iq) )*factor
     end do; end do
  end do; end do
  return
end subroutine p2pinteraction

function kernal(x,y)
  use workplace, only : pi,one
  implicit none

  real*8 x,y
  complex*16 kernal

  real*8 f1,f2,del,nu

  del=0.0
  nu=0.0+1.e-8
  
  f1 = ( x*x + 1 ) * ( y*y + 1 )
  f2 = pi * exp( -nu * abs( x - y + del/nu ) ) * ( 1 + nu * abs( x - y + del/nu ) )
  kernal = exp( -one * f2 ) / f1
      
  return
end function

subroutine checkLa()
  implicit none

  integer nx,ix,iy
  real*8 xmax,dx,x,y
  complex*16 kernal

  xmax=100
  nx=100
  dx=xmax/nx

  open(10,file='check.dat')

  do ix=-nx,nx
    x=dx*ix
    do iy=-nx,nx
	  y=dx*iy
      write(10,'(4e16.8)')x,y,kernal(x,y)
	end do
  end do
  close(10)
  
  return
end subroutine

function LaChen(nu,del)
  use workplace, only : pi,one
  implicit none

  real*8 nu,del
  complex*16 LaChen

  integer nx,ix,iy
  real*8 xmax,dx,x,y,f1,f2
  complex*16 f3,kernal
  
  xmax=100
  dx=0.1
  nx=xmax/dx
  LaChen=0
  if(abs(nu)<1.e-8)nu=nu+1.e-8
  do ix=-nx,nx
    x=dx*ix
    do iy=-nx,nx
	  y=dx*iy
      f1 = ( x*x + 1 ) * ( y*y + 1 )
      f2 = pi * exp( -nu * abs( x - y + del/nu) ) * ( 1 + nu * abs( x - y + del/nu) )
      f3 = exp( -one * f2 ) / f1
	  LaChen = LaChen + f3 * dx * dx
	end do
  end do

  LaChen=LaChen/pi/pi
  return
end function

subroutine ChenWei()
  implicit none

  integer nx,ix,iy
  real*8 xmax,dx,x,y,theta,f,p
  complex*16 LaChen,La

  open(10,file='f.dat')
  open(11,file='p.dat')
  open(12,file='la.dat')
  
  xmax=0.5
  nx=10
  dx=xmax/nx
  do ix=0,nx
    x=ix*dx
	do iy=0,nx
	  y=iy*dx
      print*,ix,iy
	  La=LaChen(x,y)
	  theta=asin(aimag(La)/abs(La))+3.1415926
	  print*,(theta)
	  f=(5-3*abs(La)*cos(theta))/8.0
	  p=(5+3*abs(La)**2)/8.0
	  print*,f,p
	  print*,La
	  write(10,*)f
	  write(11,*)p
	  write(12,'(2e16.8)')La
	end do
  end do
  
  close(10)
  close(11)
  close(12)

  return
end subroutine
