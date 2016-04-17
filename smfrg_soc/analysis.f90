subroutine FermiSurf(norb,model)
  use standard_derived_types
  implicit none
  integer norb
  type (modelblock) :: model

  integer ix,iy,ia,ispin,iispin,idim,iidim,isign,natom,iband
  real*8 pi,kv(2),Sx,Sy
  complex*16 hk(norb,norb),one
  real*8 ek(norb),eta,dos
  logical propergauge

  pi=asin(1.d0)*2; one=cmplx(0,1); eta=0.05
  natom=model%natom

  propergauge=model%propergauge
  model%propergauge=.true. 
  open(10,file='Ak.dat')
  do iy=-40,40;  do ix=-40,40; kv=(/ix,iy/)*1./40
     call gethk(norb,kv,hk,model)
     call ZHEIGEN(norb,hk,ek)
	 do iband=1,norb
	    Sx=0; Sy=0
		do ispin=1,2; do ia=1,model%natom; idim=ia+(ispin-1)*natom
		   iispin=2; if(ispin==2)iispin=1
		   isign=1; if(ispin>iispin)isign=-1
		   iidim=ia+(iispin-1)*natom
		   Sx=Sx+real( conjg(hk(idim,iband))*hk(iidim,iband) )
		   Sy=Sy-aimag( conjg(hk(idim,iband))*hk(iidim,iband)*isign )
		end do; end do
	    dos=-aimag( 1/(one*eta-ek(iband)) )/pi
		write(10,100)dos,Sx,Sy
	 end do
  end do; end do
100 format(1x,3f20.6)
  close(10)
  model%propergauge=propergauge
  return
end subroutine FermiSurf

subroutine BandStructure(norb,model)
  use standard_derived_types
  integer norb
  type (modelblock) :: model

  integer i
  real*8 pi,kv(2),K(2),M(2),G(2),dk(2)
  complex*16 hk(norb,norb),one
  real*8 ek(norb)

  pi=asin(1.d0)*2
 
  open(10,file='band.dat')
  
  G=(/0,0/); M=(/1.,1./); K=(/1.,0./)

  !cut-1: G-M
  dk=M/100
  do i=0,99; kv=dk*i
     call gethk(norb,kv,hk,model)
	 call ZHEIGEN(norb,hk,ek)
     write(10,100)ek
  end do

  !cut-2: M-K
  dk=(K-M)/100
  do i=0,99; kv=M+dk*i
     call gethk(norb,kv,hk,model)
	 call ZHEIGEN(norb,hk,ek)
     write(10,100)ek
  end do

  !cut-3: K-G
  dk=(G-K)/100
  do i=0,100; kv=K+dk*i
     call gethk(norb,kv,hk,model)
	 call ZHEIGEN(norb,hk,ek)
     write(10,100)ek
  end do

  print*,'band structure output to file'
  close(10)
100 format(1x,4f16.6)
  return
end subroutine BandStructure


subroutine analysis(la,isave)
  use workplace
  implicit none
  real*8 la
  logical isave

  integer i,info(3),iq(2),iS
  real*8 q(2),S(2),Vmax,Si,u1,u2
  integer ig,iqg
  real*8 findmax
  complex*16 work(ndim,ndim)

  real*8 t1,t2

  Vmax=findmax(ndim,nq,P,info); iq(1)=abs(info(3))
  !output%sc=info(1:2);

  q=0; call search(q,iq(1),qmesh); i=iq(1)
  work=P(:,:,i)
  if(usegroup)then
    work=0
    do ig=1,ng; iqg=indexqg(ig,i)
      work=work+P(:,:,iqg)
    end do
    work=work/ng
  end if
  work=work*0.5
  call channel(la,'P',q,ndim,work,S(1),isave)
  open(37,file='ppmode.dat')
  read(37,*)iS,Si
  iS=1;Si=0
  do i=1,ndim
    read(37,*)u1,u2
    if(u1*u1+u2*u2>Si)then
      Si=u1*u1+u2*u2
      iS=i
    end if
  end do
  close(37)
  output%sc=iS    
  
  Vmax=max(Vmax,findmax(ndim,nq,-C,info)); 
  S(2)=0; iq(2)=1;
  do i=1,nq
    call channel(la,'D',qv(1:2,i),ndim,-C(:,:,i),Si,.false.)
    if(Si>S(2))then
      S(2)=Si
      iq(2)=i
    end if
  end do
  
  call channel(la,'D',qv(1:2,iq(2)),ndim,-C(:,:,iq(2)),Si,isave)
  call phchannel(la,qv(1:2,iq(2)),ndim,-C(:,:,iq(2)))

  !q=0; call search(q,i,qmesh); 
  !call channel('D',q,ndim,-C(:,:,i),output%fm,isave)

  output%Vmax(1)=Vmax; output%Vmax(2:3)=S
  output%q(:,1)=0; output%q(:,2)=qv(1:2,iq(2))

  return
end subroutine analysis

  subroutine channel(la,name,q,nA,A,Smax,isave)
    use workplace, only : model,Mformconfig,Lformconfig,Mformfunc,one,pi
    implicit none
    character*1 name
  	real*8 q(2),la
    integer nA,iS,iSmax
	  real*8 Smax,Si
    complex*16, dimension (nA,nA) :: A
	  logical isave

    complex*16, dimension (nA,nA) :: B,U,V
	  real*8 Sig(nA)
    integer idim,jdim,iform

    character*3 Mform,Lform
    character*2 Sform
    type (Mformconfig), pointer :: Mi
    type (Lformconfig), pointer :: Li
	  real*8, pointer, dimension (:,:) :: ra
    complex*16 phase

    B=A; call ZSVD( nA,nA,B,U,Sig,V )

    Smax=-1.d10; iSmax=0
    do iS=1,nA
      Si=Sig(iS)*real( sum( U(:,iS)*V(iS,:) ) )
      if(Si>Smax)then; Smax=Si; iSmax=iS; end if
    end do
    
    !gauge transform to proper gauge
    do idim=1,nA
      Mi=>Mformfunc(idim); ra=>model%ra; phase=exp(one*pi*sum(q*ra(1:2,Mi%atom)))
      U(idim,:)=U(idim,:)*phase; V(:,idim)=V(:,idim)*conjg(phase)
    end do
    
    if(name=='P')then
      open(11,file='ppmode.dat')
        write(11,*)iSmax,Sig(iSmax)*real( sum( U(:,iSmax)*V(iSmax,:) ) )
        write(11,300)U(:,iSmax)
        call GapForm_Band(U(:,iSmax),1)
      close(11)
    end if

    if(.not.isave)return

    if(name=='P')then
      open(10,file='ppMode.txt')
    else if(name=='D')then
      open(10,position='append',file='phMode.txt')
    end if

    write(10,*)'energy scale=',la
    write(10,*)'  ordering momentum/pi='; write(10,*)q
    write(10,*)' '
    do iS=1,nA; write(10,*)'   iS, S =',iS,Sig(iS)*real( sum( U(:,iS)*V(iS,:) ) )
      write(10,*)'  iM ia   S  Lform    rbasis       Uvec              Vvec'
      do idim=1,nA
        Mi=>Mformfunc(idim); Mform=Mi%form; Sform=Mi%Sform
        Li=>Mi%Lformfunc; Lform=Li%form
        write(10,100)idim,Mi%atom,Sform,Lform,Li%r(:,1),U(idim,iS),V(iS,idim)
      end do
      write(10,*)' ' 
    end do; close(10)
100 format(1x, i4, i3, 2(' ',a3,' '),3f5.1, 4f9.4)

    if(name=='D')then 
      open(11,file='phmode.dat')
        write(11,*)iSmax,Sig(iSmax)*real( sum( U(:,iSmax)*V(iSmax,:) ) )
        write(11,*)q
        write(11,300)U(:,iSmax)
      close(11)
      call RealSpaceTexture('phd',U(:,iSmax),q,1)         
    end if
300 format(1x,2f20.6)
    return
end subroutine channel

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
    open(10,file='phmode.dat')
    read(10,*)i1,S
    read(10,*)q
    do idim=1,ndim
       read(10,*)u1,u2
       U(idim)=cmplx(u1,u2)
    end do; close(10)
  end if

  open(10,file='phtext.dat')
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

   integer i1,iS,idim
   real*8 u1,u2,S,Smax

   integer ikx,iky,ik,iband
   real*8 vk(2),Q
   complex*16 pair(norb,norb),hk(norb,norb),hkbar(norb,norb)
   real*8 eval(norb),evalmin

   type (Mformconfig), pointer :: Mi
   type (Lformconfig), pointer :: Li
   integer iLform,iorb1,iorb2
   real*8 eta,dos(norb)
   complex*16 fki,overlap,formk

   integer ispin,iispin,ia,iidim,isign
   logical propergauge

   integer error
   
   eta=0.1

   if(info==0)then
     open(10,file='ppmode.dat')
 	 read(10,*)i1,S
	 do idim=1,ndim
	    read(10,*)u1,u2
	    U(idim)=cmplx(u1,u2)
	 end do; close(10)
   end if

   open(10,file='gapform.dat')
   open(11,file='fs4cal.dat')
   propergauge=model%propergauge
   model%propergauge=.false.   !improper gauge forced since fk is in this gauge
   do while(.true.)
	 read(11,*,iostat=error)vk
	 if(error/=0)exit
   !do iky=-40,40; do ikx=-40,40; vk=(/ikx,iky/)*1./40.
   !do ik=1,nk; vk=kv(1:2,ik)
      call gethk(norb,vk,hk,model)
	  call ZHEIGEN(norb,hk,eval)
	  dos=-aimag(1/(one*eta-eval))/pi
	  
	  !prepare transpose of U(-k), computed from the Kramers pair:
	  !   |-k> = sig_2 conjg( |k> ),  sig_2 = 2nd Pauli matrix in spin space
      do iband=1,norb
	     do ispin=1,2; iispin=2; if(ispin==2)iispin=1
            isign=1; if(ispin>iispin)isign=-1
            do ia=1,natom; idim=ia+(ispin-1)*natom; iidim=ia+(iispin-1)*natom
			   hkbar(idim,iband)=conjg(hk(iidim,iband))*isign
		 end do; end do
	  end do
	  hkbar=transpose(hkbar)
	  
	  pair=0  
      do idim=1,ndim; Mi=>Mformfunc(idim)
         iLform=Mi%Lform; Li=>Mi%Lformfunc
         iorb1=Mi%orb; iorb2=Mi%nborb
         fki=formk(vk,Li)     
	     pair(iorb1,iorb2) = pair(iorb1,iorb2) + fki*U(idim)
	  end do

	  
	  pair=0; call pairinggk(norb,vk,pair,model)   
      
	  !gap function for annihilation operators in the orbit and momentum space
	  pair = conjg( transpose(pair) )   
	  
      !gap function for annihilation operators in the band space
      pair = matmul( hkbar, matmul(pair,hk) )
	  !do iband=1,norb; write(10,100)vk,pair(iband,iband),dos(iband); end do
      iband=1; evalmin=abs(eval(1))
	  do i1=2,norb; if(abs(eval(i1))<evalmin)then; iband=i1; evalmin=abs(eval(i1)); end if; end do
      write(10,100)vk,pair(iband,iband),dos(iband)
	!end do; end do
	end do
	close(10); close(11)
100 format(1x,5f15.6)
    model%propergauge=propergauge
    return
end subroutine GapForm_Band

subroutine BdG_Band(U,info)
   use workplace
   implicit none

   complex*16 U(ndim)
   integer info

   integer i1,iS,idim
   real*8 u1,u2,S,Smax

   integer ikx,iky,ik,iband
   real*8 vk(2)
   complex*16 pair(norb,norb),hk(norb,norb)
   complex*16 Ak(2*norb,2*norb)
   real*8 eval(2*norb)

   type (Mformconfig), pointer :: Mi
   type (Lformconfig), pointer :: Li
   integer iLform,iorb1,iorb2,norb2
   real*8 eta,dos(norb)
   complex*16 fki,overlap,formk

   integer ispin,iispin,ia,iidim,isign
   logical propergauge
   
   eta=0.1; norb2=2*norb

   if(info==0)then
     open(10,file='ppmode.dat')
 	 read(10,*)i1,S
	 do idim=1,ndim
	    read(10,*)u1,u2
	    U(idim)=cmplx(u1,u2)
	 end do; close(10)
   end if

   open(10,file='BdGek.dat')
   propergauge=model%propergauge
   model%propergauge=.false.   !improper gauge forced since fk is in this gauge
   do iky=-40,40; do ikx=-40,40; vk=(/ikx,iky/)*1./40.
   !do ik=1,nk; vk=kv(1:2,ik)
      call gethk(norb,vk,hk,model); Ak(1:norb,1:norb)=hk
      call gethk(norb,-vk,hk,model); Ak(norb+1:norb2,norb+1:norb2)=-conjg(hk)	  
	  pair=0  
      do idim=1,ndim; Mi=>Mformfunc(idim)
         iLform=Mi%Lform; Li=>Mi%Lformfunc
         iorb1=Mi%orb; iorb2=Mi%nborb
         fki=formk(vk,Li)     
	     pair(iorb1,iorb2) = pair(iorb1,iorb2) + fki*U(idim)
	  end do
      Ak(1:norb,norb+1:norb2)=Pair; Ak(norb+1:norb2,1:norb)=conjg(transpose(pair))
      call ZHEIGEN(norb2,Ak,eval)
	  do iband=1,norb2; write(10,100)vk,eval(iband); end do
   end do; end do
   close(10)
100 format(1x,3f15.6)
   return
end subroutine BdG_Band

subroutine BdGSlab(model,edgemodel)
use standard_derived_types
implicit none
type (modelblock) :: model
type (edgemodelblock) :: edgemodel

!purpose: find the eigen states of a sample with open boundary along normal direction along nb

!arguments imported from model:
!     nb:  normal direction of the open boundary
!     Nq:  number of conserved transverse q-vectors 
!     Lb:  length of the sample slab along the normal direction
!     norb: number of effective orbits in the normal state (including spin)
!     nambu: dimension of the Nambu space
!     ndim:  ndim=nambu*Lb, hilbert space dimension of the open-boundary slice

!notice:
!    if nambu=2*norb, the program deals with a superconducting sample
!    if nambu=norb,   the program deals with a normal-state sample 

integer nb,nq,Lb,norb,nambu,ndim
complex*16, allocatable :: hq(:,:,:),h(:,:)
real*8, allocatable :: eval(:)

integer i,idim
real*8 q(2)
character*5 sample
logical periodic

!internal arguments:
!     hq(:,:,b) = sum_kb hk(kb,q) exp(i*kb*b), effective nambu-space matrix on a bond b along nb-th direction
!     h(:,:) : nambu-space matrix for the open boundary slice along nb-th direction, with conserved q, constructed from hq. 
!     eval(:) : eigen values of h

sample=edgemodel%sample; periodic=edgemodel%periodicslab
nb=edgemodel%nb; nq=edgemodel%nq; Lb=edgemodel%Lb
norb=edgemodel%norb; nambu=edgemodel%nambu; ndim=edgemodel%ndim

allocate(hq(nambu,nambu,-4:4),h(ndim,ndim),eval(ndim))

open(45,file='eq.dat')

do i=-Nq,Nq
   print*,i
   if(nb==1)then; q=0.5*i*edgemodel%kb(1:2)/Nq; else; q=0.5*i*edgemodel%ka(1:2)/Nq; end if
   call gethq(nb,q,norb,nambu,hq,model,edgemodel)
   call gethr(nambu,hq,Lb,ndim,h,periodic)
   call ZHEIGEN(ndim,h,eval)
   do idim=1,ndim
      write(45,100)eval(idim), sum( abs(h(1:nambu,idim))**2 ), sum( abs(h(ndim-nambu+1:ndim,idim))**2 ) 
   end do
end do
!close(15)

100 format(1x,3f20.10)
close(45)
return
end subroutine BdGSlab

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

  return
end function findmax


subroutine outputdos()
  use workplace

  real*8 dos,wL,wU,dw,w,eta
  integer n,i,ix,iy,ival

  wL=-3; wU=3; n=401; dw=(wU-wL)/(n-1)
  eta=1.e-1

  open(10,file='dos.dat')
  do i=1,n; w=wL+(i-1)*dw; dos=0
     do iy=-nky,nky-1; do ix=-nkx,nkx-1
		dos = dos - sum( aimag( 1/(one*eta+w-ek(:,ix,iy)) ) )
	 end do; end do
	 dos=dos/(2*nkx*2*nky*norb)
	 write(10,*)w,dos
  end do
  close(10)

  return
end subroutine outputdos

function dopinglevel()
  use workplace
  implicit none

  integer ix,iy,iorb
  real*8 dopinglevel
  real*8 fermi

  dopinglevel=0
  do iy=-nky,nky-1; do ix=-nkx,nkx-1; do iorb=1,norb
     dopinglevel=dopinglevel+fermi( ek(iorb,ix,iy) )
  end do; end do; end do

  dopinglevel=dopinglevel/(2*nkx*2*nky)
  return
end function dopinglevel

function fermi(ek)
  implicit none
  real*8 fermi,ek

  real*8 temp; data temp/1.d-3/

  fermi=1/(1+exp(-abs(ek)/temp))
  if(ek>0)fermi=1-fermi

  return
end function fermi






