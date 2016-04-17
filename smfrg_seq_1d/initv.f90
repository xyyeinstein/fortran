subroutine initV()
  use workplace
  implicit none

  integer idim,iOform,iLform,ip,npi
  integer jdim,jOform,jLform,jp,npj

  type (Oformconfig), pointer :: Oi,Oj
  type (Lformconfig), pointer :: Li,Lj
  type (Mformconfig), pointer :: Mi,Mj
	
  real*8 fOi,fOj,fLi,fLj,addp,addc,addd
  integer orb(4),iorb,jorb
 
  if(model%sample=='CuO2')then; call initV_CuO2(); return; end if

  allocate(P(ndim,ndim,nq),C(ndim,ndim,nq),D(ndim,ndim,nq))
  allocate(Ploc(ndim,ndim),CLoc(ndim,ndim),Dloc(ndim,ndim))
  if(appendbreakpoint)return

  P=0; C=0; D=0

  do jdim=1,ndim; Mj=>Mformfunc(jdim)
  jOform=Mj%Oform; if(Mj%Lform/=1)cycle 
  Oj=>Oformfunc(jOform); npj=Oj%npair
  do jp=1,npj; fOj=Oj%f(jp); orb(3:4)=Oj%pair(:,jp)

     do idim=1,ndim; Mi=>Mformfunc(idim)
	 iOform=Mi%Oform; if(Mi%Lform/=1)cycle 
     if(Mi%atom/=Mj%atom)cycle
     Oi=>Oformfunc(iOform); npi=Oi%npair
	 do ip=1,npi; fOi=Oi%f(ip); orb(1:2)=Oi%pair(:,ip)	 
	 
	 !pairing channel
	 addp=0
	 if(orb(1)==orb(2).and.orb(3)==orb(4))then
	   if(orb(1)==orb(3))then; addp=addp-model%U; else; addp=addp-model%JH; end if
	 end if
	 !contact terms from C and D
	 if(orb(1)==orb(4).and.orb(2)==orb(3).and.orb(1)/=orb(2))addp=addp-model%JH
	 if(orb(1)==orb(3).and.orb(2)==orb(4).and.orb(1)/=orb(2))addp=addp-model%Uab

	 !crossing channel
	 addc=0
	 if(orb(1)==orb(2).and.orb(3)==orb(4))then
	   if(orb(1)==orb(3))then; addc=addc-model%U; else; addc=addc-model%JH; end if
	 end if
	 !contact terms from P and D
	 if(orb(1)==orb(4).and.orb(2)==orb(3).and.orb(1)/=orb(2))addc=addc-model%JH
	 if(orb(1)==orb(3).and.orb(2)==orb(4).and.orb(1)/=orb(2))addc=addc-model%Uab

     !direct channel
	 addd=0
	 if(orb(1)==orb(2).and.orb(3)==orb(4))then
	   if(orb(1)==orb(3))then; addd=addd-model%U; else; addd=addd-model%Uab; end if
	 end if
	 !contact terms from P and C
	 if(orb(1)==orb(4).and.orb(2)==orb(3).and.orb(1)/=orb(2))addd=addd-model%JH
	 if(orb(1)==orb(3).and.orb(2)==orb(4).and.orb(1)/=orb(2))addd=addd-model%JH
	   
	 P(idim,jdim,:)=P(idim,jdim,:)+addp*fOi*fOj
	 C(idim,jdim,:)=C(idim,jdim,:)+addc*fOi*fOj
	 D(idim,jdim,:)=D(idim,jdim,:)+addD*fOi*fOj    

  end do; end do
  end do; end do

  Ploc=P(:,:,1); Cloc=C(:,:,1); Dloc=D(:,:,1)

  if(abs(model%Vnn)>1.e-5)call addVnn()
  if(abs(model%Jnn)>1.e-5)call addJnn()

  return
end subroutine initV


subroutine initV_CuO2()
  use workplace
  implicit none

  integer idim,iOform,iLform,ip,npi
  integer jdim,jOform,jLform,jp,npj

  type (Oformconfig), pointer :: Oi,Oj
  type (Lformconfig), pointer :: Li,Lj
  type (Mformconfig), pointer :: Mi,Mj
	
  real*8 fOi,fOj,fLi,fLj,addp,addc,addd
  integer orb(4),iorb,jorb
  integer iq
  real*8 q(2)
 
  allocate(P(ndim,ndim,nq),C(ndim,ndim,nq),D(ndim,ndim,nq))
  if(appendbreakpoint)return

  P=0; C=0; D=0

  do jdim=1,ndim; Mj=>Mformfunc(jdim)
  jOform=Mj%Oform; if(Mj%Lform/=1)cycle 
  Oj=>Oformfunc(jOform); npj=Oj%npair
  do jp=1,npj; fOj=Oj%f(jp); orb(3:4)=Oj%pair(:,jp); if(orb(3)/=orb(4))cycle

     do idim=1,ndim; Mi=>Mformfunc(idim)
	 iOform=Mi%Oform; if(Mi%Lform/=1)cycle 
     Oi=>Oformfunc(iOform); npi=Oi%npair
	 do ip=1,npi; fOi=Oi%f(ip); orb(1:2)=Oi%pair(:,ip); if(orb(1)/=orb(2))cycle	 
	 
	 do iq=1,nq; q=qv(1:2,iq)
	    addd=0
	    if(orb(1)==orb(3))then
	      if(orb(1)==1)then; addd=addd-model%Udd-2*sum(cos(q*pi))*model%Vdd; else; addd=addd-model%Upp; end if
	    end if

        if((orb(1)==1.and.orb(3)==2).or.(orb(1)==2.and.orb(3)==1))addd=addd-model%Vpd*2*cos(q(1)*pi*0.5)
        if((orb(1)==1.and.orb(3)==3).or.(orb(1)==3.and.orb(3)==1))addd=addd-model%Vpd*2*cos(q(2)*pi*0.5)
		if((orb(1)==2.and.orb(3)==3).or.(orb(1)==3.and.orb(3)==2))addd=addd-model%Vpp*4*cos(q(1)*pi*0.5)*cos(q(2)*pi*0.5)
        D(idim,jdim,iq)=D(idim,jdim,iq)+addd*fOi*fOj
	 end do
  end do; end do
  end do; end do

  call contact(P,C,D)
  return
end subroutine initV_CuO2


subroutine addVnn()
  use workplace
  implicit none

  integer idim,iOform,iLform,ip,npi
  integer jdim,jOform,jLform,jp,npj

  type (Oformconfig), pointer :: Oi,Oj
  type (Lformconfig), pointer :: Li,Lj
  type (Mformconfig), pointer :: Mi,Mj
	
  real*8 fOi,fOj,fLi,fLj,addp,addc,addd
  integer orb(4),iorb,jorb
	 
  complex*16, dimension (ndim,ndim,nq) :: dP,dC,dD
  real*8 q(2),bond(3),basis(3)
  integer i,nbond
  logical hexagonal

  integer findatom

  hexagonal = ( abs(sum(model%a*model%b))>1.e-5 )
  basis=model%a; if(natom>1)basis=model%ra(:,2)-model%ra(:,1)
  nbond=4
  if(hexagonal)nbond=6
  if(sum(abs(model%kb))==0)nbond=2

  dP=0; dC=0; dD=0

  do jdim=1,ndim
    Mj=>Mformfunc(jdim)
    jOform=Mj%Oform
	if(Mj%Lform/=1)cycle 
    Oj=>Oformfunc(jOform)
	npj=Oj%npair
    do jp=1,npj
	  fOj=Oj%f(jp)
	  orb(3:4)=Oj%pair(:,jp)
	  if(orb(3)/=orb(4))cycle
      bond=basis
	  do i=1,nbond 
        if(i>1)then
          if(sum(abs(model%kb))==0)call rotate_g(bond,hexagonal)
          call rotate_g(bond,hexagonal)
	    end if
        do idim=1,ndim
          Mi=>Mformfunc(idim)
          if(findatom(model%ra(:,Mi%atom)+bond)/=Mj%atom)cycle
	      iOform=Mi%Oform
		  if(Mi%Lform/=1)cycle 
          Oi=>Oformfunc(iOform)
	      npi=Oi%npair
	      do ip=1,npi
		    fOi=Oi%f(ip)
	        orb(1:2)=Oi%pair(:,ip)
	        if(orb(1)/=orb(2))cycle

		    dD(idim,jdim,:) = dD(idim,jdim,:) -model%Vnn*fOi*fOj*exp(one*pi*( qv(1,:)*bond(1)+qv(2,:)*bond(2) ) )

	      end do
		end do
	  end do
    end do
  end do

  call contact(dP,dC,dD)
  P=P+dP; C=C+dC; D=D+dD

  return
end subroutine addVnn

subroutine addJnn()
  use workplace
  implicit none

  integer idim,iOform,iLform,ip,npi
  integer jdim,jOform,jLform,jp,npj

  type (Oformconfig), pointer :: Oi,Oj
  type (Lformconfig), pointer :: Li,Lj
  type (Mformconfig), pointer :: Mi,Mj
	
  real*8 fOi,fOj,fLi,fLj,addp,addc,addd
  integer orb(4),iorb,jorb
	 
  complex*16, dimension (ndim,ndim,nq) :: dP,dC,dD
  real*8 q(2),bond(3),basis(3)
  integer i,nbond
  logical hexagonal
  integer findatom

  hexagonal= ( abs(sum(model%a*model%b))>1.e-5 )
  basis=model%a; if(natom>1)basis=model%ra(:,2)-model%ra(:,1)
  nbond=4; if(hexagonal)nbond=6

  dP=0; dC=0; dD=0

  do jdim=1,ndim; Mj=>Mformfunc(jdim)
  jOform=Mj%Oform; if(Mj%Lform/=1)cycle 
  Oj=>Oformfunc(jOform); npj=Oj%npair
  do jp=1,npj; fOj=Oj%f(jp); orb(3:4)=Oj%pair(:,jp); if(orb(3)/=orb(4))cycle

     bond=basis
	 do i=1,nbond; if(i>1)call rotate_g(bond(1:2),hexagonal)
     do idim=1,ndim; Mi=>Mformfunc(idim)
        if(findatom(model%ra(:,Mi%atom)+bond)/=Mj%atom)cycle
	    iOform=Mi%Oform; if(Mi%Lform/=1)cycle 
        Oi=>Oformfunc(iOform); npi=Oi%npair
	    do ip=1,npi; fOi=Oi%f(ip); orb(1:2)=Oi%pair(:,ip); if(orb(1)/=orb(2))cycle

	 	   dC(idim,jdim,:) = dC(idim,jdim,:)-model%Jnn*fOi*fOj*exp(one*pi*( qv(1,:)*bond(1)+qv(2,:)*bond(2) ) )

		end do
	 end do; end do
  end do; end do

  call contact(dP,dC,dD)
  P=P+dP; C=C+dC; D=D+dD

  return
end subroutine addJnn

subroutine rotate_g(bond,hexagonal)
  implicit none
  real*8 bond(2)
  logical hexagonal

  if(hexagonal)then; call rotate60(bond); else; call rotate90(bond); end if

  return
end subroutine rotate_g

subroutine HexaddVL2L
  use workplace
  implicit none

  integer idim,jdim
  
  type (Mformconfig), pointer :: Mi,Mj
  type (Oformconfig), pointer :: Oi,Oj
	 
  complex*16, dimension (ndim,ndim,nq) :: dP,dC,dD
  real*8 q(2),bond(3),basis(3)
  integer i,j,ni,nj,iOform,jOform
  complex*16 add

  dP=0; dC=0; dD=0

  do jdim=1,ndim; Mj=>Mformfunc(jdim); if(Mj%Lform/=1)cycle
  do idim=1,ndim; Mi=>Mformfunc(idim); if(Mi%Lform/=1)cycle
     if(Mi%atom==Mj%atom)cycle
	 if(sum(abs(model%ra(1:2,Mi%atom)-model%ra(1:2,Mj%atom)))>1.e-5)cycle
	 iOform=Mi%Oform; jOform=Mj%Oform
	 Oi=>Oformfunc(iOform); Oj=>Oformfunc(jOform)
	 ni=Oi%npair; nj=Oj%npair
	 do i=1,ni; if(Oi%pair(1,i)/=Oi%pair(2,i))cycle
	 do j=1,nj; if(Oj%pair(1,j)/=Oj%pair(2,j))cycle
        dD(idim,jdim,:)=dD(idim,jdim,:)- 0*Oj%f(j)*Oi%f(i)  !model%VL2L  VL2L not yet implemented
	 end do; end do
  end do; end do

  call contact(dP,dC,dD)
  P=P+dP; C=C+dC; D=D+dD

  return
end subroutine HexaddVL2L
