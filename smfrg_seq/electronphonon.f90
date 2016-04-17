subroutine ElectronPhonon(La,Peph,Ceph)
  use workplace
  implicit none
  real*8 La
  complex*16, dimension (ndim,ndim,nq) :: Peph,Ceph,Deph

  type (Mformconfig), pointer :: Mi,Mj
  type (Oformconfig), pointer :: Oi,Oj
  integer idim,iOform,ip,npi,jdim,jOform,jp,npj
  real*8 velocity,Veph,w0,w1,fOi,fOj
  real*8 eq(nq),Dq(nq)

  !notice: this is a temper subroutine for electron phonon coupling. it assumes local A1g phonon with steep dispersion
  !        and an effective attraction Veph independent of the phonon momentum 

  !call ElectronPhonon_2Fe(La,Peph,Ceph); return

  Veph=model%Veph; w0=0.1; velocity=0.025; w1=velocity
  eq=w0*w0+velocity*velocity*(4-2*cos(qv(1,:)*pi)-2*cos(qv(2,:)*pi))   !eq=wq*wq
  Dq=0.5*Veph*sqrt( w1*w1/( eq + La*La ) )
  
  Deph=0

  do jdim=1,ndim; Mj=>Mformfunc(jdim)
  jOform=Mj%Oform; if(Mj%Lform/=1)cycle; if(Mj%form/='A1g')cycle
  Oj=>Oformfunc(jOform); npj=Oj%npair
  do jp=1,npj; fOj=Oj%f(jp)

     do idim=1,ndim; Mi=>Mformfunc(idim)
	 iOform=Mi%Oform; if(Mi%Lform/=1)cycle; if(Mi%form/='A1g')cycle 
     Oi=>Oformfunc(iOform); npi=Oi%npair
	 do ip=1,npi; fOi=Oi%f(ip)	 
	 
     !direct channel
     Deph(idim,jdim,:)=Deph(idim,jdim,:)+Dq*fOi*fOj

  end do; end do
  end do; end do

  Peph=0; Ceph=0
  if(useprojectors)then; call contact_clever(Peph,Ceph,Deph); else; call contact(Peph,Ceph,Deph); end if

  return
end subroutine ElectronPhonon



subroutine ElectronPhonon_2Fe(La,Peph,Ceph)
  use workplace
  implicit none
  real*8 La
  complex*16, dimension (ndim,ndim,nq) :: Peph,Ceph,Deph

  type (Mformconfig), pointer :: Mi,Mj
  type (Oformconfig), pointer :: Oi,Oj
  integer idim,iOform,ip,npi,jdim,jOform,jp,npj
  real*8 velocity,Veph,w0,w1,fOi,fOj
  real*8 q1Fe(2,nq),eq(nq),Dq(nq),DqQ(nq)

  !notice: this is a temper subroutine for electron phonon coupling for 2Fe/cell. 
  
  Veph=model%Veph; w0=0.002; velocity=0.025; w1=velocity
  
  !rotate q in small BZ by 45 degree (or rotate axes of small BZ by -45 degrees), and divide by rt(2), then q becomes a vector within MBZ of the large BZ
  q1Fe(1,:)=(qv(1,:)-qv(2,:))/2; q1Fe(2,:)=(qv(1,:)+qv(2,:))/2

  eq=w0*w0+velocity*velocity*(4-2*cos(q1Fe(1,:)*pi)-2*cos(q1Fe(2,:)*pi))   !eq=wq*wq
  Dq=0.5*Veph*sqrt( w1*w1/( eq + La*La ) )

  !get q beyond RBZ in the large BZ
  q1Fe(1,:)=q1Fe(1,:)+1; q1Fe(2,:)=q1Fe(2,:)+1

  eq=w0*w0+velocity*velocity*(4-2*cos(q1Fe(1,:)*pi)-2*cos(q1Fe(2,:)*pi))   !eq=wq*wq
  DqQ=0.5*Veph*sqrt( w1*w1/( eq + La*La ) )

  Deph=0

  do jdim=1,ndim; Mj=>Mformfunc(jdim)
  jOform=Mj%Oform; if(Mj%Lform/=1)cycle
  Oj=>Oformfunc(jOform); npj=Oj%npair
  do jp=1,npj; fOj=Oj%f(jp); if(Oj%pair(1,jp)/=Oj%pair(2,jp))cycle

     do idim=1,ndim; Mi=>Mformfunc(idim)
	 iOform=Mi%Oform; if(Mi%Lform/=1)cycle
     Oi=>Oformfunc(iOform); npi=Oi%npair
	 do ip=1,npi; fOi=Oi%f(ip); if(Oi%pair(1,ip)/=Oi%pair(2,ip))cycle	 
	 
     !direct channel
	 if(Mi%atom==Mj%atom)then
       Deph(idim,jdim,:)=Deph(idim,jdim,:)+0.5*(Dq+DqQ)*fOi*fOj
	 else
       Deph(idim,jdim,:)=Deph(idim,jdim,:)+0.5*(Dq-DqQ)*fOi*fOj
     end if

  end do; end do
  end do; end do

  Peph=0; Ceph=0
  if(useprojectors)then; call contact_clever(Peph,Ceph,Deph); else; call contact(Peph,Ceph,Deph); end if

  return
end subroutine ElectronPhonon_2Fe
