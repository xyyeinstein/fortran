subroutine bubble(La,vq,Xpp,Xph)
  use workplace
  implicit none

  real*8 La,vq(2)
  complex*16, dimension (ndim,ndim) :: Xpp,Xph

  integer ik   
  real*8 vk(2),sk
  real*8 fOi,fOj
  complex*16 fki,fkj,add,factor,formk
  complex*16, dimension (norb,norb) :: Gw,G_w,Gqw,Gq_w

  integer idim,iOform,iLform,npi,ip
  integer iorb1,iorb2,jorb1,jorb2,i,j

  type (Mformconfig), pointer :: Mi
  type (Oformconfig), pointer :: Oi
  complex*16, dimension (norb*norb,norb*norb) :: X,Y
  complex*16 fdag(ndim,norb*norb),f(norb*norb,ndim)
  integer nnorb

  nnorb=norb*norb

  Xpp=0; Xph=0
  
  do ik=1,nk; vk=kv(1:2,ik); sk=kv(3,ik)
    
     call green(La,vk,Gw,G_w)
     call green(La,vk+vq,Gqw,Gq_w)

     fdag=0
     do idim=1,ndim; Mi=>Mformfunc(idim)	 
     iOform=Mi%Oform; iLform=Mi%Lform; fki=formk(-vk,iLform)
	 if(QuickSearch.and.iLform/=1)cycle
     Oi=>Oformfunc(iOform); npi=Oi%npair
     do ip=1,npi; fOi=Oi%f(ip); iorb1=Oi%pair(1,ip); iorb2=Oi%pair(2,ip)
     if(natom>1)call modifyindex(iorb1,Mi%atom,iorb2,Mi%nbatom)
	    i=iorb1+(iorb2-1)*norb
		fdag(idim,i)=fdag(idim,i)+fOi*fki
	 end do; end do

	 f=conjg(transpose(fdag))

     do jorb1=1,norb; do jorb2=1,norb; j=jorb1+(jorb2-1)*norb
     do iorb1=1,norb; do iorb2=1,norb; i=iorb1+(iorb2-1)*norb
        X(i,j) =  Gqw(iorb1,jorb1)*conjg(Gw(iorb2,jorb2)) + Gq_w(iorb1,jorb1)*conjg(G_w(iorb2,jorb2)) 	    
        if(BCSflow)cycle
	    Y(i,j) =  Gqw(iorb1,jorb1)*Gw(jorb2,iorb2) + Gq_w(iorb1,jorb1)*G_w(jorb2,iorb2) 
	 end do; end do; end do; end do
	 
	 X = X*sk;	 Xpp = Xpp + matmul( fdag, matmul( X, f ) )
	 if(BCSflow)cycle
	  
	 Y = Y*sk;   Xph = Xph + matmul( fdag, matmul( Y, f ) )
  
  end do
  return
end subroutine bubble
























subroutine bubble_old(La,vq,Xpp,Xph)
  use workplace
  implicit none

  real*8 La,vq(2)
  complex*16, dimension (ndim,ndim) :: Xpp,Xph

  integer ik   
  real*8 vk(2),sk
  real*8 fOi,fOj
  complex*16 fki,fkj,add,factor,formk
  complex*16, dimension (norb,norb) :: Gw,G_w,Gqw,Gq_w

  integer idim,iOform,iLform,npi,ip
  integer jdim,jOform,jLform,npj,jp
  integer iorb1,iorb2,jorb1,jorb2

  type (Mformconfig), pointer :: Mi,Mj
  type (Oformconfig), pointer :: Oi,Oj

  Xpp=0; Xph=0
  
  do ik=1,nk; vk=kv(1:2,ik); sk=kv(3,ik)
    
     call green(La,vk,Gw,G_w)
     call green(La,vk+vq,Gqw,Gq_w)

     do jdim=1,ndim; Mj=>Mformfunc(jdim)
     jOform=Mj%Oform; jLform=Mj%Lform; fkj=formk(vk,jLform)
	 if(QuickSearch.and.jLform/=1)cycle
     Oj=>Oformfunc(jOform); npj=Oj%npair
     do jp=1,npj; fOj=Oj%f(jp); jorb1=Oj%pair(1,jp); jorb2=Oj%pair(2,jp)
     if(natom>1)call modifyindex(jorb1,Mj%atom,jorb2,Mj%nbatom)

     do idim=1,ndim; Mi=>Mformfunc(idim)	 
     iOform=Mi%Oform; iLform=Mi%Lform; fki=formk(-vk,iLform)
	 if(QuickSearch.and.iLform/=1)cycle
     Oi=>Oformfunc(iOform); npi=Oi%npair
     do ip=1,npi; fOi=Oi%f(ip); iorb1=Oi%pair(1,ip); iorb2=Oi%pair(2,ip)
     if(natom>1)call modifyindex(iorb1,Mi%atom,iorb2,Mi%nbatom)

        factor=fki*fkj*fOi*fOj*sk

        add=  Gqw(iorb1,jorb1)*conjg(Gw(iorb2,jorb2)) + Gq_w(iorb1,jorb1)*conjg(G_w(iorb2,jorb2)) 
	    Xpp(idim,jdim) = Xpp(idim,jdim) + add*factor

        if(BCSflow)cycle
	    add =  Gqw(iorb1,jorb1)*Gw(jorb2,iorb2) + Gq_w(iorb1,jorb1)*G_w(jorb2,iorb2) 
	    Xph(idim,jdim) = Xph(idim,jdim) + add*factor

	 end do; end do; end do; end do
  
  end do
  return
end subroutine bubble_old


