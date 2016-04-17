subroutine bubble(La,vq,Xpp,Xph)
  use workplace
  implicit none

  real*8 La,vq(2)
  complex*16, dimension (ndim,ndim) :: Xpp,Xph

  integer ik   
  real*8 vk(2),sk
  real*8 fOi,fOj
  complex*16 fki,fkj,add,factor,formk
  complex*16, dimension (norb,norb) :: Gw,G_w,Gqw,Gq_w,G_kw,G_k_w

  integer idim,iOform,iLform,npi,ip
  integer jdim,jOform,jLform,npj,jp
  integer iorb1,iorb2,jorb1,jorb2

  type (Mformconfig), pointer :: Mi,Mj
  type (Lformconfig), pointer :: Li,Lj

  Xpp=0; Xph=0
  
  do ik=1,nk; vk=kv(1:2,ik); sk=kv(3,ik)
    
     call green(La,vk,Gw,G_w);  call green(La,-vk,G_kw,G_k_w)
     call green(La,vk+vq,Gqw,Gq_w)

     do jdim=1,ndim; Mj=>Mformfunc(jdim)
     jLform=Mj%Lform; fkj=formk(vk,Mj%Lformfunc)
     jorb1=Mj%orb; jorb2=Mj%nborb
	 if(QuickSearch.and.jLform/=1)cycle

     do idim=1,ndim; Mi=>Mformfunc(idim)	 
     iLform=Mi%Lform; fki=formk(vk,Mi%Lformfunc)
     iorb1=Mi%orb; iorb2=Mi%nborb
	 if(QuickSearch.and.iLform/=1)cycle

        factor=conjg(fki)*fkj*sk

        add=  Gqw(iorb1,jorb1)*G_k_w(iorb2,jorb2) + Gq_w(iorb1,jorb1)*G_kw(iorb2,jorb2) 
	    Xpp(idim,jdim) = Xpp(idim,jdim) + add*factor

	    add =  Gqw(iorb1,jorb1)*Gw(jorb2,iorb2) + Gq_w(iorb1,jorb1)*G_w(jorb2,iorb2) 
	    Xph(idim,jdim) = Xph(idim,jdim) + add*factor
 
	 end do; end do
  
  end do
  return
end subroutine bubble


