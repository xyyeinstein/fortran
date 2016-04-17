subroutine singularflow(La,dLa)
  use workplace
  implicit none
  
  real*8 La,dLa
  
  integer iq
  real*8 q(2)

  complex*16, dimension (ndim,ndim,nq) :: dP,dC,dD
  complex*16, dimension (ndim,ndim) :: Xpp,Xph

  integer ig,iqg
  logical reducible
  integer idim,jdim,igdim,jgdim
  complex*16 phasePP,phasePH

  dP=0; dC=0; dD=0

  do iq=1,nq; q=qv(1:2,iq)
  
     if(usegroup.and.reducible(q,group))cycle

	 call bubble(La,q,Xpp,Xph)
	 Xpp=Xpp*dLa/twopi; Xph=Xph*dLa/twopi
	 
	 dP(:,:,iq)=0.5*matmul( P(:,:,iq), matmul(Xpp,P(:,:,iq)) )
	 dC(:,:,iq)=matmul( C(:,:,iq), matmul(Xph,C(:,:,iq)) )
	 !dD(:,:,iq)=-matmul( D(:,:,iq), matmul(Xph,D(:,:,iq)) ) 
	 !dD(:,:,iq)=matmul( C(:,:,iq)-D(:,:,iq), matmul(Xph,D(:,:,iq)) ) &
	 !         + matmul( D(:,:,iq), matmul(Xph,C(:,:,iq)-D(:,:,iq)) )
	 
	 if( .not. usegroup )cycle
     do ig=2,ng; iqg=indexqg(ig,iq); if(iqg==iq)cycle
	    dP(:,:,iqg)=0; dC(:,:,iqg)=0; dD(:,:,iqg)=0       
		do jdim=1,ndim; jgdim=indexgroup(ig,jdim)
		do idim=1,ndim; igdim=indexgroup(ig,idim)
 	       phasePP=ppXg(ig,idim)*conjg(ppXg(ig,jdim))
		   phasePH=phXg(ig,idim)*conjg(phXg(ig,jdim))
			  			  		      
			  !Notice of symmetry operation:
			  !    Vg = X(g) V X(g)^T,          X(g)^T = X(-g).
			  !    Vg_(m,n) = sum_(m',n') X(g)_(m,m') V_(m',n') X(g)_(n,n')

		   dP(idim,jdim,iqg)=dP(idim,jdim,iqg)+dP(igdim,jgdim,iq)*phasePP
		   dC(idim,jdim,iqg)=dC(idim,jdim,iqg)+dC(igdim,jgdim,iq)*phasePH
	      !dD(idim,jdim,iqg)=dD(idim,jdim,iqg)+dD(iidim,jjdim,iq)*phasePH
		end do; end do
	 end do

  end do

  dD=-dC
  if(useprojectors)then; call contact_clever(dP,dC,dD); else; call contact(dP,dC,dD); end if 
  
  P=P+dP; C=C+dC; D=D+dD
  call PauliExclusion()
  return
end subroutine singularflow

subroutine PauliExclusion()
  use workplace
  implicit none

  type (Mformconfig), pointer :: Mi,Mj
  integer idim,jdim

  !antisymmetrize particle-hole interactions
  C=(C-D)/2; D=-C

  !antisymmetrize pairing interactions, for on-site forms only
  do idim=1,ndim; Mi=>Mformfunc(idim); if(Mi%Lform/=1)cycle
  do jdim=1,ndim; Mj=>Mformfunc(jdim); if(Mj%Lform/=1)cycle
     if(Mi%atom/=Mj%atom)cycle
     if(Mi%spin/=Mj%nbspin)cycle
	 if(Mi%nbspin/=Mj%spin)cycle
     P(idim,:,:)=( P(idim,:,:)-P(jdim,:,:) )/2
     P(jdim,:,:)=-P(idim,:,:)
  end do; end do

  do idim=1,ndim; Mi=>Mformfunc(idim); if(Mi%Lform/=1)cycle
  do jdim=1,ndim; Mj=>Mformfunc(jdim); if(Mj%Lform/=1)cycle
     if(Mi%atom/=Mj%atom)cycle
     if(Mi%spin/=Mj%nbspin)cycle
	 if(Mi%nbspin/=Mj%spin)cycle
     P(:,idim,:)=( P(:,idim,:)-P(:,jdim,:) )/2
     P(:,jdim,:)=-P(:,idim,:)
  end do; end do

  return
end subroutine PauliExclusion
