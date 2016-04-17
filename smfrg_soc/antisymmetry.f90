subroutine antisymmetrize()
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
end subroutine antisymmetrize
