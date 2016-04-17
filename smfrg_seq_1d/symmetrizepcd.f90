subroutine symmetrizePCD_R(iR,iatom,jatom,PR,CR,DR)
  use workplace
  implicit none
  integer iR,iatom,jatom
  complex*16, dimension (ndim,ndim,nR) :: PR,CR,DR

  !purpose: use V(:,:,iR) with iatom,jatom as the source to get V(:,:,iRg) by symmetry, where iRg is the 
  !         image of iR under group g. 
  !notice:  the way how symmetry operation is implemented is described @ symmetrizePCD_q

  integer ig,iRg,iag,jag,idim,jdim,i,j,iidim,jjdim,igdim(8),jgdim(8)
  real*8 phasei(8),phasej(8),phase
  integer info(3),ni,nj
  integer nonzero  

  do ig=2,ng; info=indexRg(:,ig,iR,iatom,jatom)
     
	 iRg=info(1); iag=info(2); jag=info(3)
	 if(iRg==iR.and.iag==iatom.and.jag==jatom)cycle
	 if(iRg<=0.or.iag<=0.or.jag<=0)cycle
     
	 !nullifiy relevant V(:,:,iRg). This is necessary since it might have been assigned by previous symmetry operations 
	 do jdim=1,ndim; if(Mformfunc(jdim)%atom/=jag)cycle
	 do idim=1,ndim; if(Mformfunc(idim)%atom/=iag)cycle
	    PR(idim,jdim,iRg)=0; CR(idim,jdim,iRg)=0; DR(idim,jdim,iRg)=0
	 end do; end do

     !get relevant V(:,:,iRg) from V(:,:,iR) by symmetry
	 do jdim=1,ndim; if(Mformfunc(jdim)%atom/=jatom)cycle; jgdim=indexgroup(:,ig,jdim); nj=nonzero(8,jgdim); phasej=Xg(:,ig,jdim)
	 do idim=1,ndim; if(Mformfunc(idim)%atom/=iatom)cycle; igdim=indexgroup(:,ig,idim); ni=nonzero(8,igdim); phasei=Xg(:,ig,idim)
        do j=1,nj; if(jgdim(j)<=0)cycle; do i=1,ni; if(igdim(i)<=0)cycle
	       iidim=igdim(i); if(Mformfunc(iidim)%atom/=iag)cycle
		   jjdim=jgdim(j); if(Mformfunc(jjdim)%atom/=jag)cycle
	       phase=phasei(i)*phasej(j)           
		   PR(iidim,jjdim,iRg)=PR(iidim,jjdim,iRg)+PR(idim,jdim,iR)*phase
		   CR(iidim,jjdim,iRg)=CR(iidim,jjdim,iRg)+CR(idim,jdim,iR)*phase
		   DR(iidim,jjdim,iRg)=DR(iidim,jjdim,iRg)+DR(idim,jdim,iR)*phase
 	    end do; end do
	 end do; end do
  end do
  return
end subroutine symmetrizePCD_R

subroutine symmetrizePCD_q(iq,Pq,Cq,Dq)
  use workplace
  implicit none
  integer iq
  complex*16, dimension (ndim,ndim,nq) :: Pq,Cq,Dq

  !purpose: use V(:,:,iq) as the source to get V(:,:,iqg) by symmetry, where iqg is the image of iq under group g. 
  !notice:
  !
  ! Under group operation g (and denote its inverse as g'), the form function transform as
  !
  !             g [f_m(k)] = f_m(g'k) = sum_n Xg(m,n) f_n (k),  g' = inv(g). 
  !
  ! In general, group symmetry requires
  !
  !             g [ sum_mn f_m (k) P_mn(q) f_n(k') ] : = sum_mn f_m(g'k) P_mn(g'q) f_n(g'k')
  ! 
  ! to be an invariant under the point group. Here k should be understood as all internal variables appearing in a 
  ! form factor. This symmetry leads to the condition, in matrix form,
  !
  !                     P(q) = Xg^T P(g'q) Xg,  g'=inv(g),  for real Xg
  !
  ! where Xg follows from the transformation law for the form factors. Now replace q by gq in the above to find
  !
  !                     P(gq) = Xg^T P(q) Xg,
  !
  ! or in components,
  !
  !                   P(m',n',qg) = sum_(m,n) Xg(m,m') P(m,n,q) Xg(n,n').
  !
  ! This is the general transformation law for interactions P, C and D, even if the form factors are not the natural
  ! irreducible representations of the point group. The only requirement is that they form an orthonormal complete set. 		
  ! For one-to-one transformation matrices (mentioned above) there is only one term in the summation for given m'n'. 
  !
  ! The fundamental form factors, i.e., each of which includes only one bond, one orbital pair, and one aotm pair, all 
  ! transform in the one-to-one fashion. So the symmetry implementation for such form factors is very straightforward.
  ! On hexagonal lattices, the transformation is at most one-to-two, and is implemented in the present package.

  integer ig,iqg,idim,jdim,i,j,iidim,jjdim,igdim(8),jgdim(8),ni,nj
  real*8 phasei(8),phasej(8),phase
  integer nonzero

  do ig=2,ng; iqg=indexqg(ig,iq); 
		if(iqg==iq)cycle

	 Pq(:,:,iqg)=0; Cq(:,:,iqg)=0; Dq(:,:,iqg)=0       
	 do jdim=1,ndim; jgdim=indexgroup(:,ig,jdim); nj=nonzero(8,jgdim); phasej=Xg(:,ig,jdim)
	 do idim=1,ndim; igdim=indexgroup(:,ig,idim); ni=nonzero(8,igdim); phasei=Xg(:,ig,idim)
        do j=1,nj; if(jgdim(j)<=0)cycle; do i=1,ni; if(igdim(i)<=0)cycle
	       iidim=igdim(i); jjdim=jgdim(j)
	       phase=phasei(i)*phasej(j)           
		   Pq(iidim,jjdim,iqg)=Pq(iidim,jjdim,iqg)+Pq(idim,jdim,iq)*phase
		   if(BCSflow)cycle
		   Cq(iidim,jjdim,iqg)=Cq(iidim,jjdim,iqg)+Cq(idim,jdim,iq)*phase
		   Dq(iidim,jjdim,iqg)=Dq(iidim,jjdim,iqg)+Dq(idim,jdim,iq)*phase
 	    end do; end do
	end do; end do
  end do
  return
end subroutine symmetrizePCD_q

function nonzero(n,info)
  implicit none
  integer nonzero,n
  integer info(n)

  integer i

  nonzero=0; do i=1,n; if(info(i)>0)nonzero=nonzero+1; end do
  
  return
end function nonzero 
