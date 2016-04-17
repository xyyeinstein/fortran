subroutine normalhk(norb,vk,hk,model)
  use standard_derived_types
  implicit none

  integer norb
  real*8 vk(2)
  complex*16 hk(norb,norb)
  type (modelblock) :: model

  !purpose: get the normal state hamiltonian in the momentum space

  !input arguments:
  !  norb: number of effective orbits (including spin)
  !  vk :  momentum vector in units of pi
  !  sample: specifying the model

  !output arguments:
  !  hk:   normal state hamiltonian

  select case (model%sample)
    case ('cuprt' ); call cupratehk(norb,vk,hk,model)
	case default; stop 'hk model undefined'
  end select

  return
end subroutine normalhk

subroutine pairinggk(norb,vk,gk,model)
  use standard_derived_types
  implicit none

  integer norb
  real*8 vk(2)
  complex*16 gk(norb,norb)
  type (modelblock) :: model

  !purpose: get the pairing gap in the momentum space
  
  !input arguments:
  !  norb: number of effective orbits (including spin)
  !  vk:   momentum in units of pi
  !  sample: specifying the model

  !output arguments:
  !  gk:  gap function in the momentum space

  select case (model%sample)
    case ('cuprt' ); call cuprategk(norb,vk,gk)
    !case ('sqr2L');  call cuprate2Lgk(norb,vk,gk)
	!case ( 'graphene'); call graphenegk(norb,vk,gk,model)
	case default; stop 'gk model undefined'
  end select

  return
end subroutine pairinggk

subroutine getBdGAk(norb,vk,nambu,Ak,model,edgemodel)
  use standard_derived_types
  implicit none

  integer norb,nambu
  real*8 vk(2)
  complex*16 Ak(nambu,nambu),hk(norb,norb),gk(norb,norb)
  type (modelblock) :: model
  type (edgemodelblock) :: edgemodel

  !purpose: get the momentum space hamiltonian in the nambu space
  
  !input arguments
  !  norb: number of effective orbits in the model (including spin)
  !  vk:   k vector in units of pi
  !  nambu: dimension of nambu space
  !  sample: specifying the model

  !output argument:
  !  ak:  momentum space hamiltonian in the nambu space

  !notice:
  !  if nambu=norb, the program deals with normal state hamiltonian
  !  if nambu=2*norb, the program deals with superconducting state


  call normalhk(norb,vk,hk,model); Ak(1:norb,1:norb)=hk
  
  if(norb==nambu)return    !return for normal state

  call normalhk(norb,-vk,hk,model); Ak(norb+1:nambu,norb+1:nambu)=-conjg(hk)	    
  call pairinggk(norb,vk,gk,model)
  Ak(1:norb,norb+1:nambu)=gk; Ak(norb+1:nambu,1:norb)=conjg(transpose(gk))

  return
end subroutine getBdgAk


subroutine gethr(nambu,hq,Lb,ndim,h,periodicslab)
  implicit none

  integer nambu,Lb,ndim
  logical periodicslab
  complex*16 hq(nambu,nambu,-4:4),h(ndim,ndim)

  !purpose: construct nambu-space hamiltonian for an open/periodic boundary slice with conserved transverse momentum

  !input arguments:
  !     nambu : dimension of nambu space
  !     hq(:,:,b): effective hopping matrix on bond b, with conserved transverse momentum
  !     Lb: length of the open slice
  !     ndim = nambu*Lb
  !     periodicslab: boundary condition 
  
  !output arguments:
  !     h(:,:) : effective nambu space hamiltonian for the open-boundary slice

  integer i,j,m,inambu,jnambu,idim,jdim

  h=0
  do i=1,Lb; do m=-4,4; j=i-m
	 if(periodicslab)then
	   j=mod(j-1+2*Lb,Lb)+1               !periodic boundary
	 else
	   if(j<1)cycle; if(j>Lb)cycle        !open boundary
	 end if
	 do jnambu=1,nambu; jdim=jnambu+(j-1)*nambu
	 do inambu=1,nambu; idim=inambu+(i-1)*nambu
	    h(idim,jdim)=h(idim,jdim)+hq(inambu,jnambu,m)
	 end do; end do
  end do; end do
  return 
end subroutine gethr

subroutine gethq(nb,q,norb,nambu,hq,model,edgemodel)
  use standard_derived_types
  implicit none

  integer nb,norb,nambu
  real*8 q(2),bond(2)
  complex*16 hq(nambu,nambu,-4:4)
  type (modelblock) :: model
  type (edgemodelblock) :: edgemodel
  !purpose: get the effective hopping matrix along the nb-th direction with a conserved transverse momentum q

  !input arguments:
  !  nb: normal direction of the boundary
  !  q:  conserved transverse momentum
  !  norb: number of effective orbits in the normal state (including spin)
  !  nambu: dimension of nambu space
  !  sample: specifying the model
  
  !output arguments:
  !  hq(:,:,b) : effective hopping matrix along all bonds b with the conserved transverse momentum q
  
  !notice:
  !  if nambu=2*norb, the program deals with a superconducting sample
  !  if nambu=norb,   the program deals with a normal-state sample 

  !notice:
  !  it is implicitly assumed that the bond length for hq is at most 2

  real*8 vk(2),pi
  integer ik,Nk,i
  complex*16 Ak(nambu,nambu),one

  Nk=8; pi=asin(1.d0)*2; one=cmplx(0,1)

  hq=0
  do ik=0,Nk-1
     if(nb==1)then
       vk=ik*model%ka(1:2)/Nk;  bond=model%a(1:2)
	 else
	   vk=ik*model%kb(1:2)/nk;  bond=model%b(1:2)
	 end if
	 vk=vk+q
     call getBdGAk(norb,vk,nambu,Ak,model,edgemodel)
	 do i=-4,4
 	    hq(:,:,i)=hq(:,:,i)+Ak*exp(one*i*sum(bond*vk)*pi)
     end do
  end do
  hq=hq/Nk

  return
end subroutine gethq