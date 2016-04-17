subroutine initV()
  use workplace
  implicit none

  integer idim,jdim
  
  type (Mformconfig), pointer :: Mi,Mj
	 
  allocate(P(ndim,ndim,nq),C(ndim,ndim,nq),D(ndim,ndim,nq))

  P=0; C=0; D=0

  do jdim=1,ndim; Mj=>Mformfunc(jdim); if(Mj%Lform/=1)cycle
  do idim=1,ndim; Mi=>Mformfunc(idim); if(Mi%Lform/=1)cycle; if(Mi%atom/=Mj%atom)cycle
	 if(Mi%spin==Mj%spin.and.Mi%nbspin==Mj%nbspin.and.Mi%spin/=Mi%nbspin)P(idim,jdim,:)=-model%U
	 if(Mi%spin==Mj%spin.and.Mi%nbspin==Mj%nbspin.and.Mi%spin/=Mi%nbspin)C(idim,jdim,:)=-model%U
	 if(Mi%spin==Mi%nbspin.and.Mj%spin==Mj%nbspin.and.Mi%spin/=Mj%spin)D(idim,jdim,:)=-model%U
  end do; end do
  
  if(abs(model%Vnn)>1.e-5)call addVnn()
  if(abs(model%VL2L)>1.e-5)call addVL2L()
  if(abs(model%Vbonding)>1.e-5)call addVbonding()
  if(abs(model%Vantibonding)>1.e-5)call addVantibonding()
  
  call init_antisymmetry()
  return
end subroutine initV

subroutine addVnn()
  use workplace
  implicit none

  integer idim,jdim
  
  type (Mformconfig), pointer :: Mi,Mj
	 
  complex*16, dimension (ndim,ndim,nq) :: dP,dC,dD
  real*8 q(2),bond(2),rv(3)
  integer i,iq
  complex*16 add
  integer findatom

  dP=0; dC=0; dD=0
  bond=(/1,0/)

  do i=1,4; if(i>1)call rotate90(bond)
  do iq=1,nq; q=qv(1:2,iq)
     add=-model%Vnn*exp(one*pi*sum(q*bond))
     do jdim=1,ndim; Mj=>Mformfunc(jdim); if(Mj%Lform/=1)cycle
     do idim=1,ndim; Mi=>Mformfunc(idim); if(Mi%Lform/=1)cycle
	    if(Mi%spin/=Mi%nbspin.or.Mj%spin/=Mj%nbspin)cycle
        rv=model%ra(:,Mi%atom); rv(1:2)=rv(1:2)+bond	    
		if(findatom(rv,model)/=Mj%atom)cycle		   
		dD(idim,jdim,iq)=dD(idim,jdim,iq)+add
	 end do; end do
  end do; end do

  call contact(dP,dC,dD)
  P=P+dP; C=C+dC; D=D+dD

  return
end subroutine addVnn


subroutine addVL2L()
  use workplace
  implicit none

  integer idim,jdim
  
  type (Mformconfig), pointer :: Mi,Mj
	 
  complex*16, dimension (ndim,ndim,nq) :: dP,dC,dD
  real*8 q(2),bond(2),rv(3)
  integer i,iq
  complex*16 add
  integer findatom

  dP=0; dC=0; dD=0
  bond=(/1,0/)

  do jdim=1,ndim; Mj=>Mformfunc(jdim); if(Mj%Lform/=1)cycle
  do idim=1,ndim; Mi=>Mformfunc(idim); if(Mi%Lform/=1)cycle
     if(Mi%atom==Mj%atom)cycle
	 if(Mi%spin/=Mi%nbspin.or.Mj%spin/=Mj%nbspin)cycle
	 dD(idim,jdim,:)=dD(idim,jdim,:)-model%VL2L
  end do; end do

  call contact(dP,dC,dD)
  P=P+dP; C=C+dC; D=D+dD

  return
end subroutine addVL2L

subroutine addVbonding()
  use workplace
  implicit none

  integer idim,jdim
  
  type (Mformconfig), pointer :: Mi,Mj
	 
  complex*16, dimension (ndim,ndim,nq) :: dP,dC,dD
  real*8 q(2),bond(2),rv(3)
  integer i,iq
  complex*16 add
  integer findatom

  dP=0; dC=0; dD=0

  do jdim=1,ndim; Mj=>Mformfunc(jdim); if(sum(abs(Mj%Lformfunc%r(1:2,1)))>1.e-5)cycle; if(Mj%spin/=Mj%nbspin)cycle
  do idim=1,ndim; Mi=>Mformfunc(idim); if(sum(abs(Mi%Lformfunc%r(1:2,1)))>1.e-5)cycle; if(Mi%spin/=Mi%nbspin)cycle
	 dD(idim,jdim,:)=dD(idim,jdim,:)-model%VBonding/4
  end do; end do

  call contact(dP,dC,dD)
  P=P+dP; C=C+dC; D=D+dD

  return
end subroutine addVbonding


subroutine addVantibonding()
  use workplace
  implicit none

  integer idim,jdim
  
  type (Mformconfig), pointer :: Mi,Mj
	 
  complex*16, dimension (ndim,ndim,nq) :: dP,dC,dD
  real*8 q(2),bond(2),rv(3)
  integer i,iq,isign
  complex*16 add
  integer findatom

  dP=0; dC=0; dD=0

  do jdim=1,ndim; Mj=>Mformfunc(jdim); if(sum(abs(Mj%Lformfunc%r(1:2,1)))>1.e-5)cycle; if(Mj%spin/=Mj%nbspin)cycle
  do idim=1,ndim; Mi=>Mformfunc(idim); if(sum(abs(Mi%Lformfunc%r(1:2,1)))>1.e-5)cycle; if(Mi%spin/=Mi%nbspin)cycle
     isign=1
	 if(Mi%atom/=Mi%nbatom)isign=-isign
	 if(Mj%atom/=Mj%nbatom)isign=-isign
	 dD(idim,jdim,:)=dD(idim,jdim,:)-isign*model%VantiBonding/4
  end do; end do

  call contact(dP,dC,dD)
  P=P+dP; C=C+dC; D=D+dD

  return
end subroutine addVantibonding


subroutine addJnn()
  use workplace
  implicit none

  integer idim,jdim
  
  type (Mformconfig), pointer :: Mi,Mj
	 
  complex*16, dimension (ndim,ndim,nq) :: dP,dC,dD
  real*8 q(2)
  integer iq
  complex*16 add

  dP=0; dC=0; dD=0
  
  do iq=1,nq; q=qv(1:2,iq)
     add=-model%Jnn*2*sum( cos(pi*q) )
     do jdim=1,ndim; Mj=>Mformfunc(jdim); if(Mj%Lform/=1)cycle
     do idim=1,ndim; Mi=>Mformfunc(idim); if(Mi%Lform/=1)cycle
	    if(Mi%atom/=Mj%atom)cycle
		if(Mi%spin/=Mj%spin)cycle
		if(Mi%nbspin/=Mj%nbspin)cycle
		dC(idim,jdim,iq)=dC(idim,jdim,iq)+add
	 end do; end do
  end do

  call contact(dP,dC,dD)
  P=P+dP; C=C+dC; D=D+dD

  return
end subroutine addJnn


subroutine init_antisymmetry()
  use workplace
  implicit none

  complex*16, dimension (ndim,ndim,nq) :: dP,dC,dD  
 
  dP=0; dD=-(C-D)/2; dC=-dD
  call contact(dP,dC,dD)  
!  if(sum(abs(dC-(C-D)))>1.e-4)stop 'dC error @ init_antisymmetry'
!  if(sum(abs(dD+(C-D)))>1.e-4)stop 'dD error @ init_antisymmetry'
  P=dP; C=dC; D=dD
!  if(sum(abs(C+D))>1.e-4)stop 'init_antisymmetry for C and D failed'
  

  dP=P
  call PauliExclusion()
  if(sum(abs(dP-P))>1.e-4)stop 'initial antisymmetry for P failed'


  return
end subroutine init_antisymmetry















subroutine init_antisymmetry_old()
  use workplace
  implicit none

  integer idim,jdim,isign
  
  type (Mformconfig), pointer :: Mi,Mj
	 
  C=(C-D)/2; D=-C

  if(natom==1)then
   
  do idim=1,ndim; Mi=>Mformfunc(idim)
     call partner(idim,isign,jdim)
     P(idim,:,:)=( P(idim,:,:)-isign*P(jdim,:,:) ) / 2
     P(jdim,:,:)=-isign*P(idim,:,:)
  end do

  do idim=1,ndim; Mi=>Mformfunc(idim)
     call partner(idim,isign,jdim)
     P(:,idim,:)=( P(:,idim,:)-isign*P(:,jdim,:) ) / 2
     P(:,jdim,:)=-isign*P(:,idim,:)
  end do
  
  else

  if(model%Vnn/=0.or.model%VL2L/=0)stop 'PauliExclusion not sufficient here'
  call PauliExclusion()   !@singularflow.f90

  end if

  P=2*P; C=2*C; D=2*D

  return

  contains
    subroutine partner(idim,isign,jdim)
	  integer idim,isign,jdim
	  integer i,icount
	  real*8 ri(3)
	  integer findatom
	  icount=0; isign=1

	  do i=1,ndim; Mj=>Mformfunc(i)
	     if(Mj%atom/=Mi%nbatom)cycle
	     if(Mj%form/=Mi%form)cycle
		 if(Mj%spin/=Mi%nbspin)cycle
		 if(Mj%nbspin/=Mi%spin)cycle
		 ri=model%ra(:,Mi%atom)+Mi%Lformfunc%r(:,1)
		 if(findatom(ri,model)/=Mj%atom)cycle
         jdim=i; icount=icount+1
	  end do
	  if(icount==0)stop 'partner not found @ antisymmetry'
      if(icount>1)stop 'too many partner found @ antisymmetry'
	  if(Mi%form=='E1u'.or.Mi%form=='E2u')isign=-1
	  return
	end subroutine partner
end subroutine init_antisymmetry_old
