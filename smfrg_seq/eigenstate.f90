subroutine basisrotate(t5d)
  implicit none
  complex*16 t5d(5,5,-2:2,-2:2)

  integer ix,iy,iorb,jorb
  complex*16 U(5,5),V(5,5)
  complex*16 one

  !purpose: transform orbits from large unit cell basis to small unit basis

  one=1  
  !one=cmplx(0,1)
  !notice: an 'one' makes xz and yz complex and hk real but blurs time-reversal symmetry
  !caution: the above gauge changes the sign of pair hopping term between xz/yz and other orbits 

  U=0;  U(1,1)=1
  U(2:3,2)=one*(/1,-1/)/sqrt(2.d0)
  U(2:3,3)=one*(/1,1/)/sqrt(2.d0)
  U(5,4)=-1; U(4,5)=1
  V=conjg( transpose(U) )

  do iy=-2,2; do ix=-2,2
     t5d(:,:,ix,iy)=matmul( matmul(V,t5d(:,:,ix,iy)), U )
  end do; end do

  return
end subroutine basisrotate


function findek(nativeband,kv,Q,norb,model)
  use standard_derived_types
  implicit none

  real*8 findek
  integer nativeband,norb
  real*8 kv(2),Q
  type (modelblock) :: model
	type (mfblock) :: mf

  complex*16 hk(norb,norb)
  real*8 eval(norb)

  call gethk(norb,kv,Q,hk,model,mf,.false.)
  call ZHEIGEN(norb,hk,eval) 
  findek=eval(nativeband)

  return
end function findek

subroutine gethk(norb,kv,Q,hk,model,mf,slvbsn)
  use standard_derived_types
  implicit none
  integer norb
  real*8 kv(2),Q
  complex*16 hk(norb,norb)
  type (modelblock) :: model
	type (mfblock) :: mf
	logical slvbsn

  select case (norb)
    case (1); if(model%sample=='sq1a1d')call square1a1dhk(norb,kv,Q,hk,model,mf,slvbsn)
              if(model%sample=='SRO1orb')call SrRuO1orbhk(norb,kv,hk,model,mf,slvbsn)
			  if(model%sample=='BiS2_1b')call BiS2hk_1band(norb,kv,hk,model,mf,slvbsn)
    case (2)
	   if(model%sample=='SRO2orb')call SrRuO2orbhk(norb,kv,hk,model,mf,slvbsn)
	   if(model%sample=='lhq')call lhqhk(norb,kv,hk,model,mf,slvbsn)
	   if(model%sample=='ladder') call openladderhk(norb,kv,Q,hk,model,mf,slvbsn)
	   if(model%sample=='BiS2')call BiS2hk(norb,kv,Q,hk,model,mf,slvbsn)
	   if(model%sample=='Raghu')call Raghuhk(norb,kv,Q,hk,model,mf,slvbsn)
	   if(model%sample=='sq2a2d')call square2a2dhk(kv,hk,model,mf,slvbsn)
	   if(model%sample=='rice')call ricehk(norb,kv,hk,model,mf,slvbsn)
    case (3)
	   if(model%sample=='SRO3orb')call SrRuO3orbhk(norb,kv,hk,model,mf,slvbsn)
	   if(model%sample=='Dagott')call Dagottohk(norb,kv,Q,hk,model,mf,slvbsn)
	   if(model%sample=='CuO2')call CuO2hk(norb,kv,hk,model,mf,slvbsn)	   
       if(model%sample=='kagme')call kagomehk(norb,kv,hk,model,mf,slvbsn)	
	case (4)
	   if(model%sample=='cmpss')call compasshk(norb,kv,hk,model,mf,slvbsn)
	   if(model%sample=='hujp')call hujphk(norb,kv,hk,model,mf,slvbsn)       
	case (5); if(model%sample=='WangFa')call WangFa5dhk(norb,kv,Q,hk,model,mf,slvbsn)
	          if(model%sample=='KFS1L')call KFeSe1L5dhk(norb,kv,Q,hk,model,mf,slvbsn)
	          if(model%sample=='KFS1Lav')call KFeSe1Lavhk(norb,kv,Q,hk,model,mf,slvbsn)
	          if(model%sample=='FeSe1Fe')call FeSe1Fehk(norb,kv,Q,hk,model,mf,slvbsn)
			  if(model%sample=='Kuroki')call Kuroki5dhk(norb,kv,Q,hk,model,mf,slvbsn)
			  if(model%sample=='FS1F')call FeSe1Fehk(norb,kv,Q,hk,model,mf,slvbsn)
			  if(model%sample=='FS1FLee')call FS1FLeehk(norb,kv,Q,hk,model,mf,slvbsn)
			  if(model%sample=='LFA1L5d')call LiFeAs1L5dhk(norb,kv,Q,hk,model,mf,slvbsn)
			  if(model%sample=='LFA2L5d')call LiFeAs2L5dhk(norb,kv,Q,hk,model,mf,slvbsn)
			  if(model%sample=='BFA1L5d')call BaFeAs1L5dhk(norb,kv,Q,hk,model,mf,slvbsn)
	case (8)  
	   if(model%sample=='hujp2')call hujp2hk(norb,kv,hk,model,mf,slvbsn)
    case (10)
	  if(model%sample=='LFA1L10d')call LiFeAs1L10dhk(norb,kv,Q,hk,model,mf,slvbsn)
	  if(model%sample=='KFS2vL')call KFeSe2vLhk(norb,kv,hk,model,mf,slvbsn)
	  if(model%sample=='K1L10d')call KFeSe1L10dhk(norb,kv,Q,hk,model,mf,slvbsn)
	  if(model%sample=='KFSPrm')call KFeSePrmhk(norb,kv,Q,hk,model,mf,slvbsn)
      if(model%sample=='FS2F')call FeSe2Fehk(norb,kv,Q,hk,model,mf,slvbsn)
	  if(model%sample=='BFA1L10d')call BFA1L10dhk(norb,kv,Q,hk,model,mf,slvbsn)
    case (20); call KFeSe2L20dhk(norb,kv,hk,model,mf,slvbsn)
	case default; stop 'model hk unknown'
  end select
  return
end subroutine gethk






