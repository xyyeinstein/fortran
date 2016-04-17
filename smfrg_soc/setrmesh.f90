subroutine setRmesh()
  use workplace
  implicit none

  integer iR,ix,iy

  nR=2*Lcontact+1; nR=nR*nR
  allocate(R(2,nR))
  
  iR=0
  do iy=-Lcontact,Lcontact; do ix=-Lcontact,Lcontact
     iR=iR+1
	 R(:,iR)=(/ix,iy/)
  end do; end do
  
  return
end subroutine setRmesh
