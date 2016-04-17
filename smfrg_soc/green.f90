subroutine Green(w,vk,Gw,G_w)
  use workplace
  implicit none
  real*8 w,vk(2)
  complex*16, dimension (norb,norb) :: Gw,G_w

  real*8 eval(norb)
  complex*16 hk(norb,norb)
  integer i

  call gethk(norb,vk,hk,model)

  Gw=-hk; do i=1,norb; Gw(i,i)=Gw(i,i)+one*w; end do
  call ZINVERT(norb,Gw)

  G_w=-hk; do i=1,norb; G_w(i,i)=G_w(i,i)-one*w; end do
  call ZINVERT(norb,G_w)

  return
end subroutine Green