subroutine eigentable()
  use workplace
  implicit none

  integer ix,iy
  real*8 vk(2)
  complex*16 hk(norb,norb)
  real*8 eval(norb)

  open(10,file='ek.dat')
  do iy=-nky,nky; do ix=-nkx,nkx; vk=(/ix,iy/)*1.4/(/nkx,nky/)
     call gethk(norb,vk,hk,model)
	 call ZHEIGEN(norb,hk,eval)
	 write(10,100)vk,eval
  end do; end do
  close(10)
100 format(1x,5f15.6)
  
  return
end subroutine eigentable

