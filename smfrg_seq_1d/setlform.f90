
subroutine setLformfunc(ngmax)
  use workplace
  use standard_derived_types
  implicit none
  integer ngmax

  character*3 form

  integer iform,jform,ig
  real*8 basis(3)
  
  if(ngmax==1)print*,'setting fundamental lattice form factors since ngmax=1 @ setLformfunc'
  allocate(Lformfunc(nLform))
  
  do iform=1,nLform
     form=Lformtable(iform)
     basis=LformBasis(:,iform)
     call formstructure(3,basis,form,Lformfunc(iform),ngmax,group)
  end do

  return
end subroutine setLformfunc



function formk(vk,iform)
  use workplace, only : pi,one,Lformfunc
  implicit none
  real*8 vk(2)
  complex*16 formk
  integer iform

  !notice: assuming the third dimension of vk, if any, is zero. 

  integer i,n
  real*8 fi,ri(2)

  formk=0
  n=Lformfunc(iform)%nr
  do i=1,n
     ri=Lformfunc(iform)%r(1:2,i)
     fi=Lformfunc(iform)%f(i)
     formk=formk+fi*exp( one*pi*sum(vk*ri) )
  end do

  return
end function formk



