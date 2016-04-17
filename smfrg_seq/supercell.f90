subroutine RBZkv(kv,a,b,kvv)
  implicit none
  real*8, dimension (2) :: kv,kvv
  real*8 a(3),b(3),c(3),da(3),db(3)
  real*8, dimension (2) :: ka,kb,work
  real*8 xa,xb
  integer ia,ib

  !notice: RBZkv depends on a,b,c and models.

  if(abs(sum(a(1:2)*b(1:2)))<1.e-6)then
    !this applies for all models defined on square lattices
    kvv=kv; where(kvv<-1)kvv=kvv+2; where(kvv>1)kvv=kvv-2; return
  end if

  c=(/0,0,1/)
  call dualvector(b,c,a,da); call dualvector(c,a,b,db)

  ka=da(1:2)*2;  kb=db(1:2)*2
  xa=mod(sum(kv*a(1:2))/2,1.d0); if(xa>0.5)xa=xa-1; if(xa<-0.5)xa=xa+1
  xb=mod(sum(kv*b(1:2))/2,1.d0); if(xb>0.5)xb=xb-1; if(xb<-0.5)xb=xb+1
  kvv=xa*ka+xb*kb

  do ia=-2,2; do ib=-2,2
     work=kvv+ia*ka+ib*kb
	 if(sum(work*work)<sum(kvv*kvv))then; kvv=work; end if
  end do; end do

  !print*,kv,kvv
  !print*,sum((kv-kvv)*a(1:2))/2,sum((kv-kvv)*b(1:2))/2
  !pause 'kv,kvv,kv-kvv ok?'

  return
end subroutine RBZkv

subroutine hexbounds()
  implicit none

  real*8 q(2)
  integer i

  q=(/4./3.,0./)
  open(10,file='hexbz.dat')
  write(10,*)q
  do i=1,6; call rotate60(q); write(10,*)q; end do
  close(10)

  return
end subroutine hexbounds

subroutine modifyindex(i1,a1,i2,a2)
  use workplace, only : norb,natom,model
  integer i1,a1,i2,a2
  if(model%sample=='CuO2')return
  i1=i1+(a1-1)*(norb/natom)
  i2=i2+(a2-1)*(norb/natom)
  return
end subroutine modifyindex 

function findatom(rv)
  use workplace, only : model
  implicit none
  integer findatom
  real*8, dimension (3) :: rv

  real*8, dimension (3) :: dr,a,b,c,ab,bc,ca
  real*8, pointer, dimension (:,:) :: ra
  integer natom,i,ia,ib,ic

  natom=model%natom
  if(natom==1)then; findatom=1; return; end if

  a=model%a; b=model%b; c=model%c
  call dualvector(a,b,c,ab)
  call dualvector(b,c,a,bc)
  call dualvector(c,a,b,ca)
  ra=>model%ra
  findatom=-1
  do i=1,natom; dr=rv-ra(:,i)
     ia=nint(sum(dr*bc)); ib=nint(sum(dr*ca))
	 ic=0;  if(model%periodiclayers)ic=nint(sum(dr*ab)) 
     dr=dr-ia*a-ib*b-ic*c
	 if(sum( abs(dr) )<1.e-5)then
		findatom=i; return
	 end if
  end do
  !print*,rv
  !stop 'atom not found'
  return
end function findatom


subroutine dualvector(a,b,c,ab)
  implicit none
  real*8, dimension (3) :: a,b,c,ab

  ab(1)=a(2)*b(3)-a(3)*b(2)
  ab(2)=a(3)*b(1)-a(1)*b(3)
  ab(3)=a(1)*b(2)-a(2)*b(1)

  ab=ab/sum(ab*c)
  return
end subroutine dualvector
     
