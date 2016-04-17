
function findatom(rv,model)
  use standard_derived_types
  implicit none
  integer findatom
  real*8, dimension (3) :: rv
  type (modelblock), target :: model

  real*8, dimension (3) :: dr,a,b,c,ab,bc,ca
  real*8, pointer, dimension (:,:) :: ra
  integer natom,i,ia,ib,ic

  natom=model%natom
  if(natom==1)then; findatom=1; return; end if

  a=model%a; b=model%b; c=model%c
  call dualvector(a,b,c,ab)
  call dualvector(b,c,a,bc)
  call dualvector(c,a,b,ca)

  findatom=-1
  ra=>model%ra
  do i=1,natom; dr=rv-ra(:,i)
     ia=nint(sum(dr*bc)); ib=nint(sum(dr*ca))
     dr=dr-ia*a-ib*b
	 if(sum( abs(dr) )<1.e-5)then
		findatom=i; return
	 end if
  end do
  !print*,rv
  !print*,ia,ib
  !stop 'atom not found at specified position'
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

subroutine RBZkv(kv,a,b,kvv)
  implicit none
  real*8, dimension (2) :: kv,kvv
  real*8 a(3),b(3),c(3),da(3),db(3)
  real*8, dimension (2) :: ka,kb,work
  real*8 xa,xb
  integer ia,ib

  c=(/0,0,1/)
  call dualvector(b,c,a,da); call dualvector(c,a,b,db)

  ka=da(1:2)*2;  kb=db(1:2)*2
  xa=mod(sum(kv*a(1:2))/2,1.d0); if(xa>0.5)xa=xa-1; if(xa<-0.5)xa=xa+1
  xb=mod(sum(kv*b(1:2))/2,1.d0); if(xb>0.5)xb=xb-1; if(xb<-0.5)xb=xb+1
  kvv=xa*ka+xb*kb

  if(abs(sum(a*b))<1.e-5)return  !applied for square lattice

  do ia=-2,2; do ib=-2,2
     work=kvv+ia*ka+ib*kb
	 if(sum(work*work)<sum(kvv*kvv))kvv=work
  end do; end do

  return
end subroutine RBZkv

