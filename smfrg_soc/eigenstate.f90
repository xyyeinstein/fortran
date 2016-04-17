
function findek(nativeband,kv,norb,model)
  use standard_derived_types
  implicit none

  real*8 findek
  integer nativeband,norb
  real*8 kv(2)
  type (modelblock) :: model

  complex*16 hk(norb,norb)
  real*8 eval(norb)

  call gethk(norb,kv,hk,model)
  call ZHEIGEN(norb,hk,eval) 
  findek=eval(nativeband)

  return
end function findek

subroutine gethk(norb,kv,hk,model)
  use standard_derived_types
  implicit none
  integer norb
  real*8 kv(2)
  complex*16 hk(norb,norb)
  type (modelblock) :: model

  select case (norb)
    case (2); if(model%sample=='cuprt')call cupratehk(norb,kv,hk,model)
	case (4); if(model%sample=='cuprate2L')call cuprate2Lhk(norb,kv,hk,model)
			  if(model%sample=='lee')call leehk(norb,kv,hk,model)
	case default; stop 'model hk unknown'
  end select
  return
end subroutine gethk




