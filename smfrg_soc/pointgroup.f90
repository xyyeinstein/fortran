function reducible(kv,group)
  logical reducible
  real*8 kv(2),kvv(2),a(3),b(3)
  character*3 group
  
  a=(/1,0,0/); b=(/0.5,0.5*sqrt(3.),0./)
  call RBZkv(kv,a,b,kvv)

  reducible=.true.
  select case (group)
	case ('C6v'); if(kvv(2)>=0.and.kvv(1)>=kvv(2)*sqrt(3.d0)*0.999)reducible=.false.
	case ('C3v'); if(kvv(1)>=abs(kvv(2))*sqrt(3.d0)*0.999)reducible=.false.
	case ('C4v'); if(kvv(2)>=0.and.kvv(1)>=kvv(2))reducible=.false.
    case ('C2v'); if(kvv(1)>=0.and.kvv(2)>=0)reducible=.false.
	case default; stop 'unknown group @ reducible'
  end select

  if(sum(kvv*kvv)<1.e-16)reducible=.false.

  return
end function reducible

subroutine GroupAction(ig,p,pg,group)
  implicit none
  integer ig
  real*8 p(3),pg(3)
  character*3 group

  select case (group)
    case ('C2v'); call C2vGroupAction(ig,p,pg)
    case ('C4v'); call C4vGroupAction(ig,p,pg)
    case ('C3v'); call C3vGroupAction(ig,p,pg)
	case ('C6v'); call C6vGroupAction(ig,p,pg)
	case default; stop 'unknown group @ GroupAction'
  end select
  return
end subroutine GroupAction

subroutine FormStructure(ndim,basis,form,formfunc,ng,group)
  use standard_derived_types
  implicit none

  integer ndim                                 !input: spatial dim 
  real*8 basis(ndim)                           !input: basis r vec
  character*3 form                             !input: form of formfunc
  type (Lformconfig) :: formfunc
  integer ng
  character*3 group

  !purpose: get the lattice form structure with a given form

  type starconfig
    integer ndim,nr                            !ndim=spatial dim, nr=# of r vecs in the star
    real*8, allocatable, dimension (:,:) :: r  !r vecs
  end type starconfig
  type (starconfig) :: star                    !internal: star group given by basis r vec

  integer ir,nr,icount
  real*8 value,formvalue

  formfunc%ndim=ndim; formfunc%form=form
  call groupstar(ndim,basis,star,ng,group)

  !counting # of r vecs at which f(r)/=0
  nr=star%nr; icount=0
  do ir=1,nr
     if(abs(formvalue(form,ndim,star%r(:,ir)))>1.e-6)icount=icount+1
  end do

  !saving r vecs at which f(r)/=0
  formfunc%nr=icount
  allocate(formfunc%r(ndim,icount),formfunc%f(icount))
  icount=0
  do ir=1,nr
     value=formvalue(form,ndim,star%r(:,ir))
	 if(abs(value)>1.e-6)then
	   icount=icount+1
       formfunc%r(:,icount)=star%r(:,ir)
	   formfunc%f(icount)=value
	 end if
  end do

  !normalize f(r)
  if(icount==0)return
  nr=icount
  formfunc%f=formfunc%f/sqrt(sum(abs(formfunc%f)**2))

  return
end subroutine formstructure

subroutine groupstar(ndim,basis,star,ng,group)
  implicit none
  integer ndim
  real*8 basis(ndim)
  type starconfig
    integer ndim,nr                            !ndim=spatial dim, nr=# of r vecs in the star
    real*8, allocatable, dimension (:,:) :: r  !r vecs
  end type starconfig
  type (starconfig) :: star
  integer ng
  character*3 group

  !purpose: find the star group dictated by the basis vec

  integer nr
  real*8 b(3)
  real*8, dimension (3,ng) :: r

  integer ig,iig
  real*8 overlap

  b=0; b(1:ndim)=basis; star%ndim=ndim

  do ig=1,ng; call GroupAction(ig,b,r(:,ig),group); end do

  do ig=2,ng
     do iig=1,ig-1
	    if(sum(r(:,iig))>=3.d10)cycle                                   !skip trashed r vecs
		if( abs( sum(r(:,iig)*r(:,ig))-sum(r(:,iig)*r(:,iig)) )<1.e-5)r(:,ig)=1.d10  !trash duplicated r vec
     end do
  end do

  !collect independent rvectors
  nr=1
  do ig=2,ng
	 if(sum(r(:,ig))>=3.d10)cycle                                       !skip trashed r vec
	 nr=nr+1
	 r(:,nr)=r(:,ig)
  end do
  
  star%nr=nr
  allocate(star%r(ndim,nr))
  star%r=r(1:ndim,1:nr)

  return
end subroutine groupstar

function formvalue(form,ndim,r)
  implicit none
  real*8 formvalue             !value of formfunc at position r
  character*3 form             !symmetry of the form func
  integer ndim                 !spatial dim of r vec
  real*8 r(ndim),x,y

  !purpose: find f(r) for a given form (unnormalized)

  formvalue=0; x=r(1); y=r(2)
  select case (form)
    case ('A1g'); formvalue=1
    case ('Ax2'); formvalue=x*x
    case ('Ay2'); formvalue=y*y
	case ('E1u'); formvalue=x
	case ('E2u'); formvalue=y
	case ('Y+X'); formvalue=(y+x)/2
	case ('Y-X'); formvalue=(y-x)/2
	case ('B1g'); formvalue=x*x-y*y
	case ('B2g'); formvalue=x*y
    case default; stop 'error: unknown form @ formvalue'
  end select
  return
end function formvalue

subroutine orbitimage(ig,group,orbit,orbitg,Xg)
  implicit none
  integer ig                    !index of group element
  character*3 orbit,orbitg      !orbit table
  character*3 group      !name of point group
  real*8 Xg              !character of iorb under group element ig, applicable for 1-dim representations

  real*8, dimension (3) :: b,bg   !vectors used to construct the image 
  integer iig,nng
  real*8  formvalue
  character*3 form

  form=orbit
  select case (form)
    case ('A1g'); orbitg=form; Xg=1
    case ('Ax2'); orbitg=form; Xg=1
	case ('Ay2'); orbitg=form; Xg=1
	case ('B1g'); orbitg=form; b=(/1,0,0/); call groupaction(ig,b,bg,group); Xg=formvalue(form,3,bg)
    case ('B2g'); orbitg=form; b=(/1,1,0/); call groupaction(ig,b,bg,group); Xg=formvalue(form,3,bg)
	case ('E1u') 
	   b=(/1,0,0/); call groupaction(ig,b,bg,group)
	   select case (group)
         case ('C4v')
		   Xg=sum(b*bg)
		   if(abs(Xg)>0.5)then; Xg=bg(1); orbitg=form; else
		   b=(/0,1,0/); call groupaction(ig,b,bg,group); Xg=bg(1); orbitg='E2u'; end if
		 case  ('C2v')				
	       orbitg='E1u'; Xg=formvalue(form,3,bg)   
		 case default; stop 'unknown group for E1u @ orbitimage'
	   end select
	case ('E2u')
	   b=(/0,1,0/); call groupaction(ig,b,bg,group)
	   select case (group)
         case ('C4v')
	       Xg=sum(b*bg)
		   if(abs(Xg)>0.5)then; orbitg=form; else
		   b=(/1,0,0/); call groupaction(ig,b,bg,group); Xg=bg(2); orbitg='E1u'; end if   
		 case ('C2v')
	       orbitg='E2u'; Xg=formvalue(form,3,bg)
		 case default; stop 'unknown group for E2u @ orbitimage'
	   end select 

	case ('Y+X') 
	   b=(/1,1,0/); call groupaction(ig,b,bg,group)
	   select case (group)
         case ('C4v')
		   Xg=formvalue('Y+X',3,bg)
		   if(abs(Xg)>0.5)then; orbitg=form; else
		   b=(/-1,1,0/); call groupaction(ig,b,bg,group); Xg=formvalue('Y+X',3,bg); orbitg='Y-X'; end if
		 case  ('C2v')				
	       orbitg='Y+X'; Xg=formvalue(form,3,bg)   
		 case default; stop 'unknown group for Y+X @ orbitimage'
	   end select
	case ('Y-X')
	   b=(/-1,1,0/); call groupaction(ig,b,bg,group)
	   select case (group)
         case ('C4v')
	       Xg=formvalue('Y-X',3,bg)
		   if(abs(Xg)>0.5)then; orbitg=form; else
		   b=(/1,1,0/); call groupaction(ig,b,bg,group); Xg=formvalue('Y-X',3,bg); orbitg='Y+X'; end if   
		 case ('C2v')
	       orbitg='Y-X'; Xg=formvalue(form,3,bg)
		 case default; stop 'unknown group for Y-X @ orbitimage'
	   end select 

	case default; stop 'unknown form @ orbitimage'
  end select 

  return
end subroutine orbitimage

subroutine SpinImage(ig,group,SpinIndexg,SpinXg)
  use standard_derived_types
  implicit none
  integer ig,SpinIndexg(2)
  character*3 group
  complex*16 SpinXg(2)

  integer iig
  real*8 theta,pi
  complex*16 one,factor

  !notice: the transformation law for spin is model dependent. 
  !        It depends on the spin-orbit interaction and the point group. 
  !        The present implementation applies for C(theta)=exp(-i*theta*sig3/2) and m_y=sig1

  !purpose: calculate the transform matrix wrt spin under group element ig of C2v
  !notice:  Xg defined as: g s(m) = Xg(m,n) s(n) where s is an arbitrary spin wave function. 
  !         Since Xg has either diag or off-diag elements, it is sufficient to store the 
  !         nonzero column of each row as,
  !                 Xg(k) = nonzero column element in k-th row
  !                 indexg(k) = index of nonzero column element in k-th row
  
  select case (group)
    case ('C2v')
      iig=mod(ig-1,2)+1
	  !pure C2 rotation = sig_3 in Lee-model
      if(iig==1)then
        SpinIndexG=(/1,2/); SpinXg=1
      else 
        SpinIndexG=(/1,2/); SpinXg=(/1,-1/)
      end if
	  !y-mirror or x-mirror = sig_1 in Lee-model
      if(iig/=ig)then
        SpinIndexG=(/2,1/)
        SpinXg=(/SpinXg(2),SpinXg(1)/)
      end if
    case ('C4v')
	  pi=asin(1.d0)*2; one=cmplx(0,1)
      iig=mod(ig-1,4)+1
	  !pure C4 rotation = exp(-i*theta*sig_3/2) in Lee-model
	  theta=pi*(iig-1)/2.; factor=exp(-one*theta/2)
      SpinIndexG=(/1,2/); SpinXg(1)=factor; spinXg(2)=conjg(factor)
	  !y-mirror or x-mirror = sig_1 in Lee-model
      if(iig/=ig)then
        SpinIndexG=(/2,1/)
        SpinXg=(/conjg(factor),factor/)
      end if

    case default; stop 'spinimage undefined @ LeeSpinImage'
  end select	  
  return
end subroutine SpinImage

function atomimage(ig,iatom,model)
  use standard_derived_types
  implicit none
  integer ig,iatom,atomimage
  type (modelblock) :: model

  integer i,ndim
  integer findatom         
  !findatom is an external function to find an atom by its position in the super unit cell
  real*8 r(3),rg(3),dr(3)

  if(iatom<=0.or.iatom>model%natom)stop 'iatom error @ atomimage4Fe'

  r=model%ra(:,iatom)
  call groupaction(ig,r,rg,model%group)
  atomimage=findatom(rg,model); if(atomimage<=0)stop 'atomimage can not be found'

  return
end function atomimage

subroutine BlochImage(ig,norb,uk,ukg,model)
  use standard_derived_types
  implicit none

  integer ig,norb
  complex*16, dimension (norb) :: uk,ukg
  type (modelblock) :: model

  integer natom,norb1atom
  integer iatom,ispin,idim
  integer iiatom,iispin,iidim
  integer SpinIndexg(2)
  complex*16 SpinXg(2)
  integer atomimage

  natom=model%natom
  norb1atom=norb/natom
  if(norb1atom>2)stop 'norb1atom undefined @ BlochImage'
  
  if(norb1atom==2)call SpinImage(ig,model%group,SpinIndexg,SpinXg)

  do iatom=1,natom
     iiatom=atomimage(ig,iatom,model)
     if(norb1atom==1)then
	   ukg(iiatom)=uk(iatom)
	 else
       do ispin=1,2; idim=iatom+(ispin-1)*natom
	      iispin=SpinIndexg(ispin)
		  iidim=iiatom+(iispin-1)*natom
		  ukg(idim)=uk(iidim)*SpinXg(ispin)
       end do
	 end if
  end do

  return
end subroutine BlochImage

subroutine C6vGroupAction(ig,rv,rg)
    implicit none
    integer ig,iig
    real*8 rv(3),rg(3)
    !group elements: (E,R1,R2,...,R5) X (E,y-refl)
	integer i

    iig=mod(ig-1,6)+1; rg=rv
    do i=1,iig-1; call rotate60(rg(1:2)); end do
    if(ig/=iig)rg(2)=-rg(2)
    return
end subroutine C6vGroupAction

subroutine C3vGroupAction(ig,rv,rg)
    implicit none
    integer ig,iig
    real*8 rv(3),rg(3)
    !group elements: (E,R1,R2) X (E,x-refl)
	integer i

    iig=mod(ig-1,3)+1; rg=rv
    do i=1,iig-1; call rotate120(rg(1:2)); end do
    if(ig/=iig)rg(1)=-rg(1)
    return
end subroutine C3vGroupAction

subroutine C2vGroupAction(ig,rv,rg)
    implicit none
    integer ig,nng,iig
    real*8 rv(3),rg(3)
    !group elements: (E,R180)*(E,m_y)
	nng=2; iig=mod(ig-1,nng)+1
    if(iig==1)then; rg=rv; else; rg=-rv; rg(3)=rv(3); end if
	!if(iig/=ig)rg(2)=-rg(2)
	if(iig/=ig)rg(1)=-rg(1)
    return
end subroutine C2vGroupAction

subroutine C4vGroupAction(ig,rv,rg)
    implicit none
    integer ig,iig
    real*8 rv(3),rg(3)

    !group elements: (E,C4,C2,C4^3) X (E,m_y)
    
	integer i

	iig=mod(ig-1,4)+1
    rg=rv
    do i=2,iig; call rotate90(rg(1:2)); end do
    if(ig/=iig)rg(1)=-rg(1)
    return
end subroutine C4vGroupAction

subroutine rotate90(r)
  implicit none
  real*8 r(2),x
  x=r(1)
  r(1)=-r(2)
  r(2)=x
  return
end subroutine rotate90

subroutine rotate45(r)
  implicit none
  real*8 r(2),work(2)
  work=r
  r(1)=work(1)-work(2)
  r(2)=work(1)+work(2)
  r=r/sqrt(2.d0)
  return
end subroutine rotate45

subroutine rotate60(r)
  implicit none
  real*8 r(2),work(2),rt3
  work=r; rt3=sqrt(3.d0)
  r(1)=0.5d0*( work(1)- work(2)*rt3 )
  r(2)=0.5d0*( rt3*work(1) + work(2) )
  return
end subroutine rotate60

subroutine rotate120(r)
  implicit none
  real*8 r(2)
  integer i
  do i=1,2; call rotate60(r); end do
  return
end subroutine rotate120
