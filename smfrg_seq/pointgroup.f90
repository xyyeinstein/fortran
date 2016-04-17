function combine(form1,form2,group)
  character*3 form1,form2,combine,group

  !purpose: combine two 1d representations to a new 1d representation
  !notice: the combined form must be 1d irreducible for square lattice

  character*3 formtable(4),biform(4,4)

  if(group=='C2v'.or.group=='D2h')then
    formtable=(/'A1g','B2g','E1u','E2u'/)
    biform(:,1)=(/'A1g','B2g','E1u','E2u'/)
    biform(:,2)=(/'B2g','A1g','E2u','E1u'/)
    biform(:,3)=(/'E1u','E2u','A1g','B2g'/)
    biform(:,4)=(/'E2u','E1u','B2g','A1g'/)
  else
    formtable=(/'A1g','A2g','B1g','B2g'/)
    biform(:,1)=(/'A1g','A2g','B1g','B2g'/)
    biform(:,2)=(/'A2g','A1g','B2g','B1g'/)
    biform(:,3)=(/'B1g','B2g','A1g','A2g'/)
    biform(:,4)=(/'B2g','B1g','A2g','A1g'/)
  end if
  iform=findform(form1); if(iform<=0)stop 'unknown 1d form @ combine'
  jform=findform(form2); if(jform<=0)stop 'unknown 1d form @ combine'

  combine=biform(iform,jform)  
  return
  contains
  function findform(form)
    integer findform
	character*3 form
	integer i
	findform=-1
	do i=1,4; if(form==formtable(i))findform=i; end do
	return
  end function findform
end function combine


function compatible(A1,A2,A,ng,group)
  character*3 A1,A2,A          
  integer ng
  character*3 group
  logical compatible

  !purpose: check whether A1*A2=A for 1d representations A1,A2 and A, under the specified group for square lattices
  !notice: For low symmetry groups, different 1d reps of higher symmetry groups may mix.
  !        Therefore compatibility is checked using group characters rather than the direct 
  !        product formular 'combine'
  !notice: compatible=.true. if A is specified as '***' (i.e., arbitrary).

  real*8, dimension (ng) :: X1,X2,X
  integer ig,indexg(2)
  real*8 Xg(2)

  compatible=.true.
  if(A=='***')return

  do ig=1,ng
     call orbitimage(ig,1,1,(/A1/),group,indexg,Xg); X1(ig)=Xg(1)
     call orbitimage(ig,1,1,(/A2/),group,indexg,Xg); X2(ig)=Xg(1)
     call orbitimage(ig,1,1,(/A/),group,indexg,Xg); X(ig)=Xg(1)
  end do

  compatible=( sum( abs(X1*X2-X) )<1.e-6 )
  return
end function compatible


function reducible(kv,group)
  use workplace, only : model
  implicit none
  logical reducible
  real*8 kv(2),kvv(2)
  character*3 group
  
  !call RBZkv(kv,model%a,model%b,kvv)
  kvv=kv;
  reducible=.true.
  select case (group)
    case ('D4h'); if(kvv(1)>=0.and.kvv(2)>=0.and.kvv(1)>=kvv(2))reducible=.false.
    case ('C4v'); if(kvv(1)>=0.and.kvv(2)>=0.and.kvv(1)>=kvv(2))reducible=.false.
    case ('D2h'); if(kvv(1)>=0.and.kvv(2)>=0)reducible=.false. !if(kvv(2)>=0.and.kvv(2)>=abs(kvv(1)))reducible=.false.
    case ('C2v'); if(kvv(1)>=0.and.kvv(2)>=0)reducible=.false. !if(kvv(2)>=0.and.kvv(2)>=abs(kvv(1)))reducible=.false.
	case ('C6v'); if(kvv(2)>=0.and.kvv(1)>=kvv(2)*sqrt(3.d0)*0.99999)reducible=.false.
	case ('C6z'); if(kvv(1)>=abs(kvv(2))*sqrt(3.d0)*0.99999)reducible=.false.
	case ('C3v'); if(kvv(1)>=abs(kvv(2))*sqrt(3.d0)*0.99999)reducible=.false.
	case ('S4_'); if(kvv(1)>=0.and.kvv(2)>=0)reducible=.false.
    case ('S4A'); if(kvv(1)>=0.and.kvv(2)>=0.and.kvv(1)>=kvv(2))reducible=.false.
	case ('S4B'); if(kvv(1)>=0.and.kvv(2)>=0.and.kvv(1)>=kvv(2))reducible=.false.
	case ('C2_'); if(kvv(1)>=0)reducible=.false.
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
    case ('D4h'); call D4hGroupAction(ig,p,pg)
	case ('D2h'); call D2hGroupAction(ig,p,pg)
	case ('C4v'); call C4vGroupAction(ig,p,pg)
	case ('C2v'); call C2vGroupAction(ig,p,pg)
	case ('C2_'); call C2GroupAction(ig,p,pg)
    case ('C3v'); call C3vGroupAction(ig,p,pg)
	case ('C6v'); call C6vGroupAction(ig,p,pg)
	case ('C6z'); call C6vGroupAction(ig,p,pg)   
	case ('S4_'); call S4GroupAction(ig,p,pg)
	case ('S4A'); call S4AGroupAction(ig,p,pg)
	case ('S4B'); call S4BGroupAction(ig,p,pg)
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
  real*8 value
  real*8 formvalue

  formfunc%ndim=ndim; formfunc%form=form
  call groupstar(ndim,basis,star,ng,group)

  !counting # of r vecs at which f(r)/=0
  nr=star%nr; icount=0
  do ir=1,nr
     if(formvalue(form,ndim,star%r(:,ir))/=0)icount=icount+1
  end do

  !saving r vecs at which f(r)/=0
  formfunc%nr=icount
  allocate(formfunc%r(ndim,icount),formfunc%f(icount))
  icount=0
  do ir=1,nr
     value=formvalue(form,ndim,star%r(:,ir))
	 if(value/=0)then
	   icount=icount+1
       formfunc%r(:,icount)=star%r(:,ir)
	   formfunc%f(icount)=value
	 end if
  end do

  !normalize f(r)
  if(icount==0)return
  nr=icount
  formfunc%f=formfunc%f/sqrt(sum(formfunc%f**2))

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
		if(sum(r(:,iig)*r(:,ig))==sum(r(:,iig)*r(:,iig)))r(:,ig)=1.d10  !trash duplicated r vec
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
  real*8 r(ndim),x,y,z

  !purpose: find f(r) for a given form (unnormalized)

  formvalue=0
  x=r(1); y=r(2); if(ndim==3)z=r(3)
  select case (form)
    case ('A1g'); formvalue=1
	case ('Odd'); formvalue=1
	   x=r(1); y=r(2); if(abs(x)>=abs(y))then; if(x<0)formvalue=-1; else; if(y<0)formvalue=-1; end if
	case ('A2g'); formvalue=x*y*( x*x - y*y )
	case ('A1u'); formvalue=z
	case ('B1g'); formvalue=x*x-y*y
	case ('B2g'); formvalue=2*x*y
	case ('E1g'); formvalue=x*z
	case ('E2g'); formvalue=y*z
	case ('E1u'); formvalue=x
	case ('E2u'); formvalue=y
	case ('E3u'); formvalue=x*(x*x-3*y*y)
    case default; stop 'error: unknown form @ formvalue'
  end select
  return
end function formvalue


subroutine orbitimage(ig,iorb,norb,orbits,group,iorbg,Xg)
  implicit none
  integer ig             !index of group element
  integer iorb           !orbit index pointing to an orbit table
  integer norb           !number of orbits in orbits
  character*3 orbits(*)  !orbit table
  character*3 group      !name of point group
  integer iorbg(2)       !the image of iorb under group element ig, applicable for 1- and 2-dim representations
  real*8 Xg(2)           !character of iorb under group element ig, applicable for 1- and 2-dim representations
  integer orbg,orbitimage_square

  if(group=='C3v'.or.group=='C6v'.or.group=='C6z')then
    call orbitimage_hex(ig,iorb,norb,orbits,group,iorbg,Xg)
  else
    orbg=orbitimage_square(ig,iorb,norb,orbits,group)
	iorbg=(/abs(orbg),0/); Xg=orbg/abs(orbg)
  end if

  return
end subroutine orbitimage

function orbitimage_square(ig,iorb,norb,orbits,group)
  implicit none
  integer ig             !index of group element
  integer iorb           !orbit index pointing to an orbit table
  integer norb           !number of orbits in orbits
  character*3 orbits(*)  !orbit table
  character*3 group      !name of point group
  integer orbitimage_square     !the image of iorb under group element ig

  !purpose: find the image of orbits(iorb) under group element ig
  !notice:  the table orbits(:) must be irreducible in the sense that a particular type of orbit can not be repeated in the table. 
  !notice:  this subprogram requires that the mapping is 1-to-1, applicable on square lattices only

  character*3 form
  integer image

  real*8, dimension (3) :: b,bg   !vectors used to construct the image 
  integer overlap                 !overlap between a vector and its image vector
  integer i,icount

  if(iorb<1.or.iorb>norb)stop 'invalid input iorb @ orbitimage'
  !if(norb<=0.or.norb>8)stop 'invalid input norb @ orbitimage'  
  !there are at most 8 orbit forms in the D4h and D2h group

  form=orbits(iorb)
  image=0                !image=0 if it can not be identified
  
  select case (form)
    case ('A1g')
	  image=iorb
    case ('A2g')
	  b=(/1,0,0/)
	  call GroupAction(ig,b,bg,group)
	  image=iorb*( bg(1)*bg(1)-bg(2)*bg(2) )	  
	  b=(/1,1,0/)
	  call GroupAction(ig,b,bg,group)
	  image=image*( bg(1)*bg(2) )
	case ('A1u')
	  b=(/0,0,1/)
	  call GroupAction(ig,b,bg,group)
	  image=iorb; if(sum(b*bg)<0)image=-iorb
	case ('B1g')
	  b=(/1,0,0/)
	  call GroupAction(ig,b,bg,group)
	  image=iorb*( bg(1)*bg(1)-bg(2)*bg(2) )	  
	case ('B2g')
	  b=(/1,1,0/)
	  call GroupAction(ig,b,bg,group)
	  image=iorb*( bg(1)*bg(2) )
	case ('E1g')   !E1g=E1u*A1u
	  b=(/1,0,0/)
      call GroupAction(ig,b,bg,group)
	  overlap=sum(b*bg)
	  if(overlap/=0)then
	    image=iorb*overlap
	  else
	    icount=0
        do i=1,norb; if(orbits(i)=='E2g')then; image=i; icount=icount+1; end if; end do
		if(icount==0)stop 'No E2g found as partner of E1g'
		if(icount>1)stop 'Too many E2g found as partner of E1g'
	    image=image*sum(bg*(/0,1,0/))
      end if
	  b=(/0,0,1/)
	  call GroupAction(ig,b,bg,group)
	  image=image*sum(b*bg)
	case ('E2g')   !E2g=E2u*A1u
	  b=(/0,1,0/)
      call GroupAction(ig,b,bg,group)
	  overlap=sum(b*bg)
	  if(overlap==1.or.overlap==-1)then
	    image=iorb*overlap
	  else
	    icount=0
        do i=1,norb; if(orbits(i)=='E1g')then; image=i; icount=icount+1; end if; end do
		if(icount==0)stop 'No E1g found as partner of E2g'
		if(icount>1)stop 'Too many E1g found as partner of E2g'
	    image=image*sum(bg*(/1,0,0/))
      end if
	  b=(/0,0,1/)
	  call GroupAction(ig,b,bg,group)
	  image=image*sum(b*bg)
	case ('E1u')
	  b=(/1,0,0/)
      call GroupAction(ig,b,bg,group)
	  overlap=sum(b*bg)
	  if(overlap/=0)then
	    image=iorb*overlap
	  else
	    icount=0
        do i=1,norb; if(orbits(i)=='E2u')then; image=i; icount=icount+1; end if; end do
		if(icount==0)stop 'No E2u found as partner of E1u'
		if(icount>1)stop 'Too many E2u found as partner of E1u'
	    image=image*sum(bg*(/0,1,0/))
      end if
	case ('E2u')
	  b=(/0,1,0/)
      call GroupAction(ig,b,bg,group)
	  overlap=sum(b*bg)
	  if(overlap/=0)then
	    image=iorb*overlap
	  else
	    icount=0
        do i=1,norb; if(orbits(i)=='E1u')then; image=i; icount=icount+1; end if; end do
		if(icount==0)stop 'No E1u found as partner of E2u'
		if(icount>1)stop 'Too many E1u found as partner of E2u'
	    image=image*sum(bg*(/1,0,0/))
      end if
	case default; print*,form; stop 'error: unknown orbit form @ orbitimage'
  end select

  orbitimage_square=image
  return
end function orbitimage_square


subroutine orbitimage_hex(ig,iorb,norb,orbits,group,iorbg,Xg)
  implicit none
  integer ig             !index of group element
  integer iorb           !orbit index pointing to an orbit table
  integer norb           !number of orbits in orbits
  character*3 orbits(*)  !orbit table
  character*3 group      !name of point group
  integer iorbg(2)       !the image of iorb under group element ig, applicable for 1- and 2-dim representations
  real*8 Xg(2)           !character of iorb under group element ig, applicable for 1- and 2-dim representations

  !purpose: find the image and character of orbits(iorb) under group element ig

  !notice: For combined group operations g = refl * R(theta), where refl is either
  !        x- or y-reflection, theta is the rotation angle, the character matrix is
  !        given by 
  !                   X(g) = X( refl ) * X( theta ).
  !    *a) for C3v, L=2, X( xrefl ) = sig_3, and its effect in X(g) is to change 
  !        the sign of the second row of X(theta).
  !    *b) for C3v, L=1, X( xrefl ) = -sig_3, so its effect in X(g) is to change
  !        the sign of the first row of X(theta).
  !    *c) for C6v, L=2, X ( yrefl ) = sig_3, so its effect in X(g) is to change
  !        the sign of the second row of X(theta), as in *a).
  !    *d) for C6v, L=1, X ( yrefl ) = sig_3, so its effect in X(g) is to change 
  !        the second row of X(theta).
  !notice: exp(3*i*theta) forms two separate 1dim representations on short bonds of
  !        hexagonal lattices, therefore we do not consider E3c and E3s as conjugated 
  !        pair of angular momentum L=3.    

  character*3 form,formi

  real*8, dimension (3) :: b,bg   !vectors used to construct the image 
  integer iig,nng
  real*8  formvalue
  integer Eorb(2),irow,iL
  real*8 twopi,theta

  if(iorb<1.or.iorb>norb)stop 'invalid input iorb @ orbitimage'

  form=orbits(iorb)

  !----------1dim representations----------
  if(form=='A1g')then; Xg=(/1,0/); iorbg=(/iorb,0/); return; end if
  if(form=='Odd')then; Xg=(/1,0/); if(ig==2)Xg=(/-1,0/); iorbg=(/iorb,0/); return; end if

  if(form=='E3u')then
     b=(/1,0,0/); call groupaction(ig,b,bg,group)
	 iorbg=(/iorb,0/); Xg=(/formvalue(form,3,bg),0.d0/); return
   end if
  !----------------------------------------

  !----------2dim representations----------------
  twopi=asin(1.d0)*4

  if(group=='C3v')then; nng=3; else if(group=='C6v')then; nng=6; end if
  
  iig=mod(ig-1,nng)+1; theta=(iig-1)*twopi/nng

  call getEorbits()
  if(irow==1)then; Xg=(/cos(theta*iL),-sin(theta*iL)/); 
  else if(irow==2)then; Xg=(/sin(theta*iL),cos(theta*iL)/); 
  else; stop 'error: E-orbits do not come in pair.'; end if

  if(iig/=ig)then
    select case (group)
	  case ('C3v')
        !    *a) for C3v, L=2, X( xrefl ) = sig_3, and its effect in X(-g) is to change 
        !        the sign of the second row of X(theta). 
        !    *b) for C3v, L=1 or 3, X( xrefl ) = -sig_3, so its effect in X(-g) is to change
        !        the sign of the first row of X(theta).
        if(iL==2.and.irow==2)Xg=-Xg
		if(iL==1.and.irow==1)Xg=-Xg
      case ('C6v')
        !    *c) for C6v, L=2, X ( yrefl ) = sig_3, so its effect in X(-g) is to change
        !        the sign of the second row of X(theta), as in *a).
        !    *d) for C6v, L=1 or 3, X ( yrefl ) = sig_3, so its effect in X(-g) is to change 
        !        the second row of X(theta).
		if(iL==2.and.irow==2)Xg=-Xg
		if(iL==1.and.irow==2)Xg=-Xg
	  case default; stop 'unknown reflection in unknown group @ orbitimage'
    end select
  end if
  return
  contains 
  subroutine getEorbits()
	select case (form)
	case ('E1u'); iorbg(1)=iorb; iorbg(2)=findorbit('E2u'); irow=1; iL=1
	case ('E2u'); iorbg(2)=iorb; iorbg(1)=findorbit('E1u'); irow=2; iL=1
	case ('B1g'); iorbg(1)=iorb; iorbg(2)=findorbit('B2g'); irow=1; iL=2
	case ('B2g'); iorbg(2)=iorb; iorbg(1)=findorbit('B1g'); irow=2; iL=2
    case default; stop 'unknown conjugated pair @ orbitimage'
	end select
    return
  end subroutine getEorbits

  function findorbit(name)
    character*3 name
	integer i,findorbit
	do i=1,norb; if(orbits(i)==name)then; findorbit=i; return; end if; end do
	stop 'orbit not found @ findorbit'
	return
  end function findorbit
end subroutine orbitimage_hex


subroutine OformImage(ig,iOform,nOform,Oformfunc,morb,orbits,group,indexOg,XOg)
  use standard_derived_types
  implicit none
  integer ig,iOform,nOform,morb
  type (Oformconfig), target :: Oformfunc(nOform)
  character*3 orbits(morb),group
  integer indexOg(4)
  real*8 XOg(4)

  !purpose: find the group image of a bi-orbital form
  !notice:  make sure that the table orbits contains irreducible orbits only 

  integer orbg(2),orbbg(2)
  real*8  Xorbg(2),Xorbbg(2),work(nOform)
  integer iorb,jOform
  type (Oformconfig), pointer :: Oi,Oj
  integer ip,jp,np,mp,orb,orbb,icount,i,j

  if(iOform<1.or.iOform>nOform)stop 'invalid iOform @ OfromImage'

  Oi=>Oformfunc(iOform); np=Oi%npair
        
  work=0
  do ip=1,np; orb=Oi%pair(1,ip); orbb=Oi%pair(2,ip)
     call orbitimage(ig,orb,morb,orbits,group,orbg,Xorbg)
	 call orbitimage(ig,orbb,morb,orbits,group,orbbg,Xorbbg)
	 do jOform=1,nOform; Oj=>Oformfunc(jOform); mp=Oj%npair    
	    do jp=1,mp; orb=Oj%pair(1,jp); orbb=Oj%pair(2,jp)
	    do i=1,2; if(orbg(i)<=0)cycle
		do j=1,2; if(orbbg(j)<=0)cycle
		   if(orb/=orbg(i).or.orbb/=orbbg(j))cycle
		   work(jOform)=work(jOform)+Oi%f(ip)*Oj%f(jp)*Xorbg(i)*Xorbbg(j)
		end do; end do; end do
	 end do
  end do

  icount=0; indexOg=0; XOg=0
  do jOform=1,nOform
     if(abs(work(jOform))<1.e-6)cycle
     icount=icount+1; if(icount>4)stop 'too many Oform images @ OformImage'
	 indexOg(icount)=jOform; XOg(icount)=work(jOform)
  end do

  if(icount==0)stop 'OformImage not found @ OformImage'

  return
end subroutine OformImage


subroutine LformImage(ig,iLform,nLform,Lformfunc,group,indexg,Xg)
  use standard_derived_types
  implicit none
  integer ig,iLform,nLform
  character*3 group
  type (Lformconfig), target :: Lformfunc(nLform)
  integer indexg(2)
  real*8 Xg(2)

  ! Purpose: find Lj so that g[ Li(r) ] = sign*Lj(r) for r on an invariant space. The information sign*j is encoded as Lformimage.
  ! notice: 1) this subprogram only applies to square lattices for which the transformation is 1-to-1
  !         2) under the group operations, it is understood that the r-vectors form an invariant space (i.e., the support for the form functions),
  !            although only a few of them are necessary in a particular lattice harmonics, such as in E1u or E2u (the form functions may be zero at 
  !            some r-vectors). 
  
  integer nr,mr,ir,jr,jL,icount
  real*8 r(3),rg(3),work(nLform)
  type (Lformconfig), pointer :: Li,Lj
  
  if(iLform<1.or.iLform>nLform)stop 'invalid iLform @ Lformimage'

  Li=>Lformfunc(iLform); nr=Li%nr

  work=0
  do ir=1,nr; r=Li%r(:,ir); call groupaction(ig,r,rg,group)	  		
  do jL=1,nLform; Lj=>Lformfunc(jL); mr=Lj%nr;  do jr=1,mr
     if(sum( abs(rg-Lj%r(:,jr)) )<1.e-5 )work(jL)=work(jL)+Li%f(ir)*Lj%f(jr)
  end do; end do; end do

  icount=0; indexg=0; Xg=0
  do jL=1,nLform
     if( abs(work(jL))<1.e-6 ) cycle
	 icount=icount+1; if(icount>2)stop 'too many Lform images @ LformImage'
	 indexg(icount)=jL; Xg(icount)=work(jL)
  end do

  if(icount==0)then
    print*,'group,ig,iLform=',group,ig,iLform
	print*,Lformfunc(iLform)%r,Lformfunc(iLform)%f
	stop 'LformImage can not be found @ LformImage'
  end if
  return
end subroutine LformImage

function atomimage(ig,iatom,model)
  use standard_derived_types
  implicit none
  integer ig,iatom,atomimage
  type (modelblock) :: model

  integer i,ndim
  integer findatom         
  !findatom is an external function to find an atom by its position in the super unit cell
  real*8 r(3),rg(3),dr(3)

  if(iatom<=0.or.iatom>model%natom)stop 'invalid iatom @ atomimage'

  r=model%ra(:,iatom)
  call groupaction(ig,r,rg,model%group)
  atomimage=findatom(rg)

  return
end function atomimage

subroutine BlochImage(ig,norb,uk,ukg,model)
  use standard_derived_types
  implicit none
  integer ig,norb
  complex*16, dimension (norb) :: uk,ukg
  type (modelblock), target :: model

  integer norb1atom,iorb,iatom,iiorb,iiatom,idim,iidim,i
  integer indexg(2,norb),atomimage
  real*8  Xg(2,norb)

  integer, allocatable, dimension (:) :: orbg

  if(model%sample=='CuO2')then; call Blochimage_CuO2(ig,norb,uk,ukg,model); return; end if

  if(norb/=model%norb)stop 'invalid norb @ Blochimage'
  if(norb==1)then; ukg=uk; return; end if

  norb1atom=norb/model%natom

  do iorb=1,norb1atom
     call orbitimage(ig,iorb,norb1atom,model%orbits,model%group,indexg(:,iorb),Xg(:,iorb))
  end do

  do iatom=1,model%natom
	 iiatom=atomimage(ig,iatom,model)		
     do iorb=1,norb1atom
	    idim=iorb+(iatom-1)*norb1atom
		do i=1,2; if(indexg(i,iorb)<=0)cycle
	       iiorb=indexg(i,iorb)
	       iidim=iiorb+(iiatom-1)*norb1atom
	       ukg(iidim)=uk(idim)*Xg(i,iorb)
		end do
	 end do
   end do
   return
end subroutine BlochImage		   

subroutine BlochImage_CuO2(ig,norb,uk,ukg,model)
  use standard_derived_types
  implicit none
  integer ig,norb
  complex*16, dimension (norb) :: uk,ukg
  type (modelblock) :: model

  integer iorb,iatom,iiatom
  integer atomimage

  integer indexg(2,norb)
  real*8 Xg(2,norb)

  do iorb=1,norb
     call orbitimage(ig,iorb,norb,model%orbits,model%group,indexg(:,iorb),Xg(:,iorb))
  end do

  do iatom=1,model%natom
	 iiatom=atomimage(ig,iatom,model)
	 if(indexg(1,iatom)/=iiatom)stop 'atom image and orbitimage mismatch @ BlochImage_CuO2'
	 if(indexg(2,iatom)>0)stop 'unexpected orbit image @ BlochImage_CuO2'
	 ukg(iiatom)=uk(iatom)*Xg(1,iatom)
   end do
   return
end subroutine BlochImage_CuO2		   

subroutine D4hGroupAction(ig,rv,rvg)
  implicit none
  integer ig
  real*8, dimension (3) :: rv,rvg

  !purpose: perform D4h group action on rv, yielding rvg
  !notice:  for general purposes, rv and rvg are 3d vecs in this subroutine

  integer iig

  if(ig<1.or.ig>16)stop 'invalid group element @ D4hGroupAction'

  iig=mod(ig-1,8)+1
  select case (iig)
    case (1); rvg=rv
	case (2); rvg(1)=-rv(2); rvg(2)=rv(1);  rvg(3)=rv(3)
	case (3); rvg(1)=rv(2);  rvg(2)=-rv(1); rvg(3)=rv(3)
	case (4); rvg(1)=-rv(1); rvg(2)=-rv(2); rvg(3)=rv(3)
	case (5); rvg(1)=rv(1);  rvg(2)=-rv(2); rvg(3)=-rv(3)
	case (6); rvg(1)=-rv(1); rvg(2)=rv(2);  rvg(3)=-rv(3)
	case (7); rvg(1)=rv(2);  rvg(2)=rv(1);  rvg(3)=-rv(3)
	case (8); rvg(1)=-rv(2); rvg(2)=-rv(1); rvg(3)=-rv(3)
    case default
  end select

  !3d inversion (I)
  if(ig/=iig)rvg=-rvg
  !acting on (x,y,0), I is not an independent group element

  return
end subroutine

subroutine S4AGroupAction(ig,rv,rvg)
  implicit none
  integer ig
  real*8, dimension (3) :: rv,rvg

  !purpose: perform D4h group action on rv, yielding rvg
  !notice:  for general purposes, rv and rvg are 3d vecs in this subroutine

  integer iig

  if(ig<1.or.ig>8)stop 'invalid group element @ D4hGroupAction'

  select case (ig)
    case (1); rvg=rv
	case (2); rvg(1)=-rv(2); rvg(2)=rv(1);  rvg(3)=-rv(3)
	case (3); rvg(1)=rv(2);  rvg(2)=-rv(1); rvg(3)=-rv(3)
	case (4); rvg(1)=-rv(1); rvg(2)=-rv(2); rvg(3)=rv(3)
	case (5); rvg(1)=rv(1);  rvg(2)=-rv(2); rvg(3)=rv(3)
	case (6); rvg(1)=-rv(1); rvg(2)=rv(2);  rvg(3)=rv(3)
	case (7); rvg(1)=rv(2);  rvg(2)=rv(1);  rvg(3)=-rv(3)
	case (8); rvg(1)=-rv(2); rvg(2)=-rv(1); rvg(3)=-rv(3)
    case default
  end select

  return
end subroutine

subroutine S4BGroupAction(ig,rv,rvg)
  implicit none
  integer ig
  real*8, dimension (3) :: rv,rvg

  if(ig<1.or.ig>8)stop 'invalid group element @ D4hGroupAction'

  select case (ig)
    case (1); rvg=rv
	case (2); rvg(1)=-rv(2); rvg(2)=rv(1);  rvg(3)=-rv(3)
	case (3); rvg(1)=rv(2);  rvg(2)=-rv(1); rvg(3)=-rv(3)
	case (4); rvg(1)=-rv(1); rvg(2)=-rv(2); rvg(3)=rv(3)
	case (5); rvg(1)=rv(1);  rvg(2)=-rv(2); rvg(3)=-rv(3)
	case (6); rvg(1)=-rv(1); rvg(2)=rv(2);  rvg(3)=-rv(3)
	case (7); rvg(1)=rv(2);  rvg(2)=rv(1);  rvg(3)=rv(3)
	case (8); rvg(1)=-rv(2); rvg(2)=-rv(1); rvg(3)=rv(3)
    case default
  end select

  return
end subroutine


subroutine D2hGroupAction(ig,p,pg)
  implicit none
  integer ig
  real*8, dimension (3) :: p,pg
  
  if(ig<1.or.ig>4)stop 'invalid ig @ D2hGroupAction'
  
!  select case (ig)
!    case (1); pg=p                    !E
!	case (2); pg=-(/p(2),p(1),p(3)/)  !Dy'
!	case (3); pg=-p                   !I
!	case (4); pg=(/p(2),p(1),p(3)/)   !I*Dy'
!  end select

  select case (ig)
    case (1); pg=p                    !E
	case (2); pg=(/-p(1),p(2),p(3)/)  !mx
	case (3); pg=(/p(1),-p(2),p(3)/)  !my
	case (4); pg=(/-p(1),-p(2),p(3)/) !c2
  end select

  return
end subroutine D2hGroupAction


subroutine C4vGroupAction(ig,rv,rg)
    implicit none
    integer ig,iig
    real*8 rv(3),rg(3)
    !group elements: (E,x-refl,y-refl,inversion) X (E,R90)
    iig=mod(ig-1,4)+1
    select case (iig)
      case (1); rg=rv
	  case (2); rg=(/-rv(1),rv(2),rv(3)/)
	  case (3); rg=(/rv(1),-rv(2),rv(3)/)
	  case (4); rg=(/-rv(1),-rv(2),rv(3)/)
    end select
    if(ig/=iig)rg=(/-rg(2),rg(1),rg(3)/)
    return
end subroutine C4vGroupAction

subroutine C2vGroupAction(ig,p,pg)
  implicit none
  integer ig
  real*8, dimension (3) :: p,pg
  
  if(ig<1.or.ig>4)stop 'invalid ig @ C2vGroupAction'
  
  select case (ig)
    case (1); pg=p                    !E
	case (2); pg=(/-p(1),p(2),p(3)/)  !mx
	case (3); pg=(/p(1),-p(2),p(3)/)  !my
	case (4); pg=(/-p(1),-p(2),p(3)/) !c2
  end select

  return
end subroutine C2vGroupAction

subroutine C2GroupAction(ig,p,pg)
  implicit none
  integer ig
real*8, dimension (3) :: p,pg
  
  if(ig<1.or.ig>2)stop 'invalid ig @ C2vGroupAction'
  
  select case (ig)
    case (1); pg=p                    !E
	case (2); pg=(/-p(1),-p(2),p(3)/) !c2
  end select

  return
end subroutine

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

subroutine S4GroupAction(ig,p,pg)
  implicit none
  integer ig
  real*8, dimension (3) :: p,pg
  if(ig<1.or.ig>4)stop 'invalid ig @ S4GroupAction'

  select case(ig)
    case (1); pg=p
	case (2); pg=(/p(2),-p(1),p(3)/)
	case (3); pg=(/-p(1),-p(2),p(3)/)
	case (4); pg=(/-p(2),p(1),p(3)/)
  end select

  return
end subroutine

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
