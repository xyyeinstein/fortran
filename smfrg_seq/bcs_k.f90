subroutine runBCS_k(cutoff,mft,modecoupled)
  use workplace
  implicit none

  real*8 cutoff
  integer nkf,ik,iq,idim
  complex*16, allocatable, dimension (:,:) :: Peff
  real*8, allocatable, dimension (:) :: eval
  logical mft,modecoupled

  !purpose: self-consistent BCS theory with a cutoff
  !notice:  mft=.false.: only eigen problem is solved
  !         mft=.true.:  both eigen problem and self-consistent theory is solved

  !judge whether band and patch info are available to run this program
  if(.not.allocated(model%npatch))then; print*, 'patch info not specified to runBCS_k'; return; end if

  !judge whether kf, lf, and pdos are already obtained. if not yet, run meshferisurface.
  if(.not.allocated(model%pdos))call meshfermisurface(model)

  nkf=model%nkf

  allocate(Peff(nkf,nkf),eval(nkf))

  call getPeff_k(nkf,model%kf,model%band,Peff,modecoupled)
  !call getPeff(nkf,model%kf,model%band,Peff)
  open(10,file='peff.dat');write(10,'(2e16.8)')Peff; close(10)
  if(.not.model%propergauge)call symmetrize(nkf,Peff,model)
  if(mft)call BCSMFT(nkf,model%pdos,Peff,cutoff)
  call BCS_k(nkf,model%pdos,Peff,eval)
  write(*,'(2e16.8)')eval(1:2)
  open(10,file='lsc.dat'); open(11,file='gapk.dat')
  !output 10 (at most) eigenvalues and eigenfuncs
  do ik=1,nkf; if(ik>10)cycle
     write(10,*)eval(ik),cutoff*exp(-1/eval(ik))
	 write(11,100)Peff(:,ik)
  end do
  close(10); close(11)

100 format(1x,2e20.6)

  !return
  
  ! Permoranchuk instability:

  call getDeff_k(nkf,model%kf,model%band,Peff,modecoupled)  
  if(.not.model%propergauge)call symmetrize(nkf,Peff,model)
  Peff=-Peff; call BCS_k(nkf,model%pdos,Peff,eval)
  
  open(10,file='lcdw.dat'); open(11,file='cdwk.dat')
  do ik=1,nkf; if(ik>10)cycle
     write(10,*)eval(ik)
	 write(11,100)Peff(:,ik)
  end do
  close(10); close(11)


  ! Spin Permoranchuk instability:

  call getCeff_k(nkf,model%kf,model%band,Peff,modecoupled)  
  if(.not.model%propergauge)call symmetrize(nkf,Peff,model)
  Peff=-Peff; call BCS_k(nkf,model%pdos,Peff,eval)

  open(10,file='lsdw.dat'); open(11,file='sdwk.dat')
  do ik=1,nkf;	if(ik>10)cycle
     write(10,*)eval(ik)
	 write(11,100)Peff(:,ik)
  end do
  close(10); close(11)

  print*,'BCS eigenvalues and eigen modes output to file'

  return
end subroutine runBCS_k

subroutine BCS_k(nk,pdos,Peff,eval)
  implicit none
  integer nk                !number of patches on all bands
  real*8 pdos(nk)           !partial density of states on patches
  real*8 eval(nk) 
  complex*16 Peff(nk,nk)    !effective pairing interaction between patches in the momentum space

  !purpose
  !  find the eigen gap functions in the BCS theory, using Peff as the effective pairing interaction


  integer i,j

  do j=1,nk; do i=1,nk; Peff(i,j)=-Peff(i,j)*sqrt(pdos(i)*pdos(j)); end do; end do
  do i=1,nk; do j=1,nk; if(abs(Peff(i,j)-conjg(Peff(j,i)))>1.e-6)print*,'hermitian error @ BCS_k'; end do; end do

  call ZHEIGEN(nk,Peff,eval);  eval=-eval

  do j=1,nk; do i=1,nk; Peff(i,j)=Peff(i,j)/sqrt(pdos(i)); end do; end do

  return
end subroutine BCS_k

subroutine BCSMFT(nk,pdos,Peff,cutoff)
  implicit none
  integer nk
  real*8 pdos(nk),cutoff
  complex*16 Peff(nk,nk)

  complex*16 work(nk,nk),del(nk),delnew(nk),delta
  real*8 factor(nk),eval(nk),lambda
  integer iseed,ik,jk,icount
  real*8 error,dt,epsilon,temperature,Pmax
  real*8 ran,arcsinh

  print*,'running BCS mean field theory in k space @ cutoff ', cutoff

  !judge the coupling constant by linear BCS
  work=Peff; call BCS_k(nk,pdos,work,eval); lambda=eval(1)
  
  !rescale Peff so that coupling constant ~ 4
  Peff = Peff * 4. / lambda
    
  !initial gap, either random or from eigenvectors
  !iseed=98435
  !do ik=1,nk; del(ik)=del(ik) + cmplx( ran(iseed)-0.5, ran(iseed)-0.5 ); end do
  del = work(:,1) + cmplx(0,1) * work(:,2); del = del *0.01

  error=1; epsilon=1.e-6; dt=0.5*lambda/4; if(dt>1)dt=1; temperature=1.e-6

  icount=0
  do while (error>epsilon)
     factor=0
	 do ik=1,nk; if(abs(del(ik))==0)cycle
	    factor(ik)=arcsinh( cutoff/abs(del(ik)) )  !log( cutoff / sqrt( temperature**2+abs(del(ik))**2 ) )  !arcsinh( cutoff/abs(del(ik)) )
	 end do
	 do ik=1,nk
        delnew(ik)=sum( Peff(ik,:)*pdos(:)*del(:)*factor(:) )
     end do
	 error=0
	 do ik=1,nk
		delta=delnew(ik)-del(ik)
		error=max(error,abs(delta))
		del(ik)=del(ik)+dt*delta
	 end do
	 icount=icount+1
	 if(icount>4000)then
	   print*,'mean field error=',error; exit
	 end if
  end do
  !recover Peff
  Peff = Peff * lambda / 4

  open(10,file='mftgapk.dat')
  write(10,100)del; close(10)
100 format(1x,2f20.10)

  print*,'BCS mean field theory results output to file.'
  return
end subroutine BCSMFT

function arcsinh(x)
  implicit none
  real*8 x,arcsinh

  !purpose: calculate the inverse-sinh function

  real*8 xx

  xx=x/sqrt(x*x+1)
  arcsinh=0.5*log( (1+xx)/(1-xx) )

  !if(abs(x-sinh(arcsinh))>1.e-6)pause ' error @ arcsinh func'

  return
end function arcsinh


subroutine getpdos_k(model)
  use standard_derived_types
  implicit none

  type (modelblock) :: model
  type (mfblock) :: mf

  integer nkf,ik,iband,norb
  real*8 kv(2),dk,epsilon,work
  complex*16 hk(model%norb,model%norb)
  real*8 ek(model%norb),ekk(model%norb)
  real*8 Q,pi

  norb=model%norb;  nkf=model%nkf
  dk=1.e-5; epsilon=1.e-8
  allocate(model%pdos(nkf))
  open(15,file='pdos.dat')
  pi=asin(1.d0)*2
  do ik=1,nkf; iband=model%band(ik)
     kv=model%kf(:,ik);  call gethk(norb,kv,Q,hk,model,mf,.false.); call ZHEIGEN(norb,hk,ek)
	 kv=model%kf(:,ik)+(/dk,0.d0/); call gethk(norb,kv,Q,hk,model,mf,.false.); call ZHEIGEN(norb,hk,ekk)
	 work=( (ekk(iband)-ek(iband))/dk )**2
     kv=model%kf(:,ik)+(/0.d0,dk/); call gethk(norb,kv,Q,hk,model,mf,.false.); call ZHEIGEN(norb,hk,ekk)
	 work = work + ( (ekk(iband)-ek(iband))/dk )**2
	 model%pdos(ik) = model%lf(ik) / ( sqrt(work) + epsilon ) /pi
	 write(15,*)model%pdos(ik)
  end do
  close(15)

  return
end subroutine getpdos_k

subroutine meshfermisurface(model) 
  use workplace, only : appendkmesh,appendfs
  use standard_derived_types 
  implicit none

  type (modelblock) :: model
	type (mfblock) :: mf

  integer norb,nband
  type (bandmesh), target :: bands(model%nband)
  type (bandmesh), pointer :: band
  type (patchmesh), pointer :: patch

  integer iband,npatch,ipatch,activeband,nkf,ik,nk,ig,nkir
  character*2 pocket
  real*8 a(3),b(3),c(3),ka(3),kb(3),kc(3),kg(2),k(2)
  real*8 smear,pi
  real*8 filling
  complex*16 hk(model%norb,model%norb)
  real*8 eval(model%norb),pdos
  logical reducible

  interface 
    subroutine meshpocket(band,model)
      use standard_derived_types
      type(bandmesh),target::band
      type(modelblock),target::model
    end subroutine
  end interface
  
  if(.true.)then
    if(appendfs)then
      nkf=sum(model%npatch(:))
	    model%nkf=nkf
	    allocate(model%kf(2,nkf),model%band(nkf),model%lf(nkf),model%pdos(nkf))
	    open(10,file='fs.dat')
	    do ik=1,nkf
        read(10,200)model%kf(:,ik),model%band(ik),model%pdos(ik)
      end do
	    close(10)
	    return
    end if
  end if

  norb=model%norb; nband=model%nband
  a=model%a; b=model%b; c=model%c
  if(abs(sum(a(1:2)*b(1:2)))<1.e-6)then
    !for all square lattices 
    a=(/1,0,0/); b=(/0,1,0/); c=(/0,0,1/)
  end if
  call dualvector(b,c,a,ka); ka=ka*2
  call dualvector(c,a,b,kb); kb=kb*2
  call dualvector(a,b,c,kc); kc=kc*2

  do iband=1,nband   
	 print*,'configuring band',iband,' ...'
	 band=>bands(iband)
     band%a=model%a; band%b=model%b; band%c=model%c
	 band%ka=ka;     band%kb=kb;     band%kc=kc
	 band%npatch=model%npatch(iband)
	 band%pocket=model%pocket(iband)
	 band%nativeband=model%nativeband(iband)
	 band%mushift=model%mushift(iband)
     call meshpocket(band,model)
	 print*,'nk on this band=', sum(band%patch(:)%nk)
  end do
  
  nkf=sum(bands(:)%npatch); model%nkf=nkf; pi=asin(1.d0)*2
  allocate(model%kf(2,nkf),model%band(nkf),model%lf(nkf),model%pdos(nkf))
  ik=0; smear=1.e-3  !smear is a smearing factor to calculate pdos on fermi surface
  do iband=1,nband; band=>bands(iband)
     do ipatch=1,band%npatch; patch=>band%patch(ipatch); ik=ik+1
	    model%kf(:,ik)=band%patch(ipatch)%kf
		model%band(ik)=band%nativeband
		model%lf(ik)=band%patch(ipatch)%lf
		!model%pdos(ik)=pdos(norb,iband,model%kf(:,ik),model%lf(ik),model)
		model%pdos(ik)=sum( patch%dS*smear/(smear**2+patch%ek**2) )/pi
	 end do
  end do
  
  open(7,file='fs.dat'); open(8,file='fschar.dat')
  do ik=1,nkf
     write(7,200)model%kf(:,ik),model%band(ik),model%pdos(ik)
     call gethk(model%norb,model%kf(:,ik),0.d0,hk,model,mf,.false.)
	 call ZHEIGEN(model%norb,hk,eval)
	 write(8,100)abs( hk(:, model%band(ik) ) )**2
  end do
  close(7); close(8)

	!open(7,file='fs.dat')
	!do while(.not.eof(7))
	!  read(7,200)k(2),iband,pdos
!		call gethk(norb,k,0.d0,hk,model,mf

200 format(1x,2f20.6,i6,f20.6)
100 format(1x,f20.6)

  return
end subroutine meshfermisurface


function pdos(norb,iband,kf,lf,model)
  use standard_derived_types
  implicit none

  integer norb,iband
  real*8 kf(2),lf,pdos
  type (modelblock) :: model
  type (mfblock) :: mf

  real*8 dk,kv(2),ek(norb),ekk(norb),work,epsilon
  complex*16 hk(norb,norb)

  work=0
  epsilon=1.e-3
  dk=1.e-8
  call gethk(norb,kf,0d0,hk,model,mf,.false.); call ZHEIGEN(norb,hk,ek)
  kv=kf+(/dk,0.d0/); call gethk(norb,kv,0d0,hk,model,mf,.false.); call ZHEIGEN(norb,hk,ekk)
  work=( (ekk(iband)-ek(iband))/dk )**2
  kv=kf+(/0.d0,dk/); call gethk(norb,kv,0d0,hk,model,mf,.false.); call ZHEIGEN(norb,hk,ekk)
  work = work + ( (ekk(iband)-ek(iband))/dk )**2
  pdos = lf / ( sqrt(work) + epsilon ) 
  
  return
end function

subroutine meshpocket(band,model)
  use standard_derived_types
  implicit none
  
  type (bandmesh), target :: band
  type (modelblock) :: model
  type (patchmesh), pointer :: patch
  
  character*2 pocket
  integer npatch,norb,nativeband,ipatch
  real*8 Stot,mushift
  
  interface
    subroutine Gpocket(npatch,band)
        use standard_derived_types
	integer npatch
	type(bandmesh),target :: band
    end subroutine
    subroutine Mpocket(npatch,band)
        use standard_derived_types
	integer npatch
	type(bandmesh),target :: band
    end subroutine   
    subroutine GMpocket(npatch,band)
        use standard_derived_types
	integer npatch
	type(bandmesh),target :: band
    end subroutine   
    subroutine XYpocket(npatch,band)
        use standard_derived_types
	integer npatch
	type(bandmesh),target :: band
    end subroutine   
    subroutine XMpocket(npatch,band)
        use standard_derived_types
	integer npatch
	type(bandmesh),target :: band
    end subroutine   
    subroutine MMpocket(npatch,band)
        use standard_derived_types
	integer npatch
	type(bandmesh),target :: band
    end subroutine   
    subroutine HXpocket(npatch,band)
        use standard_derived_types
	integer npatch
	type(bandmesh),target :: band
    end subroutine 
    subroutine GKpocket(npatch,band)
        use standard_derived_types
	integer npatch
	type(bandmesh),target :: band
    end subroutine      
    subroutine GKpocket_(npatch,band)
        use standard_derived_types
	integer npatch
	type(bandmesh),target :: band
    end subroutine   
  end interface  
  
  npatch=band%npatch
  pocket=band%pocket
  nativeband=band%nativeband
  mushift=band%mushift

  if(npatch==0)return        !a band is ignored if npatch=0 in this band

  select case (pocket)
    case ('G_')
	  if(mod(npatch,4)/=0)stop 'invalid npatch'
	  call Gpocket(npatch,band)
	case ('M_')
	  if(mod(npatch,4)/=0)stop 'invalid npatch'
	  call Mpocket(npatch,band)
	case ('GM')
	  if(mod(npatch,8)/=0)stop 'invalid npatch'
	  call GMpocket(npatch,band)
    case ('XY')
	  if(mod(npatch,8)/=0)stop 'invalid npatch'
	  call XYpocket(npatch,band)
    case ('XM')
	  if(mod(npatch,12)/=0)stop 'invalid npatch'
	  call XMpocket(npatch,band)
    case ('MM')
	  !if(mod(npatch,12)/=0)stop 'invalid npatch'
	  call MMpocket(npatch,band)
	case ('HX')
	  call HXpocket(npatch,band)
	case ('GK')
	  call GKpocket_(npatch,band)
	case default; stop 'unknown pocket'
  end select
  do ipatch=1,npatch; patch=>band%patch(ipatch)
	 call meshpatch(patch,nativeband,model)
	 call eigenpatch(patch,nativeband,mushift,model)
  end do

  Stot=0
  do ipatch=1,npatch; Stot=Stot+sum(band%patch(ipatch)%dS); end do
  do ipatch=1,npatch; band%patch(ipatch)%dS=band%patch(ipatch)%dS/Stot; end do

  return
end subroutine meshpocket

subroutine meshpatch(patch,nativeband,model)
  use standard_derived_types
  implicit none

  integer nativeband
  type (patchmesh) :: patch
  type (modelblock) :: model
  type (meshconfig) :: line(2)

  integer ntri,itri,ik,nk
  real*8, dimension (2) :: kc,km,dk,kstart,kend

  integer ng,ip,ig
  real*8 wmax,wmin,S
  real*8 area,dL

  ntri=patch%ntri

  ng=12; wmax=9; wmin=1.e-5
  !ng=8; wmax=3; wmin=1.e-4
  !ng=10; wmax=10; wmin=1.e-5

  nk=0
  do itri=1,ntri
     kc=patch%kc(:,itri); km=patch%km
     call linemesh(ng,wmax,wmin,kc,km,line(itri),epfunc)
	 nk=nk+line(itri)%npoints
  end do
  
  patch%nk=nk; allocate(patch%kv(2,nk),patch%dS(nk),patch%dl(nk))

  ik=0
  do itri=1,ntri
     kc=patch%kc(:,itri); km=patch%km; S=area(patch%kb(:,1+(itri-1)*3:itri*3))
	 dL=sqrt( sum( ( patch%kb(:,2+(itri-1)*3) - patch%kb(:,itri*3) )**2 ) )
     do ig=1,ng; do ip=1,line(itri)%g(ig)%np        
	    if(line(itri)%g(ig)%info(ip)/=0)cycle
	    ik=ik+1
		patch%kv(:,ik)=kc+(km-kc)*line(itri)%g(ig)%p(1,ip)
		patch%dS(ik)=line(itri)%g(ig)%p(1,ip)*line(itri)%g(ig)%dp*2*S
		patch%dl(ik)=line(itri)%g(ig)%p(1,ip)*dL
     end do; end do
  end do

  return
      contains
      function epfunc(p)
        implicit none
        real*8 epfunc,p,Q,findek
		real*8, dimension (2) :: kv
		kv=kc+(km-kc)*p
        epfunc=findek(nativeband,kv,Q,model%norb,model)       
        return
      end function epfunc

end subroutine meshpatch

  
subroutine eigenpatch(patch,nativeband,mushift,model)
  use standard_derived_types
  implicit none

  integer nativeband
  type (patchmesh) :: patch
  type (modelblock) :: model
  real*8 mushift

  integer norb
  integer nk,ik,itri
  real*8 ek,ekmin,Q
  real*8, dimension (2) :: kc,km,dk,kv

  real*8 findek,fermi

  norb=model%norb

  nk=patch%nk
  allocate(patch%ek(nk),patch%fk(nk))

  ekmin=1.e10
  do ik=1,nk
     kv=patch%kv(:,ik)
     ek=findek(nativeband,kv,Q,norb,model)
	 patch%ek(ik)=ek
	 patch%fk(ik)=fermi(ek,1.e-3)
	 if(abs(ek-mushift)<ekmin)then
	   ekmin=abs(ek-mushift)
	   patch%kf=kv; patch%lf=patch%dl(ik)
	 end if
  end do
  return
end subroutine eigenpatch



subroutine Gpocket(npatch,band)
  use standard_derived_types
  implicit none

  integer npatch
  type (bandmesh) ,target :: band

  real*8, dimension (2) :: kc1,ki,kc2,kf 
  real*8, allocatable, dimension (:,:) :: p,dp

  integer n4,i,ig,icount
  type (patchmesh), pointer :: patch
  real*8 angle

  !information in a patch
  !integer ntri
  !real*8  kc(2,2),km(2),kb(2,6)

  if(mod(npatch,4)/=0)stop 'invalid npatch'
  band%npatch=npatch; allocate(band%patch(npatch))

  n4=npatch/4
  allocate(p(2,n4),dp(2,n4))

  kc1=0; kc2=(/1,1/); ki=(/1,0/); kf=(/0,1/)
  call tiltpatch(ki,kf,n4,p,dp)
  do i=1,n4; patch=>band%patch(i)
     patch%ntri=2;	 patch%kc(:,1)=kc1
	 patch%kc(:,2)=kc2; patch%km=p(:,i) 
	 patch%kb(:,1)=kc1; patch%kb(:,2)=p(:,i)+dp(:,i)/2; patch%kb(:,3)=p(:,i)-dp(:,i)/2
	 patch%kb(:,4)=kc2; patch%kb(:,5)=p(:,i)+dp(:,i)/2; patch%kb(:,6)=p(:,i)-dp(:,i)/2
  end do 
  icount=n4
  
  do ig=1,3
     call rotate(kc1); call rotate(kc2)
     do i=1,n4; call rotate(p(:,i)); call rotate(dp(:,i)); end do
     do i=1,n4; icount=icount+1; patch=>band%patch(icount)
        patch%ntri=2;	 patch%kc(:,1)=kc1
	    patch%kc(:,2)=kc2; patch%km=p(:,i) 
	    patch%kb(:,1)=kc1; patch%kb(:,2)=p(:,i)+dp(:,i)/2; patch%kb(:,3)=p(:,i)-dp(:,i)/2
	    patch%kb(:,4)=kc2; patch%kb(:,5)=p(:,i)+dp(:,i)/2; patch%kb(:,6)=p(:,i)-dp(:,i)/2
     end do
  end do

  do i=1,npatch; patch=>band%patch(i)
     patch%angle=angle(patch%km)
  end do

  return
end subroutine Gpocket

subroutine Mpocket(npatch,band)
  use standard_derived_types
  implicit none

  integer npatch
  type (bandmesh), target :: band

  integer ipatch,i
  real*8 shift(2)
  type (patchmesh), pointer :: patch
  real*8 angle
  
  interface
    subroutine Gpocket(npatch,band)
      use standard_derived_types
      integer npatch
      type(bandmesh),target :: band
    end subroutine
  end interface

  if(mod(npatch,4)/=0)stop 'invalid npatch @ M pocket'

  call Gpocket(npatch,band)
  
  do ipatch=1,npatch; patch=>band%patch(ipatch)
	 shift=0
	 where(patch%km+(/1,1/)>1)shift=-2; where(patch%km+(/1,1/)<-1)shift=2
	 patch%km=patch%km+(/1,1/)+shift
     do i=1,2; patch%kc(:,i)=patch%kc(:,i)+(/1,1/)+shift; end do
	 do i=1,6; patch%kb(:,i)=patch%kb(:,i)+(/1,1/)+shift; end do
  end do

  do i=1,npatch; patch=>band%patch(i)
     patch%angle=angle(patch%km-(/1,1/))
  end do

  return
end subroutine Mpocket

subroutine GMpocket(npatch,band)
  use standard_derived_types
  implicit none

  integer npatch
  type (bandmesh), target :: band
  type (patchmesh), pointer :: patch

  real*8, dimension (2) :: kc,ki,km,kf,shift
  real*8, allocatable, dimension (:,:) :: p,dp
  real*8 angle

  integer n2,n8,i,ig
  integer version

  version=2
 
  band%npatch=npatch; allocate(band%patch(npatch))
  if(mod(npatch,8)/=0)stop 'invalid npatch @ GM pocket'
  n2=npatch/2; n8=npatch/8
  allocate(p(2,n2),dp(2,n2))
  
  kc=0; ki=(/1,0/); kf=(/0,1/); km=(/0.5,0.5/)
  call radiatepatch(ki,km,kf,n8,p,dp,version)

  do ig=2,4
     call rotate(kc)
	 p(:,n8*(ig-1)+1 : n8*ig)=p(:,n8*(ig-2)+1 : n8*(ig-1) )
	 dp(:,n8*(ig-1)+1 : n8*ig)=dp(:,n8*(ig-2)+1 : n8*(ig-1) )
     do i=1,n8; call rotate( p(:,i+n8*(ig-1)) ); call rotate( dp(:,i+n8*(ig-1)) ); end do
  end do

  do i=1,n2; patch=>band%patch(i)
     patch%ntri=1
     patch%kc(:,1)=kc; patch%km=p(:,i)
	 patch%kb(:,1)=kc; patch%kb(:,2)=p(:,i)+dp(:,i)/2; patch%kb(:,3)=p(:,i)-dp(:,i)/2
	 patch%angle=angle(p(:,i)-kc)
  end do

  kc=(/1,1/)
  do i=1,n2; p(:,i)=p(:,i)+kc; end do

  do i=1,n2; patch=>band%patch(i+n2) 
     shift=0
     where(p(:,i)>1)shift=-2; where(p(:,i)<-1)shift=2
     patch%ntri=1
     patch%kc(:,1)=kc+shift; patch%km=p(:,i)+shift
	 patch%kb(:,1)=kc+shift; patch%kb(:,2)=p(:,i)+shift+dp(:,i)/2; patch%kb(:,3)=p(:,i)+shift-dp(:,i)/2
     patch%angle=angle(patch%km-kc)
  end do

  return
end subroutine GMpocket

subroutine XYpocket(npatch,band)
  use standard_derived_types
  implicit none

  integer npatch
  type (bandmesh), target :: band

  integer ipatch,i
  real*8 shift(2)
  type (patchmesh), pointer :: patch
  real*8 angle

  interface
    subroutine GMpocket(npatch,band)
        use standard_derived_types
	integer npatch
	type(bandmesh),target :: band
    end subroutine
  end interface

  if(mod(npatch,8)/=0)stop 'invalid npatch @ XYpocket'

  call GMpocket(npatch,band)
  
  do ipatch=1,npatch; patch=>band%patch(ipatch)
     shift=0
	 where(patch%km+(/1,0/)>1)shift=-2
	 where(patch%km+(/1,0/)<-1)shift=2
     patch%km=patch%km+(/1,0/)+shift
     do i=1,2; patch%kc(:,i)=patch%kc(:,i)+(/1,0/)+shift; end do
	 do i=1,6; patch%kb(:,i)=patch%kb(:,i)+(/1,0/)+shift; end do
	 patch%angle=angle(patch%km-patch%kc(:,1))   !only kc(:,1) is used in a 1-triangle patch
  end do

  return
end subroutine XYpocket

subroutine XMpocket(npatch,band)
  use standard_derived_types
  implicit none

  integer npatch,n1,n2
  type (bandmesh), target :: band,b1,b2

  integer ipatch,i
  real*8 shift(2),v(2)
  type (patchmesh), pointer :: patch,p1,p2
  real*8 angle,twopi

  interface
    subroutine GMpocket(npatch,band)
        use standard_derived_types
	integer npatch
	type(bandmesh),target :: band
    end subroutine
    subroutine Mpocket(npatch,band)
        use standard_derived_types
	integer npatch
	type(bandmesh),target :: band
    end subroutine
  end interface

  if(mod(npatch,12)/=0)stop 'invalid npatch @ XMpocket'

  n1=npatch*2/3; n2=npatch-n1; twopi=asin(1.d0)*2
  
  call GMpocket(n1,b1); shift=(/1,1/)
  do ipatch=1,n1; patch=>b1%patch(ipatch)
     v=patch%km; call rotate45(v); patch%km=v/sqrt(2.d0)+shift
	 do i=1,2; v=patch%kc(:,i); call rotate45(v); patch%kc(:,i)=v/sqrt(2.d0)+shift; end do
     do i=1,6; v=patch%kb(:,i); call rotate45(v); patch%kb(:,i)=v/sqrt(2.d0)+shift; end do
     patch%angle=mod(patch%angle+twopi/8,twopi)
  end do

  do ipatch=1,n1; patch=>b1%patch(ipatch)
     shift=0
	 where(patch%km>1)shift=-2
	 where(patch%km<-1)shift=2
     patch%km=patch%km+shift
     do i=1,2; patch%kc(:,i)=patch%kc(:,i)+shift; end do
	 do i=1,6; patch%kb(:,i)=patch%kb(:,i)+shift; end do
  end do

  call Mpocket(n2,b2)  
  do ipatch=1,n2; patch=>b2%patch(ipatch)
	 v=patch%km; call rotate45(v); patch%km=v/sqrt(2.d0)
     do i=1,2; v=patch%kc(:,i); call rotate45(v); patch%kc(:,i)=v/sqrt(2.d0); end do
	 do i=1,6; v=patch%kb(:,i); call rotate45(v); patch%kb(:,i)=v/sqrt(2.d0); end do
	 patch%angle=mod(patch%angle+twopi/8,twopi)
  end do

  band%npatch=npatch
  band%pocket='XM'
  allocate(band%patch(npatch))
  
  do ipatch=1,n1; patch=>band%patch(ipatch); p1=>b1%patch(ipatch)
     patch%ntri=p1%ntri; patch%km=p1%km
	 patch%kc=p1%kc; patch%kb=p1%kb; patch%angle=p1%angle
  end do

  do ipatch=1,n2; patch=>band%patch(ipatch+n1); p2=>b2%patch(ipatch)
     patch%ntri=p2%ntri; patch%km=p2%km
	 patch%kc=p2%kc; patch%kb=p2%kb; patch%angle=p2%angle
  end do

  return
end subroutine XMpocket


subroutine MMpocket(npatch,band)
  use standard_derived_types
  implicit none

  integer npatch,n1,n2
  type (bandmesh), target :: band,b1,b2

  integer ipatch,i
  real*8 shift(2),v(2)
  type (patchmesh), pointer :: patch,p1,p2
  real*8 angle,twopi

  interface
    subroutine GMpocket(npatch,band)
        use standard_derived_types
	integer npatch
	type(bandmesh),target :: band
    end subroutine
  end interface


  n1=npatch/2; n2=npatch-n1; twopi=asin(1.d0)*2
  if(mod(n1,4)/=0.or.mod(n2,4)/=0)stop 'invalid GM pocket @ MMpocket'
 
    call GMpocket(n1,b1); shift=(/1,1/)
  do ipatch=1,n1; patch=>b1%patch(ipatch)
     v=patch%km; call rotate45(v); patch%km=v/sqrt(2.d0)+shift
	 do i=1,2; v=patch%kc(:,i); call rotate45(v); patch%kc(:,i)=v/sqrt(2.d0)+shift; end do
     do i=1,6; v=patch%kb(:,i); call rotate45(v); patch%kb(:,i)=v/sqrt(2.d0)+shift; end do
     patch%angle=mod(patch%angle+twopi/8,twopi)
  end do

  do ipatch=1,n1; patch=>b1%patch(ipatch)
     shift=0
	 where(patch%km>1)shift=-2
	 where(patch%km<-1)shift=2
     patch%km=patch%km+shift
     do i=1,2; patch%kc(:,i)=patch%kc(:,i)+shift; end do
	 do i=1,6; patch%kb(:,i)=patch%kb(:,i)+shift; end do
  end do

  call GMpocket(n2,b2)  
  do ipatch=1,n2; patch=>b2%patch(ipatch)
	 v=patch%km; call rotate45(v); patch%km=v/sqrt(2.d0)
     do i=1,2; v=patch%kc(:,i); call rotate45(v); patch%kc(:,i)=v/sqrt(2.d0); end do
	 do i=1,6; v=patch%kb(:,i); call rotate45(v); patch%kb(:,i)=v/sqrt(2.d0); end do
	 patch%angle=mod(patch%angle+twopi/8,twopi)
  end do

  band%npatch=npatch
  band%pocket='MM'
  allocate(band%patch(npatch))
  
  do ipatch=1,n1; patch=>band%patch(ipatch); p1=>b1%patch(ipatch)
     patch%ntri=p1%ntri; patch%km=p1%km
	 patch%kc=p1%kc; patch%kb=p1%kb; patch%angle=p1%angle
  end do

  do ipatch=1,n2; patch=>band%patch(ipatch+n1); p2=>b2%patch(ipatch)
     patch%ntri=p2%ntri; patch%km=p2%km
	 patch%kc=p2%kc; patch%kb=p2%kb; patch%angle=p2%angle
  end do

  return
end subroutine MMpocket


subroutine HXpocket(npatch,band)
  use standard_derived_types
  implicit none

  integer npatch
  type (bandmesh), target :: band

  real*8, dimension (2) :: kc1,ki,kc2,kf 
  real*8, allocatable, dimension (:,:) :: p,dp

  integer n6,i,ig,icount
  type (patchmesh), pointer :: patch
  real*8 angle

  !information in a patch
  !integer ntri
  !real*8  kc(2,2),km(2),kb(2,6)

  if(mod(npatch,6)/=0)stop 'invalid npatch'
  band%npatch=npatch; allocate(band%patch(npatch))

  n6=npatch/6
  allocate(p(2,n6),dp(2,n6))

  kc1=0; kc2=(/2./3,2./sqrt(3.)/); ki=(/1.,1./sqrt(3.)/); kf=(/0.,2./sqrt(3.)/)
  call tiltpatch(ki,kf,n6,p,dp)
  do i=1,n6; patch=>band%patch(i)
     patch%ntri=2;	 patch%kc(:,1)=kc1
	 patch%kc(:,2)=kc2; patch%km=p(:,i) 
	 patch%kb(:,1)=kc1; patch%kb(:,2)=p(:,i)+dp(:,i)/2; patch%kb(:,3)=p(:,i)-dp(:,i)/2
	 patch%kb(:,4)=kc2; patch%kb(:,5)=p(:,i)+dp(:,i)/2; patch%kb(:,6)=p(:,i)-dp(:,i)/2
  end do 
  icount=n6
  
  do ig=1,5
     call rotate60(kc1); call rotate60(kc2)
     do i=1,n6; call rotate60(p(:,i)); call rotate60(dp(:,i)); end do
     do i=1,n6; icount=icount+1; patch=>band%patch(icount)
        patch%ntri=2;	 patch%kc(:,1)=kc1
	    patch%kc(:,2)=kc2; patch%km=p(:,i) 
	    patch%kb(:,1)=kc1; patch%kb(:,2)=p(:,i)+dp(:,i)/2; patch%kb(:,3)=p(:,i)-dp(:,i)/2
	    patch%kb(:,4)=kc2; patch%kb(:,5)=p(:,i)+dp(:,i)/2; patch%kb(:,6)=p(:,i)-dp(:,i)/2
     end do
  end do

  do i=1,npatch; patch=>band%patch(i)
     patch%angle=angle(patch%km)
  end do

  return
end subroutine HXpocket

subroutine GKpocket(npatch,band)
  use standard_derived_types
  implicit none

  integer npatch
  type (bandmesh), target :: band

  real*8, dimension (2) :: kc,ki,kf,km 
  real*8, allocatable, dimension (:,:) :: p,dp

  integer n6,n12,i,k,ig,nL,nU
  type (patchmesh), pointer :: patch
  real*8 angle

  integer version

  !information in a patch
  !integer ntri
  !real*8  kc(2,2),km(2),kb(2,6)
  
  version=2

  if(mod(npatch,12)/=0)stop 'invalid npatch @ GKpocket'
  band%npatch=npatch; allocate(band%patch(npatch))

  n6=npatch/6; n12=npatch/12
  allocate(p(2,n12),dp(2,n12))

  kc=0; ki=(/1.,-1./sqrt(3.)/); kf=(/1.,1./sqrt(3.)/); km=(ki+kf)/2
  call radiatepatch(ki,km,kf,n12,p,dp,version)
  do i=1,n12; patch=>band%patch(i)
     patch%ntri=1;	 patch%kc(:,1)=kc; patch%km=p(:,i) 
	 patch%kb(:,1)=kc; patch%kb(:,2)=p(:,i)+dp(:,i)/2; patch%kb(:,3)=p(:,i)-dp(:,i)/2
  end do 

  kc=(/4./3,0./); ki=(/1.,1./sqrt(3.)/); kf=(/1.,-1./sqrt(3.)/); km=(ki+kf)/2
  call radiatepatch(ki,km,kf,n12,p,dp,version)
  do i=1,n12; patch=>band%patch(i+n12)
     patch%ntri=1;	patch%kc(:,1)=kc; patch%km=p(:,i) 
	 patch%kb(:,1)=kc; patch%kb(:,2)=p(:,i)+dp(:,i)/2; patch%kb(:,3)=p(:,i)-dp(:,i)/2
  end do 
  
  do ig=1,5; nL=1+n6*ig; nU=n6*(ig+1)
     band%patch(nL:nU)=band%patch(nL-n6:nU-n6)
     do i=1,n6; patch=>band%patch(nL+i-1)
        call rotate60(patch%kc(:,1)); call rotate60(patch%km)
	 	do k=1,3; call rotate60(patch%kb(:,k)); end do
     end do
  end do

  do i=1,npatch; patch=>band%patch(i)
     patch%angle=angle(patch%km)
  end do

  return
end subroutine GKpocket


subroutine GKpocket_(npatch,band)
  use standard_derived_types
  implicit none

  integer npatch
  type (bandmesh), target :: band

  real*8, dimension (2) :: kc,ki,kf,km 
  real*8, allocatable, dimension (:,:) :: p,dp

  integer n2,n6,n12,i,k,ig,nL,nU
  type (patchmesh), pointer :: patch
  real*8 angle

  integer version

  !information in a patch
  !integer ntri
  !real*8  kc(2,2),km(2),kb(2,6)
  
  version=2

  if(mod(npatch,12)/=0)stop 'invalid npatch @ GKpocket'
  band%npatch=npatch; allocate(band%patch(npatch))

  n2=npatch/2; n6=npatch/6; n12=npatch/12
  allocate(p(2,n12),dp(2,n12))

  !G pocket
  kc=0; ki=(/1.,-1./sqrt(3.)/); kf=(/1.,1./sqrt(3.)/); km=(ki+kf)/2
  call radiatepatch(ki,km,kf,n12,p,dp,version)
  do i=1,n12; patch=>band%patch(i)
     patch%ntri=1;	 patch%kc(:,1)=kc; patch%km=p(:,i) 
	 patch%kb(:,1)=kc; patch%kb(:,2)=p(:,i)+dp(:,i)/2; patch%kb(:,3)=p(:,i)-dp(:,i)/2
  end do 

  do ig=1,5; nL=1+n12*ig; nU=n12*(ig+1)
     band%patch(nL:nU)=band%patch(nL-n12:nU-n12)
     do i=1,n12; patch=>band%patch(nL+i-1)
        call rotate60(patch%kc(:,1)); call rotate60(patch%km)
	 	do k=1,3; call rotate60(patch%kb(:,k)); end do
     end do
  end do

  !K pockets
  kc=(/4./3,0./); ki=(/1.,1./sqrt(3.)/); kf=(/1.,-1./sqrt(3.)/); km=(ki+kf)/2
  call radiatepatch(ki,km,kf,n12,p,dp,version)
  do i=1,n12; patch=>band%patch(i+n2)
     patch%ntri=1;	patch%kc(:,1)=kc; patch%km=p(:,i) 
	 patch%kb(:,1)=kc; patch%kb(:,2)=p(:,i)+dp(:,i)/2; patch%kb(:,3)=p(:,i)-dp(:,i)/2
  end do 
  
  do ig=1,5; nL=1+n12*ig+n2; nU=n12*(ig+1)+n2
     band%patch(nL:nU)=band%patch(nL-n12:nU-n12)
     do i=1,n12; patch=>band%patch(nL+i-1)
        call rotate60(patch%kc(:,1)); call rotate60(patch%km)
	 	do k=1,3; call rotate60(patch%kb(:,k)); end do
     end do
  end do

  !angle of patches
  do i=1,npatch; patch=>band%patch(i)
     patch%angle=angle(patch%km)
  end do

  return
end subroutine GKpocket_
  
subroutine rotate(r)
  implicit none
  real*8 r(2),x
  x=r(1)
  r(1)=-r(2)
  r(2)=x
  return
end subroutine rotate


subroutine RadiatePatch(ki,km,kf,np,p,dp,version)
  implicit none
  integer np
  real*8, dimension (2) :: ki,km,kf
  real*8, dimension (2,np) :: p,dp
  integer version

  integer nh,n
  real*8 dk(2)

  select case (version)
    case (2)
      dk=(kf-ki)/np
      do n=1,np
         p(:,n)=ki+dk*(n-0.5)
	     dp(:,n)=dk
      end do
    case (1)
      if(mod(np,2)/=0)stop 'np must be even @ radiatepatch'
      nh=np/2

      dk=(km-ki)/nh
      do n=1,nh
         p(:,n)=ki+dk*(n-0.5)
	     dp(:,n)=dk
      end do

      dk=(kf-km)/nh
      do n=1,nh
         p(:,n+nh)=km+dk*(n-0.5)
	     dp(:,n)=dk
      end do
	case default; stop 'unknown version of radiatepatch'
  end select
  return
end subroutine RadiatePatch

subroutine TiltPatch(ki,kf,np,p,dp)
  implicit none
  integer np
  real*8, dimension (2) :: ki,kf
  real*8, dimension (2,np) :: p,dp

  integer n
  real*8 dk(2)

  dk=(kf-ki)/np
  do n=1,np
     p(:,n)=ki+dk*(n-0.5)
	 dp(:,n)=dk
  end do

  return
end subroutine TiltPatch

function withinpatch(vk,vb)
  implicit none
  logical withinpatch
  real*8 vk(2),vb(2,3)
 
  real*8, dimension(2) :: c,m,n,cm,vm
  integer i,ii

  !purpose: check whether a point vk is within a triangle specified the three points in vb

  !internal variables:
  ! c: center of mass of the triangle
  ! m: mid point of a side
  ! n: normal vector of a side, directed outward
  ! cm: vector from c to m
  ! vm: vector from vk to m
  
  !method:
  ! let n be a vector orthogonal to m-side. the direction of n is fixed by requiring n*cm>0.
  ! if vm*n>0 for all sides then v is inside the triangle
  ! if vm*n=0 for m-side, then v is on this side
  ! if vm*n<0 for any side, then v is out of the triangle
 
  c=0; do i=1,3; c=c+vb(:,i); end do; c=c/3
  
  withinpatch=.true.
  do i=1,3; ii=mod(i,3)+1
     m=vb(:,i)-vb(:,ii); n(1)=m(2); n(2)=-m(1)  !construct normal vector
     m=vb(:,i)+vb(:,ii); m=m/2                  !mid point of the side
	 cm=m-c; vm=m-vk                            !vectors from c to m, and vk to m
	 if(sum(n*cm)<0)n=-n                        !fix the direction of n (outward)
	 if(sum(vm*n)<=0)then
	   withinpatch=.false.; return              !if projection of vm on n is negative, vk is outside of the triangle
     end if
  end do

  return
end function withinpatch


function angle(r)
  implicit none
  real*8 angle,r(2)
  
  real*8 rr(2),pi

  rr=mod(r,2.d0); where(rr<-1)rr=rr+2; where(rr>1)rr=rr-2
  
  if(sum(abs(rr))==0)then; angle=0; return; end if

  pi=asin(1.d0)*2
  if(rr(1)==0)then
    if(rr(2)>0)then; angle=pi/2; else; angle=3*pi/2; end if
  else
    angle=atan(rr(2)/rr(1))
	if(rr(1)<0)then; angle=angle+pi; end if
  end if
 
  if(angle<0)angle=angle+pi*2
  return
end function angle
