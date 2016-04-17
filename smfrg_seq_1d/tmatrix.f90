subroutine TMatrix()
  use workplace, only : one,nambu,norb,ndim,pi
  implicit none
  
  integer is,idim,nrcut,nrv,icut,ircut,ir,iorb,iw,nw,m,n,nvs,nvm,icount
  real*8 s,u1,u2,Vm,Vs,wmax,wmin,eta,vr(2),Akw,Akw_,delta,bw,vmax,vmin,dw
  real*8,allocatable :: wcut(:),rv(:,:)
  complex*16 U(ndim),w,ZDET
  complex*16,dimension(nambu,nambu) :: T
  complex*16,allocatable :: Gr(:,:,:,:)
  complex*16 :: Grr(nambu,nambu)
  logical appendGr,etaw,logw
  character*3 fileid,digitize
  
  delta=0.05
  eta=0.01 
  wmax=0.15
  nw=60
  nrv=0
  appendGr=.true.    
  
  Vs=-0.01
  Vm=0
  !set up w-mesh:  (-2*nw:-nw-1)&(nw+1:2*nw): linear; (-nw:-1)&(1:nw): log; (0): 0
  allocate(wcut(-2*nw:2*nw))
  do iw=-2*nw,2*nw
    if(nw==0)then
      wcut(0)=0.d0
      exit
    endif
  wcut(iw)=wmax*iw/2/nw
  end do
  open(10,file='w.dat')
    write(10,'(1e16.8)')wcut
  close(10)

  !setup r-mesh
  allocate(rv(2,-nrv:nrv))
  open(10,file='rv.dat')
  do ir=-nrv,nrv
    rv(:,ir)=(/ir,0/)*1.d0
    write(10,'(2e16.8)')rv(:,ir)
  end do
  close(10)

  !real pairmode
  open(10,file='pairmode.dat')
    read(10,*)is,s
    do idim=1,ndim
    read(10,*)u1,u2
    U(idim)=cmplx(u1,u2)
  end do
  close(10)
  call BdGBand(delta,U)

  !prepare Gr
  allocate(Gr(nambu,nambu,-2*nw:2*nw,-nrv:nrv))
  if(.not.appendGr)call BdGeigentable(delta,U)
  
  icount=0
  do iw=-2*nw,2*nw
    icount=icount+1
    w=wcut(iw)+one*eta
    fileid=digitize(icount)
    open(10,file='./iron_impurity/nfasc/'//'Gr.'//fileid)
    !open(10,file='Gr.'//fileid)
	if(appendGr)then
      do ir=-nrv,nrv
        read(10,'(1e16.8)')Gr(:,:,iw,ir)
      end do
    else
      print*,'iw=',iw
      do ir=-nrv,nrv
        vr=rv(:,ir)
        call getBdGGr(w,vr,Gr(:,:,iw,ir))
        write(10,'(1e16.8)')Gr(:,:,iw,ir)
      end do
    end if
    close(10)
  end do

  !output G00
  open(10,file='G00.dat')
  do m=1,nambu
  do n=1,nambu
    write(10,'(2e16.8)')Gr(m,n,0,0)
  end do
  end do
  close(10)

  !caculate los at iw
  if(.false.)then
  iw=0; ir=0
  open(10,file='losw0.dat')
  do m=1,501; vm=10.0/500*(m-1)-5.0
  do n=1,501; vs=10.0/500*(n-1)-5.0
    call initT(vm,vs,T)
	T=T-Gr(:,:,iw,0)
	call ZINVERT(nambu,T)
    Grr=Gr(:,:,iw,0)+matmul( Gr(:,:,iw,-ir), matmul( T, Gr(:,:,iw,ir) ) )
    Akw=0
    do iorb=1,norb
      Akw=Akw-1/pi*imag(Grr(iorb,iorb))
    end do
    call initT(vm,vs,T)
	T=T-Gr(:,:,-iw,0)
	call ZINVERT(nambu,T)
	Grr=Gr(:,:,-iw,0)+matmul( Gr(:,:,-iw,-ir), matmul( T, Gr(:,:,-iw,ir) ) )
    do iorb=norb+1,nambu
      Akw=Akw-1/pi*imag(Grr(iorb,iorb))
    end do
    write(10,*)Akw
  end do
  end do
  close(10)
  end if

  !calculating LOS
  open(10,file='akw.dat')
  open(12,file='akw_.dat')
  open(11,file='T.dat')
  do iw=-2*nw,2*nw
    call initT(Vm,Vs,T)  !get V-1
    T=T-Gr(:,:,iw,0)
    call ZINVERT(nambu,T)
    write(11,'(1e16.8)')T
    do ir=-nrv,nrv
      print*,'iw,ir',iw,ir
      Grr=Gr(:,:,iw,0)+matmul( Gr(:,:,iw,-ir), matmul( T, Gr(:,:,iw,ir) ) )
      Akw=0
      Akw_=0
      do iorb=1,norb
	    Akw=Akw-1/pi*imag(Grr(iorb,iorb))
        !Akw=Akw+real(Grr(iorb,iorb))
      end do
      write(10,'(1e16.8)')Akw       !A(w) for spin up
      do iorb=norb+1,nambu
	    Akw_=Akw_-1/pi*imag(Grr(iorb,iorb))
		!Akw_=Akw_+real(Grr(iorb,iorb))
      end do
      write(12,'(1e16.8)')Akw_      !A(-w) for spin down
    end do
  end do
  close(10)
  close(11)
  close(12)

  return
end subroutine

subroutine initT(Vm,Vs,T)
  use workplace
  implicit none

  real*8 Vm,Vs
  complex*16 T(nambu,nambu)

  integer iorb
  real*8 eta

  T=0; if(Vm<1.e-8.or.Vs<1.e-8.or.(Vm-Vs)<1.e-8)eta=1.e-8
  do iorb=1,norb
    T(iorb,iorb)=Vm+Vs+eta
	if(nambu/=norb)T(iorb+norb,iorb+norb)=Vm-Vs+eta
  end do
  call ZINVERT(nambu,T)

  return
end subroutine


subroutine getBdGGr(w,vr,Gr)   !Gr0(w,ri,rj)=Gr0(w,ri-rj)
  use workplace
  implicit none

  real*8 vr(2),dk,k(2)
  complex*16 w,Gr(nambu,nambu)

  integer ik
  complex*16 Gk(nambu,nambu,nk)

  call BdGGk(w,Gk)
  Gr=0
  do ik=1,nk
    k=kv(1:2,ik)
	dk=kv(3,ik)
    Gr=Gr+Gk(:,:,ik)*exp(one*pi*sum(k*vr))*dk
  end do

  return
end subroutine

subroutine BdGGk(w,Gk)
  use workplace
  implicit none
  
  complex*16 w,Gk(nambu,nambu,nk)
  
  integer ik,iorb,jorb,m

  Gk=0
  do ik=1,nk
	do iorb=1,nambu
	do jorb=1,nambu
      do m=1,nambu
	    Gk(iorb,jorb,ik)=Gk(iorb,jorb,ik)+BdGAk(iorb,m,ik)*conjg(BdGAk(jorb,m,ik))/(w-BdGek(m,ik))
	  end do
	end do
	end do
  end do

  return
end subroutine

subroutine getBdGGr_(w,vr,Gr)   !Gr0(w,ri,rj)=Gr0(w,ri-rj)
  use workplace
  implicit none

  real*8 vr(2)
  complex*16 w,Gr(nambu,nambu)

  integer ik,ikk,ikir,ig
  real*8 dk,k(2),kg(2)
  complex*16 Gk(nambu,nambu,nk)
  logical reducible

  call BdGGk_(w,Gk)
  Gr=0
  ikir=0
  ikk=0
  do ik=1,nk
    k=kv(1:2,ik)
	dk=kv(3,ik)
	if(reducible(k,group))cycle
	ikir=ikir+1
	ikk=ikk+1
	Gr=Gr+Gk(:,:,ikir)*exp(one*pi*sum(k*vr))*dk
    do ig=2,ng
	  call groupaction(ig,k,kg,group)
	  ikk=ikk+1
	  Gr=Gr+Gk(:,:,ikk)*exp(one*pi*sum(kg*vr))*dk
	end do
	ikir=ikir+ng-1
  end do
  print*,ikk,ikir
  return
end subroutine

subroutine BdGGk_(w,Gk)
  use workplace
  implicit none
  
  complex*16 w,Gk(nambu,nambu,nk)
  
  integer ik,iorb,jorb,m,ig,ikg,iorbg(2),jorbg(2)
  integer ikir,ikk
  real*8 k(2),kg(2),iXg(2),jXg(2)
  logical reducible

  Gk=0
  ikk=0
  ikir=0
  do ik=1,nk
    k=kv(1:2,ik)
	if(reducible(k,group))cycle
	ikir=ikir+1
	ikk=ikk+1
	do iorb=1,nambu
    do jorb=1,nambu
      do m=1,nambu
	    Gk(iorb,jorb,ikir)=Gk(iorb,jorb,ikir)+BdGAk(iorb,m,ik)*conjg(BdGAk(jorb,m,ik))/(w-BdGek(m,ik))
	  end do
	end do
	end do
	do ig=2,ng
	  ikk=ikk+1
	  call groupaction(ig,k,kg,group)
	  do iorb=1,norb
	  do jorb=1,norb
	    call orbitimage(ig,iorb,norb,model%orbits,group,iorbg,iXg)
        call orbitimage(ig,jorb,norb,model%orbits,group,jorbg,jXg)
		Gk(iorbg(1),jorbg(1),ikk)=Gk(iorb,jorb,ikir)*iXg(1)*jXg(1)
	  end do
	  end do
	end do
	ikir=ikir+ng-1
	print*,ikk,ikir; pause
  end do
  print*,ikk,ikir; pause
  return
end subroutine

subroutine searchkv(ik,k)
  use workplace
  implicit none

  integer ik
  real*8 k(2)

  integer ikk
  real*8 dis,distance

  ik=0
  dis=100
  do ikk=1,nk
    if(distance(kv(1:2,ikk)-k)<1.e-6)then
	  ik=ikk
	  return
	end if
  end do

  return
end subroutine
  
subroutine BdGeigentable(delta,U)
  use workplace
  implicit none

  real*8 delta
  complex*16 U(ndim)

  integer ik
  real*8 k(2)
  complex*16 BdGhk(nambu,nambu)

  print*,'Setting BdG eigentable ...'
  allocate(BdGek(nambu,nk))
  allocate(BdGAk(nambu,nambu,nk))
  do ik=1,nk
    k=kv(1:2,ik)
	call getBdGHk(k,delta,U,BdGHk)
    call ZHEIGEN(nambu,BdGhk,BdGek(:,ik))
	BdGAk(:,:,ik)=BdGhk
  end do

  print*,'BdG eigentable done!'
  
  return
end subroutine


subroutine BdGband(delta,U)
  use workplace
  implicit none

  real*8 delta
  complex*16 U(ndim)

  integer is,ikx,iky,idim,nkx1,nky1
  real*8 u1,u2,s,kx,ky,k(2),BdGeval(nambu)
  complex*16 BdGhk(nambu,nambu)
  integer stat1
  
  
  print*,'Calculating BdG band...'
  nkx1=40; nky1=nkx1
  open(11,file='bdgband.dat')
  do ikx=1,nkx1+1; 
    kx=-1+2.d0/nkx1*(ikx-1)
	do iky=1,nky1+1
	  ky=-1+2.d0/nky1*(iky-1)
	  k=(/kx,ky/)
	  call getBdGhk(k,delta,U,BdGhk); !pause
	  call ZHEIGEN(nambu,BdGhk,BdGeval)
	  write(11,'(1e16.8)')BdGeval
	end do
  end do
  close(11)
  open(10,file='bdggap.dat')
  open(11,file='fs.dat')
  do while(.true.)
    read(11,*,iostat=stat1)k
	call getBdGhk(k,delta,U,BdGhk)
    call ZHEIGEN(nambu,BdGhk,BdGeval)
	write(10,'(1e16.8)')BdGeval(norb+1)
	if(stat1/=0)exit
  end do
  close(10)
  close(11)

  print*,'BdG band OK'

  return
end subroutine

subroutine getBdGHk(k,delta,U,BdGHk)
  use workplace
  implicit none
 
  real*8 k(2),delta
  complex*16 BdGhk(nambu,nambu)
  complex*16 U(ndim)
  
  integer is,iorb
  real*8 Q,eval(norb)
  complex*16 hk(norb,norb),gk(norb,norb)

  call gethk(norb,k,Q,hk,model,MF,.false.)
  call getgk(k,gk,U)
  gk=gk*delta
  
  BdGhk=0
  
  if(model%sample=='***')then
    call ZHEIGEN(norb,hk,eval)
    gk=matmul(conjg(transpose(hk)),matmul(gk,hk))
    do iorb=1,norb
      BdGhk(iorb,iorb)=eval(iorb)
    end do
    BdGhk(norb+1:nambu,norb+1:nambu)=-BdGhk(1:norb,1:norb)
  else
    BdGhk(1:norb,1:norb)=hk
    BdGhk(norb+1:nambu,norb+1:nambu)=-conjg(hk)
  end if
  BdGhk(1:norb,norb+1:nambu)=conjg(transpose(gk))
  BdGhk(norb+1:nambu,1:norb)=gk

  return
end subroutine
 
subroutine Tmatrix__()
  implicit none

  integer iw,nw,ip
  real*8 pi,delta,La,eta,V1,V2,akw,akw_,rho0
  real*8,dimension(2) :: mu,zeta,rho1
  complex*16 one,w,a
  real*8,dimension(2,2) :: sig0,sig1,sig3
  real*8,dimension(4,4) :: ga00,ga10,ga03,ga13,ga30,ga33,ga31
  complex*16,dimension(2,2,2) :: geh
  complex*16,dimension(4,4) :: T,g0
  integer Vim,sc
  
  sig0=0
  sig0(1,1)=1
  sig0(2,2)=1

  sig1=0
  sig1(1,2)=1
  sig1(2,1)=1
   
  sig3=sig1
  sig3(2,2)=-1

  ga00=0
  ga00(1:2,1:2)=sig0
  ga00(3:4,3:4)=sig0

  ga03=0
  ga03(1:2,1:2)=sig0
  ga03(3:4,3:4)=-sig0

  ga10=0
  ga10(1:2,1:2)=sig1
  ga10(3:4,3:4)=sig1

  ga13=0
  ga13(1:2,1:2)=sig1
  ga13(3:4,3:4)=-sig1

  ga30=0
  ga30(1:2,1:2)=sig3
  ga30(3:4,3:4)=sig3

  ga33=0
  ga33(1:2,1:2)=sig3
  ga33(3:4,3:4)=-sig3

  ga31=0
  ga31(1:2,3:4)=sig3
  ga31(3:4,1:2)=sig3

  pi=3.1415926
  one=cmplx(0.0,1.0)

  rho0=0.5
  rho1(1)=-1.5
  rho1(2)=rho1(1)
  zeta(1)=1
  zeta(2)=-1
  mu(1)=0.2
  mu(2)=0.2
  La=1
  delta=0.03
  eta=0.001
  V1=-1.5
  V2=0.5*V1
  sc=0
  Vim=0
  open(10,file='akw.dat')
  open(11,file='akw_.dat')
  open(12,file='w.dat')
  do iw=-99,99
    w=iw*2*delta/100
    write(12,'(1e16.8)')real(w)
    w=w+one*eta
	
	!get g0(0,w)
	a=sc*delta**2-w**2
	a=sqrt(a)
	do ip=1,2
	  geh(:,:,ip) = - pi*(rho0+zeta(ip)*rho1(ip)*mu(ip))*(w*sig0+sc &
	                                   !*delta*sig1)/a & 
						               *zeta(ip)*delta*sig1)/a &
	                + zeta(ip)*rho1(ip)*(La+mu(ip)-pi*a)*sig3 &
				    - log(mu(ip)/La)*( -rho1(ip)*w*sig0 - sc*rho1(ip) &
				                       !*delta*sig1 & 
						               *zeta(ip)*rho1(ip)*delta*sig1 &
				    + (zeta(ip)*rho0+rho1(ip)*mu(ip))*sig3 )
    end do
    g0=0
	g0(1:2,1:2)=geh(:,:,1)
	g0(3:4,3:4)=geh(:,:,2)
	
	!get T-matrix
	T=0
	T(1:2,1:2)=V1*sig3
	T(3:4,3:4)=V1*sig3
	T(1:2,3:4)=-V2*sig3
    T(3:4,1:2)=-V2*sig3
    T=T/(V1**2-V2**2)

	!get G(0,w)
	T=T-g0
	call ZINVERT(4,T)
	if(Vim==0)T=0
	g0=g0+matmul(g0,matmul(T,g0))

	!get DOS
    Akw=0
    Akw_=0
    Akw=-1/pi*(imag(g0(1,1))+imag(g0(3,3)))
	Akw_=-1/pi*(imag(g0(2,2))+imag(g0(4,4)))
    write(10,'(1e16.8)')Akw       !A(w) for spin up
    write(11,'(1e16.8)')Akw_      !A(-w) for spin down
  end do
  close(10)
  close(11)
  close(12)

  return
end subroutine

subroutine Tmatrix_()
  implicit none

  integer iw,nw,ig
  real*8 pi,Dt,La,eta,V1,V2,akw,akw_,wmax
  real*8,dimension(2) :: mu,zt,rho0,rho1,rho2,al,bt,ga
  complex*16 one,w,a
  real*8,dimension(2,2) :: sig0,sig1,sig3
  complex*16,dimension(2) :: int0,int1,int2,int3
  complex*16,dimension(2,2,2) :: geh
  complex*16,dimension(4,4) :: T,g0
  integer Vim,sc

  sig0=0
  sig0(1,1)=1
  sig0(2,2)=1

  sig1=0
  sig1(1,2)=1
  sig1(2,1)=1
   
  sig3=sig1
  sig3(2,2)=-1

  pi=3.1415926
  one=cmplx(0.0,1.0)

  rho0=0.3
  rho1(1)=0.
  rho1(2)=-0.
  rho2=0.0
  zt(1)=1
  zt(2)=-1
  mu(1)=0.2
  mu(2)=0.2
  Dt=0.03
  eta=0.001
  La=3
  
  sc=1
  Vim=1
  V1=0.1
  V2=V1*0.5
  
  al=rho0+rho1*mu+rho2*mu*mu
  bt=-rho1-2*rho2*mu
  ga=rho2
  wmax=Dt*2
  Dt=Dt*sc
  
  open(10,file='akw.dat')
  open(11,file='akw_.dat')
  open(12,file='w.dat')
  do iw=-99,99
    w=iw*wmax/100
    write(12,'(1e16.8)')real(w)
    w=w+one*eta
	a=sc*Dt**2-w**2
	!print*,a; pause
	a=sqrt(a)
	!print*,a; pause
    
	!get the integrals
    int0=-pi/a
	int1=log(La/mu) !; print*,int1; pause
	int2=-La-mu+pi*a
	int3=-0.5*(La**2+mu**2)-a**2*int1

	!get g0(0,w)
	do ig=1,2
	  geh(:,:,ig) = ( al(ig)*w*int0(ig)+bt(ig)*w*int1(ig)+ga(ig)*w*int2(ig) )*sig0 &
	             + ( al(ig)*zt(ig)*Dt*int0(ig)+bt(ig)*zt(ig)*Dt*int1(ig)+ga(ig)*zt(ig)*Dt*int2(ig) )*sig1 &
				 + ( al(ig)*zt(ig)*int1(ig)+bt(ig)*zt(ig)*int2(ig)+ga(ig)*zt(ig)*int3(ig) )*sig3
	  !geh(:,:,ig) = al(ig)*w*int0(ig)*sig0
	end do
    g0=0
	g0(1:2,1:2)=geh(:,:,1)
	g0(3:4,3:4)=geh(:,:,2)
	
	!get T-matrix
	T=0
	T(1:2,1:2)=V1*sig3
	T(3:4,3:4)=V1*sig3
	T(1:2,3:4)=-V2*sig3
    T(3:4,1:2)=-V2*sig3
    T=T/(V1**2-V2**2)

	!get G(0,w)
	T=T-g0
	call ZINVERT(4,T)
	!print*,T; pause
	if(Vim==0)T=0
	g0=g0+matmul(g0,matmul(T,g0))

	!get DOS
    Akw=0
    Akw_=0
    Akw=-1/pi*(imag(g0(1,1))+imag(g0(3,3)))
	Akw_=-1/pi*(imag(g0(2,2))+imag(g0(4,4)))
    write(10,'(1e16.8)')Akw       !A(w) for spin up
    write(11,'(1e16.8)')Akw_      !A(-w) for spin down
  end do
  close(10)
  close(11)
  close(12)

  return
end subroutine


	!g0=-(rhoh*alphah+rhoe*alphae)*w/a*ga00-(rhoh*alphah-rhoe*alphae)*delta/a*ga10+(rhoh*betah+rhoe*betae)*ga30 &
	!   -(rhoh*alphah-rhoe*alphae)*w/a*ga03-(rhoh*alphah+rhoe*alphae)*delta/a*ga13+(rhoh*betah-rhoe*betae)*ga33
	!g0=g0*0.5
    !g0=-pi*rho0*w/a*ga00-pi*rho0*delta/a*ga13+rho0*alpha*ga30+rho0*beta*ga33	  
