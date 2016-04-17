subroutine meanfield()
  use workplace
  implicit none
  
	call filling2mu(mf%slaveboson)
  call initialize()
  call selfconsistent(model,MF) 
  call deltaonfs()

  return
end subroutine

subroutine initialize()
  use workplace
  implicit none
  
  integer ibond,iorb,jorb,ix,itemp,jtemp
	real*8 bond(2)
  real*8,dimension(MF%nbond,norb1atom,norb1atom) :: ranx,rand
  
	!allocate the order parameters
	do ibond=1,MF%nbond
    allocate(MF%chi(ibond)%fo(norb1atom,norb1atom))
		allocate(MF%delta(ibond)%fo(norb1atom,norb1atom))
		allocate(MF%mag(ibond)%fo(norb1atom,norb1atom))
		MF%chi(ibond)%fo=0
		MF%delta(ibond)%fo=0
		MF%mag(ibond)%fo=0
	end do

	!initialize 
	call random_seed()
	call random_number(ranx)
	call random_number(rand)
	ranx=0.3*ranx
	rand=0.3*(rand-0.5)
  do ibond=1,MF%nbond
	  if(MF%dochi)MF%chi(ibond)%fo=ranx(ibond,:,:)
		if(MF%dosc)MF%delta(ibond)%fo=rand(ibond,:,:)
		if(MF%domag.and.sum(abs(MF%mag(ibond)%bond))<1.e-5)MF%mag(ibond)%fo=rand(ibond,:,:)
	end do 
  print*,'mean field initialization OK!'

  !read data from breakpoint
  if(mf%breakpoint.or.mf%finished)then
		!initialize from file
  	open(10,file='chi.input')
		open(11,file='delta.input')
		open(12,file='mag.input')
		open(13,file='mu.input')
		do ibond=1,MF%nbond
	  	do iorb=1,norb1atom
			do jorb=1,norb1atom
		  	read(10,100)MF%chi(ibond)%fo(iorb,jorb)
		  	read(11,100)MF%delta(ibond)%fo(iorb,jorb)
				read(12,*)MF%mag(ibond)%fo(iorb,jorb)
	  	end do
			end do
  	end do
    read(13,*)model%mu
100 format(2e16.8)
 	 	close(10)
		close(11)
		close(12)
		close(13)
  endif

  call symmetrizationmf()

	return
end subroutine

subroutine selfconsistent(model,mf)
  use standard_derived_types
  implicit none
  
  type (modelblock) :: model
  type (mfblock) :: mf

  integer norb1atom,ibond,iorb,jorb,icount,nl
  real*8 errorx,errord,errorm,error,eta,errmin,dp
  type (orderparameter),allocatable,dimension(:) :: workx,workd,workm

  norb1atom=model%norb/model%natom
  allocate(workx(mf%nbond),workd(mf%nbond),workm(mf%nbond))
	nl=mf%nl
	eta=mf%eta
	errmin=mf%errmin
	icount=0
  do while(.not.mf%finished)
	  icount=icount+1
    errorx=0
    errord=0
		errorm=0
    workx=mf%chi
    workd=mf%delta
		workm=mf%mag
		call neworderparameter(model%norb,mf%nambu,nl,model,mf,dp)
    call symmetrizationmf()
		do ibond=1,mf%nbond
      do iorb=1,norb1atom
      do jorb=1,norb1atom
			  if(mf%dochi)then
        	errorx=errorx+abs(workx(ibond)%fo(iorb,jorb)-mf%chi(ibond)%fo(iorb,jorb))
        	mf%chi(ibond)%fo(iorb,jorb)=(1-eta)*workx(ibond)%fo(iorb,jorb) &
          	+eta*mf%chi(ibond)%fo(iorb,jorb)
				end if
				if(mf%dosc)then
				  errord=errord+abs(workd(ibond)%fo(iorb,jorb)-mf%delta(ibond)%fo(iorb,jorb))
        	mf%delta(ibond)%fo(iorb,jorb)=(1-eta)*workd(ibond)%fo(iorb,jorb) &
          	+eta*mf%delta(ibond)%fo(iorb,jorb)
				end if
				if(mf%domag)then
					errorm=errorm+abs(workm(ibond)%fo(iorb,jorb)-mf%mag(ibond)%fo(iorb,jorb))
					mf%mag(ibond)%fo(iorb,jorb)=(1-eta)*workm(ibond)%fo(iorb,jorb) &
				  	+eta*mf%mag(ibond)%fo(iorb,jorb)
				end if
      end do
      end do
    end do
		error=max(errorx,errord)
		error=max(error,errorm)
		error=max(error,abs(dp))
		print*,icount,'error=',error
		if(error<errmin)exit
		if(mod(icount,mf%nstep)==0)then
			call outputmf(.false.)
			!call bdgmf(-2.d0,2.d0,200,.true.)
			!call deltaonfs()
		  print*,'results output to file!'
		end if

  end do
  call outputmf(.true.)
  call bdgmf(-4.d0,4.d0,400,.false.)
	deallocate(workx,workd)
  
  return
end subroutine

subroutine tunemu()
  use workplace
	implicit none

  integer icount
  real*8 emax,emin,dope,mu2dopemf

  emax=5
  emin=-5
	icount=0
  do while(.true.)
	  icount=icount+1
    model%mu=(emax+emin)/2
		dope=mu2dopemf()
		if(dope>abs(MF%nhole))emin=model%mu
		if(dope<abs(MF%nhole))emax=model%mu
		if(abs(dope-abs(MF%nhole))<mf%errmin4mu)exit
		!print*,icount,model%mu,dope
  end do
  print*,'mu tuned by doping level:',model%mu 
	return
end subroutine 
  
function mu2dopemf()
  use workplace
	implicit none

  real*8 mu2dopemf
	integer nx,ny,ix,iy,idim,iidim,jdim,jjdim,kdim
	real*8 kx,ky,k(2),ekmf(mf%nambu),dope,fk,fermi
	complex*16 Hmf(mf%nambu,mf%nambu)
  logical reducible

  dope=0
  nx=mf%nl
	ny=nx
	do ix=1,nx
  do iy=1,ny
    kx=2.d0/nx*(ix-1)
    ky=2.d0/ny*(iy-1)
		k=(/kx,ky/)
		call gethkmf(norb,mf%nambu,k,Hmf,model,mf,mf%slaveboson)
		call ZHEIGEN(mf%nambu,Hmf,ekmf)
		do idim=1,norb 
		  iidim=idim+norb
		  do kdim=1,mf%nambu
		    fk=fermi(ekmf(kdim),mf%temp)
		    dope=dope+fk*conjg(Hmf(idim,kdim))*Hmf(idim,kdim)+(1-fk)*Hmf(iidim,kdim)*conjg(Hmf(iidim,kdim))
	    end do
    end do
  end do
	end do
		
  dope=dope/nx/ny/model%natom
	dope=1-dope
  mu2dopemf=dope

	return
end function

subroutine neworderparameter(norb,nambu,nl,model,mf,dp)
  use standard_derived_types
  implicit none

  integer nambu,nl
  type (modelblock) :: model
  type (mfblock) :: mf

  integer norb,nx,ny,ix,iy,idim,jdim,kdim,ibond,iatom,jatom,iorb,jorb,norb1atom,iidim,jjdim
  real*8 pi,kx,ky,kv(2),ekmf(nambu),bond(2),fermi,dope,fk(nambu),dp
  complex*16 workx,workd,workm,phase
	complex*16 one,Hmf(nambu,nambu)
	complex*16,dimension(norb,norb,nl,nl) ::  fxkup,fxkdn,fdk,fmkup,fmkdn
 
  if(mf%tunemu)then
	  print*,'tuning mu ...'
		call tunemu()
	end if
  norb1atom=model%norb/model%natom
  pi=2*asin(1.d0)
  one=cmplx(0.d0,1.d0)
  nx=nl
  ny=nl
  dope=0
	fxkup=0
	fxkdn=0
	fdk=0
	fmkup=0
	fmkdn=0
  do ix=1,nx
  do iy=1,ny
    kx=2.d0/nx*(ix-1)
    ky=2.d0/ny*(iy-1)
    kv=(/kx,ky/)
		call gethkmf(norb,nambu,kv,Hmf,model,mf,mf%slaveboson)
		call ZHEIGEN(nambu,Hmf,ekmf)
		do kdim=1,nambu
		  fk(kdim)=fermi(ekmf(kdim),mf%temp)
		end do
    do iorb=1,norb
		  do jorb=1,norb
			  fxkup(iorb,jorb,ix,iy)=sum(fk*conjg(Hmf(iorb,:))*Hmf(jorb,:))
				fxkdn(iorb,jorb,ix,iy)=sum((1-fk)*Hmf(iorb+norb,:)*conjg(Hmf(jorb+norb,:)))
				fdk(iorb,jorb,ix,iy)=sum((1-fk)*Hmf(iorb,:)*conjg(Hmf(jorb+norb,:)))
				if(mf%domag)then
				  if(mf%nambu==norb*4)then
				    fmkup(iorb,jorb,ix,iy)=sum(fk*conjg(Hmf(iorb,:))*Hmf(jorb+2*norb,:))
						fmkdn(iorb,jorb,ix,iy)=sum((1-fk)*Hmf(iorb+norb,:)*conjg(Hmf(jorb+3*norb,:)))
					else
				  	fmkup(iorb,jorb,ix,iy)=fxkup(iorb,jorb,ix,iy)
						fmkdn(iorb,jorb,ix,iy)=fxkdn(iorb,jorb,ix,iy)
					end if
				end if
	    end do
    end do
		do iorb=1,norb
		  dope=dope+fxkup(iorb,iorb,ix,iy)+fxkdn(iorb,iorb,ix,iy)	
		end do
  end do
  end do
	dope=dope/nx/ny/model%natom
	dope=1-dope
	dp=dope-abs(mf%nhole)
  print*,'mu',model%mu,'dp:',dp
	model%mu=model%mu+dp*abs(mf%nhole)
  if(abs(dp)>mf%errmin4mu)then
	  mf%tunemu=.true.
  else 
	  mf%tunemu=.false.
	end if
	do ibond=1,mf%nbond
    iatom=mf%chi(ibond)%atom
    jatom=mf%chi(ibond)%nbatom
		bond=mf%chi(ibond)%bond
    do iorb=1,norb1atom
    do jorb=1,norb1atom
      idim=iorb+(iatom-1)*norb1atom
      jdim=jorb+(jatom-1)*norb1atom
      workx=0
      workd=0
			workm=0
      do ix=1,nx
      do iy=1,ny
        kx=2.d0*pi/nx*(ix-1)
        ky=2.d0*pi/ny*(iy-1)
        kv=(/kx,ky/)
				phase=exp(-one*sum(kv*bond)) !-model%ra(1:2,iatom)+model%ra(1:2,jatom)))
        if(mf%dochi)workx=workx+fxkup(idim,jdim,ix,iy)*phase &
				                       +fxkdn(idim,jdim,ix,iy)*conjg(phase)
				if(mf%dosc)workd=workd+fdk(idim,jdim,ix,iy)*conjg(phase) &
				                      +fdk(jdim,idim,ix,iy)*phase
				if(mf%domag)workm=workm+fmkup(idim,jdim,ix,iy)*phase &
				                       -fmkdn(idim,jdim,ix,iy)*conjg(phase)
      end do
      end do
      mf%chi(ibond)%fo(iorb,jorb)=workx/nx/ny
      mf%delta(ibond)%fo(iorb,jorb)=workd/nx/ny
			mf%mag(ibond)%fo(iorb,jorb)=workm/nx/ny/2
			if(abs(mf%J(ibond,iorb,jorb))<1.e-5)then
			  mf%delta(ibond)%fo(iorb,jorb)=0
				mf%chi(ibond)%fo(iorb,jorb)=0
			end if
			if(sum(abs(bond))>1.e-5)then
			  mf%mag(ibond)%fo(iorb,jorb)=0
			end if
    end do
    end do
  end do

  return
end subroutine

subroutine symmetrizationmf()
  use workplace
  implicit none

  integer iatom,ibond,idim,jbond,jdim,findbond,jatom
  real*8 bondi(2),bondj(2), m

  !caution: only apply to lhq mod
  !Inversion symmetry: general
  do iatom=1,model%natom
    do ibond=1,mf%nbond1atom
      idim=ibond+(iatom-1)*mf%nbond1atom
      bondi=mf%bonds(:,ibond)
      bondj=-bondi
      jbond=findbond(bondj,mf)
      jdim=jbond+(iatom-1)*mf%nbond1atom
      mf%chi(jdim)%fo=conjg(transpose(mf%chi(idim)%fo))
      mf%delta(jdim)%fo=conjg(transpose(mf%delta(idim)%fo))
    end do
  end do
  !special for lhq model
  do ibond=10,18
    if(ibond==10)then
      m=real(mf%mag(10)%fo(1,1))
      mf%mag(10)%fo=-mf%mag(1)%fo
    end if
    if(ibond>=11.and.ibond<=14)then
      mf%chi(ibond)%fo=mf%chi(ibond-9)%fo
      mf%delta(ibond)%fo=mf%delta(ibond-9)%fo
    end if
    if(ibond==15.or.ibond==16)then
      mf%chi(ibond)%fo=mf%chi(ibond-7)%fo
      mf%delta(ibond)%fo=mf%delta(ibond-7)%fo
    end if
    if(ibond==17.or.ibond==18)then
      mf%chi(ibond)%fo=mf%chi(ibond-11)%fo
      mf%delta(ibond)%fo=mf%delta(ibond-11)%fo
    end if
  end do
 
  if(model%sample=='lhq2')then
	  do iatom=3,8
			jatom=mod(iatom+1,2)+1
			do ibond=1,9
				idim=ibond+(iatom-1)*mf%nbond1atom
				jdim=ibond+(jatom-1)*mf%nbond1atom
				mf%chi(idim)%fo=mf%chi(jdim)%fo
				mf%delta(idim)%fo=mf%delta(jdim)%fo
			end do
		end do
	end if
  return
end subroutine

function findbond(bond,mf)
  use standard_derived_types
  implicit none

  integer findbond
  real*8 bond(2)
  type(mfblock) :: mf

  integer ibond
  real*8 bondi(2)

  do ibond=1,mf%nbond1atom
    bondi=mf%bonds(:,ibond)
    if(sum(abs(bondi-bond))<1.e-5)then
      findbond=ibond
      exit
    end if
  end do
  return
end function
   
subroutine gethkmf(norb,nambu,kv,Hmf,model,mf,slvbsn)
  use standard_derived_types
  implicit none
  
  integer norb,nambu
  real*8 kv(2)
  complex*16 Hmf(nambu,nambu)
  type (modelblock) :: model
  type (mfblock) :: mf
	logical slvbsn

  integer idim,jdim,norb1atom,nbond1atom,iorb,jorb,iatom,jatom,ibond
	complex*16 h(norb)
  complex*16 hk(norb,norb),xk(norb,norb),dk(norb,norb)
  
  Hmf=0

	call gethk(norb,kv,0.d0,hk,model,mf,slvbsn)
	call getxkdk(norb,kv,xk,dk,model,mf)

	if(mf%dochi)hk=hk+xk
  Hmf(1:norb,1:norb)=hk
	Hmf(norb+1:norb*2,norb+1:norb*2)=-hk
  
	if(mf%dosc)then
	  Hmf(1:norb,norb+1:norb*2)=conjg(transpose(dk))
	  Hmf(norb+1:norb*2,1:norb)=dk
	end if
	
	if(mf%domag)then
		call gethmag(norb,h,model,mf)
    if(mf%nambu==norb*2)then
		  do idim=1,norb
			  hmf(idim,idim)=hmf(idim,idim)+h(idim)*0.5
				hmf(idim+norb,idim+norb)=hmf(idim+norb,idim+norb)+h(idim)*0.5
			end do
		end if
		!for Qsdw=(pi,pi)
		if(mf%nambu==norb*4)then
			kv=kv+mf%qsdw
			call gethk(norb,kv,0.d0,hk,model,mf,slvbsn)
	    call getxkdk(norb,kv,xk,dk,model,mf)
    
		  if(mf%dochi)hk=hk+xk
		  Hmf(2*norb+1:3*norb,2*norb+1:3*norb)=hk
	    Hmf(3*norb+1:4*norb,3*norb+1:4*norb)=-hk

      if(mf%dosc)then
	      Hmf(2*norb+1:3*norb,3*norb+1:norb*4)=conjg(transpose(dk))
	      Hmf(3*norb+1:norb*4,2*norb+1:3*norb)=dk
	    end if

	  	do idim=1,norb
		    Hmf(idim,idim+2*norb)=h(idim)*0.5
				Hmf(idim+norb,idim+3*norb)=h(idim)*0.5
		    Hmf(idim+2*norb,idim)=conjg(hmf(idim,idim+2*norb))
				Hmf(idim+3*norb,idim+norb)=conjg(hmf(idim+norb,idim+3*norb))
		  end do
	  end if
	end if
  
  do idim=1,nambu
  do jdim=1,nambu
  	if(abs(Hmf(idim,jdim)-conjg(Hmf(jdim,idim)))>1.e-5)then
	  print*,idim,jdim
	  print*,Hmf(idim,jdim),Hmf(jdim,idim)
	  read*
	end if
  end do
  end do

  return
end subroutine

subroutine gethmag(norb,h,model,mf)
  use standard_derived_types
	implicit none

	integer norb
	complex*16 h(norb)
	type (modelblock) :: model
	type (mfblock) :: mf

	integer norb1atom,nbond1atom
  integer iatom,jatom,iorb,jorb,idim,jdim,ibond
	real*8 bond(2),qm(2),pi
	complex*16 m(norb),one

  pi=2*asin(1.d0)
	one=cmplx(0.d0,1.d0)

  norb1atom=model%norb/model%natom
	nbond1atom=mf%nbond1atom
	m=0
	do iatom=1,model%natom
	  ibond=1+(iatom-1)*nbond1atom
		do iorb=1,norb1atom
		  idim=iorb+(iatom-1)*norb1atom
			m(idim)=mf%mag(ibond)%fo(iorb,iorb)
		end do
	end do
	h=0
	do iatom=1,model%natom
	 	do iorb=1,norb1atom
	  	idim=iorb+(iatom-1)*norb1atom
			do ibond=1,mf%nbond
				if(mf%mag(ibond)%atom/=iatom)cycle
				bond=mf%mag(ibond)%bond
				qm=mf%qsdw
				jatom=mf%mag(ibond)%nbatom
				bond=bond !-model%ra(1:2,iatom)+model%ra(1:2,jatom)
			  do jorb=1,norb1atom
				  jdim=jorb+(jatom-1)*norb1atom
					h(idim)=h(idim)+mf%J(ibond,iorb,jorb)*m(jdim)*exp(one*pi*sum(qm*bond))
				end do
			end do
			!h(idim)=4*mf%J(2,1,1)*mf%mag(1)%fo(1,1)
		end do
	end do
	
	return
end subroutine

subroutine getxkdk(norb,kv,xk,dk,model,mf)
  use standard_derived_types
  implicit none

  integer norb
  real*8 kv(2)
  complex*16 xk(norb,norb),dk(norb,norb)
  type (modelblock) :: model
  type (mfblock) :: mf

  integer norb1atom,iatom,jatom,ibond,iorb,jorb
  real*8 pi,vk(2),bond(2)
  complex*16 one
  complex*16,allocatable,dimension(:,:) :: workx,workd
  !complex*16 workx,workd
 
  pi=2*asin(1.d0)
  one=cmplx(0.d0,1.d0)
  vk=kv*pi

  xk=0
  dk=0
  norb1atom=norb/model%natom
  allocate(workx(norb1atom,norb1atom),workd(norb1atom,norb1atom))

  do iatom=1,model%natom
  do jatom=1,model%natom
    workx=0
    workd=0
	  do ibond=1,mf%nbond
		  bond=mf%chi(ibond)%bond   !-model%ra(1:2,iatom)+model%ra(1:2,jatom)
	    if(mf%chi(ibond)%atom/=iatom.or.mf%chi(ibond)%nbatom/=jatom)cycle
      workx=workx-3.d0/8*mf%J(ibond,:,:)*conjg(mf%chi(ibond)%fo) &
        * exp(-one*sum(vk*bond))
      workd=workd+3.d0/8*mf%J(ibond,:,:)*conjg(mf%delta(ibond)%fo) & 
        * exp(-one*sum(vk*bond))
	  end do
		xk(norb1atom*(iatom-1)+1:norb1atom*iatom,norb1atom*(jatom-1)+1:norb1atom*jatom)=workx
    dk(norb1atom*(iatom-1)+1:norb1atom*iatom,norb1atom*(jatom-1)+1:norb1atom*jatom)=workd
  end do
  end do
  
	do iorb=1,model%norb
	  do jorb=1,iorb-1
	    xk(iorb,jorb)=conjg(xk(jorb,iorb))
	  end do
    xk(iorb,iorb)=real(xk(iorb,iorb))
	end do
	  
  deallocate(workx,workd)

  return
end subroutine

subroutine bdgmf(wmin,wmax,nw,slvbsn)
  use workplace, only : mf,model,norb,nambu,pi,one
	implicit none
  
	integer nw
  real*8 wmax,wmin
  logical slvbsn

	integer ik,ikx,iky,nx,ny,iorb,iw,nw_,mw,icount
	integer,allocatable :: windex(:)
	real*8 k(2),G(2),X(2),M(2),eval(mf%nambu),kx,ky
	real*8 w,dos(nw),dosup(nw),dosdown(nw),dw,eta
	real*8 wup,wdown,filling,dw_,mu0
	real*8,allocatable,dimension(:) :: w_,dos_
	complex*16 Hmf(mf%nambu,mf%nambu),gkw(mf%nambu),gkw_(mf%nambu)

  eta=1.e-2
	G=(/0,0/)
	X=(/1,0/)*1.d0
	M=(/1,1/)*1.d0
  open(10,file='bdgband.dat')
	do ik=0,99
	  k=G+ik*(M-G)/100
		call gethkmf(norb,mf%nambu,k,Hmf,model,mf,slvbsn)
		call ZHEIGEN(mf%nambu,Hmf,eval)
		write(10,'(1e16.8)')eval
	end do
	do ik=0,99
	  k=M+ik*(X-M)/100
		call gethkmf(norb,mf%nambu,k,Hmf,model,mf,slvbsn)
		call ZHEIGEN(mf%nambu,Hmf,eval)
		write(10,'(1e16.8)')eval
	end do
	do ik=0,100
	  k=X+ik*(G-X)/100
    call gethkmf(norb,mf%nambu,k,Hmf,model,mf,slvbsn)
		call ZHEIGEN(mf%nambu,Hmf,eval)
		write(10,'(1e16.8)')eval
	end do
  close(10)
  print*,'BdG band output to file'

	open(10,file='bdgdos.dat')
	nx=400
	ny=nx
	dos=0
	dosup=0
	dosdown=0
	dw=(wmax-wmin)/nw
	do ikx=1,nx
	  kx=2.d0*ikx/nx
	  do iky=1,ny
		  ky=2.d0*iky/ny
			k=(/kx,ky/)
			call gethkmf(norb,mf%nambu,k,Hmf,model,mf,slvbsn)
			call ZHEIGEN(mf%nambu,Hmf,eval)
      do iw=1,nw
			  w=wmin+dw*(iw-1)
			  gkw=1.d0/(w+one*eta-eval)
			  gkw_=1.d0/(-w+one*eta-eval)
				do iorb=1,norb
				  dosup(iw)=dosup(iw)-1/pi*sum(Hmf(iorb,:)*conjg(Hmf(iorb,:))*imag(gkw))
					dosdown(iw)=dosdown(iw)-1/pi*sum(Hmf(iorb+norb,:)*conjg(Hmf(iorb+norb,:))*imag(gkw_))
				end do
			end do
		end do
	end do
	dosup=dosup/nx/ny/norb
	dosdown=dosdown/nx/ny/norb
	dos=dosup+dosdown
	do iw=1,nw
	  w=wmin+(iw-1)*dw
	  write(10,'(4e16.8)')w,dosup(iw),dosdown(iw),dos(iw)
	end do
  close(10)
  print*,'dos output to file'		  
  
	return
end subroutine


function mu2filling(slvbsn)
  use workplace
  implicit none

  real*8 mu2filling
  logical slvbsn

  integer nx,ny,ix,iy,iorb
  real*8 kx,ky,k(2),fermi,eval(norb),temp
  complex*16 hk(norb,norb)
  

  temp=1.e-5
  nx=100
  ny=100
  mu2filling=0
  do ix=1,nx
  do iy=1,ny
    kx=2.d0*(ix-1)/nx
		ky=2.d0*(iy-1)/ny
		k=(/kx,ky/)
		call gethk(norb,k,0.d0,hk,model,MF,slvbsn)
    !call gethk(norb,k,hk,model,MF)
		call ZHEIGEN(norb,hk,eval)
    do iorb=1,norb
	  	mu2filling=mu2filling+fermi(eval(iorb),temp)
		end do
  end do
  end do
  mu2filling=2*mu2filling/nx/ny/model%natom
  return
end function

subroutine filling2mu(slvbsn)
  use workplace
  implicit none

  logical slvbsn

  integer icount
  real*8 emax,emin,filling,mu2filling

  emax=model%mu+5
  emin=model%mu-5
	icount=0
  print*,'tuning mu ...'
  do while(.true.)
	  icount=icount+1
    model%mu=(emax+emin)/2
		filling=mu2filling(slvbsn)
    print*,icount,model%mu,filling
		if(filling<model%filling)emin=model%mu
		if(filling>model%filling)emax=model%mu
		if(abs(filling-model%filling)<1.e-5)exit
  end do
  print*,'model%mu determined by nhole:',model%mu

  return
end subroutine

subroutine outputmf(finished)
  use workplace
	implicit none

  logical finished

  integer ibond,iorb,jorb,iatom,jatom
	real*8 mfenergy,bond(2)
	character*3 fileid

	open(10,file='mfresult.dat')
	open(11,file='chi.input')
	open(12,file='delta.input')
	open(13,file='mag.input')
	open(14,file='mu.input')
    write(10,*)'X-',mf%dochi,'D-',mf%dosc
    write(10,*)'M-',mf%domag,'Qsdw',mf%qsdw
	  write(10,*)'temperature:',MF%temp
 	  write(10,*)'lattice size:',MF%nl,'*',MF%nl
	  write(10,*)'errmin,errmin4mu,eta:',MF%errmin,mf%errmin4mu,MF%eta
	  write(10,*)'nhole: ',MF%nhole
		write(10,*)'nbond: ',MF%nbond
		write(10,*)'J: '
		do ibond=1,MF%nbond
		  iatom=MF%chi(ibond)%atom
			jatom=MF%chi(ibond)%nbatom
			bond=MF%chi(ibond)%bond
		  do iorb=1,norb1atom
			do jorb=1,norb1atom
			  write(10,100)bond,iatom,jatom,iorb,jorb,MF%J(ibond,iorb,jorb)
      end do
			end do
		end do
	  write(10,*)' '
		write(10,*)'Chi:'
		do ibond=1,MF%nbond
	  	iatom=MF%chi(ibond)%atom
			jatom=MF%chi(ibond)%nbatom
			bond=MF%chi(ibond)%bond
      do iorb=1,norb1atom
			do jorb=1,norb1atom
		  	write(10,100)bond,iatom,jatom,iorb,jorb,MF%chi(ibond)%fo(iorb,jorb)
				write(11,'(2e16.8)')MF%chi(ibond)%fo(iorb,jorb)
		  end do
		  end do
    end do
		write(10,*)' '
		write(10,*)'delta:'
		do ibond=1,MF%nbond
	  	iatom=MF%chi(ibond)%atom
			jatom=MF%chi(ibond)%nbatom
			bond=MF%chi(ibond)%bond
      do iorb=1,norb1atom
			do jorb=1,norb1atom
		    write(10,100)bond,iatom,jatom,iorb,jorb,MF%delta(ibond)%fo(iorb,jorb)
				write(12,'(2e16.8)')MF%delta(ibond)%fo(iorb,jorb)
	    end do
		  end do
    end do
		write(10,*)' '
		write(10,*)'m:'
		do ibond=1,MF%nbond
	  	iatom=MF%chi(ibond)%atom
			jatom=MF%chi(ibond)%nbatom
			bond=MF%chi(ibond)%bond
      do iorb=1,norb1atom
			do jorb=1,norb1atom
		    write(10,100)bond,iatom,jatom,iorb,jorb,MF%mag(ibond)%fo(iorb,jorb)
				write(13,'(2e16.8)')MF%mag(ibond)%fo(iorb,jorb)
	    end do
		  end do
    end do
		write(14,*)model%mu

100 format(2f10.5, 4i4, 2e16.8)
    if(finished)write(10,*)mfenergy(model%norb,mf%nambu,model,mf)
	close(10)
  close(11)
	close(12)
	close(13)
	close(14)

	return
end subroutine

subroutine deltaonfs()
  use workplace
	implicit none

	integer ik
	real*8 k(2),eval(norb),Q
  complex*16 hk(norb,norb),xk(norb,norb),dk(norb,norb)


  open(10,file='gapform.dat')
	do ik=1,model%nkf
	  k=model%kf(:,ik)
		call gethk(norb,k,Q,hk,model,MF,.false.)
		call ZHEIGEN(norb,hk,eval)
    call getxkdk(norb,k,xk,dk,model,MF)
    dk=matmul( conjg(transpose(hk)), matmul( dk, hk ) )
		write(10,'(4e16.8)')k,dk(model%band(ik),model%band(ik))
	end do
	close(10)

	return
end subroutine

function mfenergy(norb,nambu,model,mf)
  use standard_derived_types
	implicit none

  real*8 mfenergy
  integer norb,nambu
	type (modelblock) :: model
	type (mfblock) :: mf

  integer ix,iy,nx,ny,iorb,jorb,idim,jdim,kdim,norb1atom
	real*8 kx,ky,k(2),Q,ekmf(nambu),fk(nambu),fkh(nambu),fermi
	complex*16 emf,fok,hk(norb,norb),Hmf(nambu,nambu)

  call filling2mu(.false.); 	
	norb1atom=model%norb/model%natom
	emf=0
	nx=mf%nl
	ny=nx
	do ix=1,nx
  do iy=1,ny
    kx=2.d0/nx*(ix-1)
    ky=2.d0/ny*(iy-1)
    k=(/kx,ky/)
		call gethk(norb,k,Q,hk,model,mf,.false.)
    call gethkmf(norb,nambu,k,Hmf,model,mf,.false.)
		call ZHEIGEN(nambu,Hmf,ekmf)
		do kdim=1,nambu
		  fk(kdim)=fermi(ekmf(kdim),mf%temp)
			fkh(kdim)=1-fk(kdim)
		end do
		emf=emf+sum(fk*ekmf)
		do iorb=1,norb
		  emf=emf+hk(iorb,iorb)
		end do
	end do
	end do

  emf=emf/nx/ny
	
	mfenergy=emf

	return
end function

