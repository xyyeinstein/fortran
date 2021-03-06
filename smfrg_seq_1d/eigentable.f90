subroutine eigentable()
  use workplace
  implicit none

  integer ix,iy,iorb,iband,i,n
  real*8 Q,vk(2)
  complex*16 hk(norb,norb)
  real*8 eval(norb),occup(norb),dos,f,eta,wL,wU,w,dw
  real*8 fermi
  
  if(sum(abs(model%kb))==0)then
    print*,'dos for ehm1d'
    call eigentable_1d()
	return
  endif
 
  Q=0
  allocate(ek(norb,-nkx:nkx,-nky:nky))
  occup=0
  do iy=-nky,nky; do ix=-nkx,nkx
    vk=0.5*ix*model%ka/nkx+0.5*iy*model%kb/nky
    call gethk(norb,vk,Q,hk,model,MF,.false.)
	call ZHEIGEN(norb,hk,eval)
	ek(:,ix,iy)=eval
	if(ix==nkx.or.(iy==nky.and.iy/=0))cycle
	do iband=1,norb; f=fermi(eval(iband),1.e-3)
	  occup=occup+f*abs(hk(:,iband))**2
	end do
  end do; end do
  occup=2*occup/((2*nkx)*(2*nky))
  
  print*,'filling=',sum(occup)
  do iorb=1,norb
     write(*,100)iorb,occup(iorb)
100  format(1x,'occupation on orbital-',i3,'=',f10.4)
  end do
  print*,'occupation ok'

  wL=-3; wU=3; n=601; dw=(wU-wL)/(n-1);  eta=0.01
  open(10,file='dos.dat')
  do i=1,n; w=wL+(i-1)*dw; dos=0
    do iy=-nky,nky; do ix=-nkx,nkx
	  if(ix==nkx.or.(iy==nky.and.iy/=0))cycle
      dos = dos - sum( aimag( 1/(one*eta+w-ek(:,ix,iy)) ) )
	end do; end do
    dos=dos/(2*nkx*2*nky*norb)
    write(10,*)w,dos
  end do
  close(10)
  print*,'dos output to file'
  
  stop
  return
end subroutine eigentable

subroutine eigentable_1d()
  use workplace
  implicit none

  integer ix,iy,iorb,iband,i,n
  real*8 Q,vk(2)
  complex*16 hk(norb,norb)
  real*8 eval(norb),occup(norb),dos,f,eta,wL,wU,w,dw
  real*8 fermi

  Q=0
  allocate(ek(norb,-nkx:nkx,0:0))
  occup=0
  do ix=-nkx,nkx
    vk=0.5*ix*model%ka/nkx
    call gethk(norb,vk,Q,hk,model,MF,.false.)
	call ZHEIGEN(norb,hk,eval)
	ek(:,ix,0)=eval
	if(ix==nkx)cycle
	do iband=1,norb
      f=fermi(eval(iband),1.d-3)
	  occup=occup+f*abs(hk(:,iband))**2
	end do
  end do
  occup=2*occup/(2*nkx)
  
  print*,'filling=',occup,sum(occup)
  do iorb=1,norb
     write(*,100)iorb,occup(iorb)
100  format(1x,'occupation on orbital-',i3,'=',f10.4)
  end do
  print*,'occupation ok'

  wL=-3; wU=3; n=601; dw=(wU-wL)/(n-1);  eta=0.01
  open(10,file='dos.dat')
  do i=1,n; w=wL+(i-1)*dw; dos=0
    do ix=-nkx,nkx-1
	  dos = dos - sum( aimag( 1/(one*eta+w-ek(:,ix,0)) ) )
	end do
    dos=dos/(2*nkx*norb)
    write(10,*)w,dos
  end do
  close(10)
  print*,'dos output to file'
  
  return
end subroutine 

