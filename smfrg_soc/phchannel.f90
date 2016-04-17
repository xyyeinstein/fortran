subroutine phchannel(la,q,nA,A)
  use workplace, only : Mformfunc,Mformconfig
  implicit none
  real*8 q(2),la
  integer nA,iS,iSS,info(nA),mode(nA)
	real*8 Smax,Si(nA),Umax
  complex*16, dimension (nA,nA) :: A
	logical isave

  complex*16, dimension (nA,nA) :: B,U,V
	real*8 Sig(nA)
  integer idim,jdim,iatom,iorb,iform

  character*3 Mform,Lform,Oform,form(nA)
  type (Mformconfig), pointer :: Mi
  real*8, pointer, dimension (:,:) :: ra
	real*8 pi,explain
  complex*16 phase,one


  B=A; call ZSVD( nA,nA,B,U,Sig,V )
  do iS=1,nA
    Si(iS)=Sig(iS)*real( sum( U(:,iS)*V(iS,:) ) )
  end do
  Si=-Si; call ascendingorder(nA,Si,info); Si=-Si
  Smax=Si(1)

  open(45,position='append',file='Sph.dat')
  open(46,position='append',file='Sph.txt')
  do iS=1,2; iSS=info(iS)
    Umax=abs(U(1,iSS)); mode(iS)=1; form(iS)=mformfunc(1)%form
    do idim=2,nA
      if(abs(U(idim,iSS))>Umax)then
        Umax=abs(U(idim,iSS)); mode(iS)=idim; form(iS)=mformfunc(idim)%form
      end if
    end do
  end do
  write(45,110)Si(1),mode(1),form(1),Si(2),mode(2),form(2)
  write(46,120)la,Si(1),mode(1),Si(2),mode(2)
110 format(1x,2(f20.10,i5,a5))
120 format(1x,1f10.5,2(f20.10,i5))
  close(45)
  close(46)

  return
end subroutine phchannel
