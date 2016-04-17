subroutine setMolecules()
  use workplace, only : model

  !if(sum(abs(model%kb))==0)then
  !  call setMolecules_1d()
  !  return
  !end if

  if( sum( model%a*model%b ) == 0 )then
    call setMolecules_square()
  else
    call setMolecules_LowSymmetry()
  end if
  
  return
end subroutine setMolecules


subroutine setMolecules_square()
  use workplace
  implicit none

  integer idim,iLform,iOform,i,i1,i2,iatom,nbatom
  integer findatom
  logical compatible
  type (Mformconfig), pointer :: Mi
  type (Oformconfig), pointer :: Oi
  type (Lformconfig), pointer :: Li

  character*3 combine
  real*8, pointer, dimension (:,:) :: ra

  if(group=='C2v'.or.group=='D2h'.or.group=='C2_')then
    print*,'Set Mform under Low Symmetry...'
	call setMolecules_LowSymmetry()
	return
  end if

  ra=>model%ra
  ndim=0
  
  do iatom=1,natom; do iOform=1,nOform; do iLform=1,nLform
     print*,ioform,ilform
     Oi=>Oformfunc(iOform); Li=>Lformfunc(iLform)
     if(SkipInterOrbitPair.and.Oi%pair(1,1)/=Oi%pair(2,1))cycle     
     if( (iLform/=1) .and. (.not.compatible(Oi%form,Li%form,guess,ng,group)) )cycle
	 nbatom=findatom(ra(:,iatom)+Li%r(:,1)); if(nbatom<=0)cycle
     ndim=ndim+1
  end do; end do; end do

  print*,'number of molecules =', ndim
  allocate(Mformfunc(ndim))
  print*,ndim 
  idim=0
  do iatom=1,natom; do iOform=1,nOform; do iLform=1,nLform
     Oi=>Oformfunc(iOform); Li=>Lformfunc(iLform)
     if(SkipInterOrbitPair.and.Oi%pair(1,1)/=Oi%pair(2,1))cycle     
     if( (iLform/=1) .and. (.not.compatible(Oi%form,Li%form,guess,ng,group)) )cycle
	 nbatom=findatom(ra(:,iatom)+Li%r(:,1)); if(nbatom<=0)cycle
     idim=idim+1; Mi=>Mformfunc(idim)
     Mi%atom=iatom; Mi%nbatom=nbatom
	 Mi%Lform=iLform; Mi%Oform=iOform
	 if(iLform==1)then
	   Mi%form=Oi%form
	 else	 
	   Mi%form=guess; if(guess=='***')Mi%form=combine(Li%form,Oi%form,group)
	 end if
  end do; end do; end do
  print*,'beforecheck'
  call checkMform()
  print*,'check ok'
  !save details of molecules
  open(10,file='Mform.txt')
  do idim=1,ndim; Mi=>Mformfunc(idim)
     write(10,*)idim,'-th molecule:'
     write(10,*)'atom,nbatom=',Mi%atom,Mi%nbatom
	 Li=>Lformfunc(Mi%Lform); Oi=>Oformfunc(Mi%Oform)
	 write(10,*)'M-',Mi%form,' L-',Li%form,' O-',Oi%form,Oi%parity
	 write(10,*)'nr in Lform=',Li%nr
	 if(Li%ndim==2)then
	   write(10,*)'    rx        ry        f(r) '
	 else if(Li%ndim==3)then
	   write(10,*)'    rx        ry       rz        f(r) '
     end if
	 do i=1,Li%nr; write(10,100)Li%r(:,i),Li%f(i); end do
	 write(10,*)'np in Oform=',Oi%npair
	 write(10,*)'    O1    O2    f(O1,O2) '
	 do i=1,Oi%npair; write(10,200)Oi%pair(:,i),Oi%f(i); end do
	 write(10,*)'      '
  end do
  close (10)
100 format(1x,4f9.4)
200 format(1x,2i6,f10.4)

  print*, 'Mform saved to file.'

  return
end subroutine setMolecules_square


subroutine setMolecules_LowSymmetry()
  use workplace
  implicit none

  integer idim,iLform,iOform,i,i1,i2,iatom,nbatom
  integer findatom
  logical compatible
  type (Mformconfig), pointer :: Mi
  type (Oformconfig), pointer :: Oi
  type (Lformconfig), pointer :: Li

  real*8, pointer, dimension (:,:) :: ra
  real*8 :: rbatom(2)

  !notice: In low symmetry groups C2v, D2h and C2, all representations are one-dimensional: A1g, B2g, E1u, E2u. 

  ra=>model%ra
  ndim=0
  
  do iatom=1,natom; do iOform=1,nOform; do iLform=1,nLform
    Oi=>Oformfunc(iOform); Li=>Lformfunc(iLform)
    if(SkipInterOrbitPair.and.Oi%pair(1,1)/=Oi%pair(2,1))cycle  
    rbatom=ra(:,iatom)+Li%r(:,1)
	nbatom=findatom(rbatom)
	 if(nbatom<=0)cycle
	 if( model%sample=='CuO2' .and. (iatom/=Oi%pair(1,1).or.nbatom/=Oi%pair(2,1)) )cycle
     ndim=ndim+1
  end do; end do; end do

  print*,'number of molecules =', ndim
  allocate(Mformfunc(ndim))
 
  idim=0
  do iatom=1,natom; do iOform=1,nOform; do iLform=1,nLform
     Oi=>Oformfunc(iOform); Li=>Lformfunc(iLform)
     if(SkipInterOrbitPair.and.Oi%pair(1,1)/=Oi%pair(2,1))cycle     
	 nbatom=findatom(ra(:,iatom)+Li%r(:,1)); if(nbatom<=0)cycle
	 if( model%sample=='CuO2' .and. (iatom/=Oi%pair(1,1).or.nbatom/=Oi%pair(2,1)) )cycle
     idim=idim+1; Mi=>Mformfunc(idim)
     Mi%atom=iatom; Mi%nbatom=nbatom
	 Mi%Lform=iLform; Mi%Oform=iOform
	 Mi%form='***'
  end do; end do; end do

  call checkMform()

  !save details of molecules
  open(10,file='Mform.txt')
  do idim=1,ndim; Mi=>Mformfunc(idim)
     write(10,*)idim,'-th molecule:'
     write(10,*)'atom,nbatom=',Mi%atom,Mi%nbatom
	 Li=>Lformfunc(Mi%Lform); Oi=>Oformfunc(Mi%Oform)
	 write(10,*)'M-',Mi%form,' L-',Li%form,' O-',Oi%form,Oi%parity
	 write(10,*)'nr in Lform=',Li%nr
	 if(Li%ndim==2)then
	   write(10,*)'    rx        ry        f(r) '
	 else if(Li%ndim==3)then
	   write(10,*)'    rx        ry       rz        f(r) '
     end if
	 do i=1,Li%nr; write(10,100)Li%r(:,i),Li%f(i); end do
	 write(10,*)'np in Oform=',Oi%npair
	 write(10,*)'    O1    O2    f(O1,O2) '
	 do i=1,Oi%npair; write(10,200)Oi%pair(:,i),Oi%f(i); end do
	 write(10,*)'      '
  end do
  close (10)
100 format(1x,4f9.4)
200 format(1x,2i6,f10.4)

  print*, 'Mform saved to file.'

  return
end subroutine setMolecules_LowSymmetry


subroutine CheckMForm()
    use workplace
    implicit none

    real*8, dimension (3) :: ri,rj

	integer idim,iOform,iLform,npi,nbi,ip,i
	integer jdim,jOform,jLform,npj,nbj,jp,j
    integer iorb1,iorb2,jorb1,jorb2

    type (Mformconfig), pointer :: Mi,Mj
	type (Lformconfig), pointer :: Li,Lj
	type (Oformconfig), pointer :: Oi,Oj

    real*8 fOi,fLi,fOj,fLj,fi(ndim)
	logical mismatch
	integer ig
	integer indexMg(8)
	real*8 XMg(8)

	do jdim=1,ndim; Mj=>Mformfunc(jdim)

	   fi=0

	   jOform=Mj%Oform; jLform=Mj%Lform
	   Oj=>Oformfunc(jOform); Lj=>Lformfunc(jLform)
	   npj=Oj%npair; do jp=1,npj; fOj=Oj%f(jp)
       jorb1=Oj%pair(1,jp); jorb2=Oj%pair(2,jp)
	   if(natom>1)call modifyindex(jorb1,Mj%atom,jorb2,Mj%nbatom)
	   nbj=Lj%nr; do j=1,nbj; fLj=Lj%f(j); rj=Lj%r(:,j)
	   
       do idim=1,ndim; Mi=>Mformfunc(idim)
	   iOform=Mi%Oform; iLform=Mi%Lform
	   Oi=>Oformfunc(iOform); Li=>Lformfunc(iLform)
	   npi=Oi%npair; do ip=1,npi; fOi=Oi%f(ip)
       iorb1=Oi%pair(1,ip); iorb2=Oi%pair(2,ip)
	   if(natom>1)call modifyindex(iorb1,Mi%atom,iorb2,Mi%nbatom)	   
	   if( mismatch((/jorb1,jorb2/),(/iorb1,iorb2/) ) )cycle
	   nbi=Li%nr; do i=1,nbi; fLi=Li%f(i); ri=Li%r(:,i)
	   if(sum( abs(ri-rj) )>1.e-5)cycle
  
          fi(idim) = fi(idim) + fOi*fLi*fOj*fLj

	   end do; end do; end do; end do; end do

	   do idim=1,ndim
	      if(idim==jdim.and.abs(fi(idim)-1)>1.e-5)pause 'normalization error @ checkMform'
		  if(idim/=jdim.and.abs(fi(idim))>1.e-5)pause 'orthogonality error @ checkMform'
	   end do
	end do
	print*,'Mform checked to be orthonormal'

    return
end subroutine CheckMForm
