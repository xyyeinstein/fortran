function formk(vk,formfunc)
  use workplace, only : pi,one
  use standard_derived_types
  implicit none

  real*8 vk(2)
  type (Lformconfig) :: formfunc
  complex*16 formk
  integer iform

  !notice: assuming the third dimension of vk, if any, is zero. 

  integer i,n
  real*8 fi,ri(2)

  formk=0
  n=formfunc%nr
  do i=1,n
     ri=formfunc%r(1:2,i)
     fi=formfunc%f(i)
     formk=formk+fi*exp( one*pi*sum(vk*ri) )
  end do

  return
end function formk

	    
subroutine setMolecules_()
  use workplace
  implicit none

  integer idim,iform,ia,ib,i,spin,nbspin
  integer findatom
  type (Mformconfig), pointer :: Mi
  type (Lformconfig), pointer :: Li

  real*8, pointer, dimension (:,:) :: ra
  real*8 rv(3)
  character*1 chr(2)

  ra=>model%ra; chr=(/'+','-'/)

  ndim=0
  do ia=1,natom; do iform=1,nLform
     rv=ra(:,ia)+Lformbasis(:,iform)
	 if(findatom(rv,model)>0)ndim=ndim+1
  end do; end do
  ndim=ndim*4

  allocate(Mformfunc(ndim))
 
  idim=0
  do ia=1,natom; do iform=1,nLform
     rv=ra(:,ia)+Lformbasis(:,iform)
	 ib=findatom(rv,model)
	 if(ib<=0)cycle
	 do spin=1,2; do nbspin=1,2; idim=idim+1; Mi=>Mformfunc(idim)
        Mi%atom=ia; Mi%nbatom=ib
		Mi%spin=spin; Mi%nbspin=nbspin; Mi%Sform=chr(spin)//chr(nbspin)
	    Mi%orb=ia+(spin-1)*natom; Mi%nborb=ib+(nbspin-1)*natom
	    Mi%Lform=iform;  Mi%form=Lformtable(iform)
        call formstructure(3,Lformbasis(:,iform),Mi%form,Mi%Lformfunc,ng,group)
     end do; end do
  end do; end do


  !save details of molecules
  open(10,file='Mform.txt')
  do idim=1,ndim; Mi=>Mformfunc(idim)
     write(10,*)idim,'-th molecule:'
     write(10,*)'atom,nbatom=',Mi%atom,Mi%nbatom
	 write(10,*)'orb,nborb=',Mi%orb,Mi%nborb
	 Li=>Mi%Lformfunc
	 write(10,*)'Sform: ',Mi%Sform
	 write(10,*)'Lform: ',Mi%form
	 write(10,*)'nr in Lform=',Li%nr
	 if(Li%ndim==2)then
	   write(10,*)'    rx        ry        f(r) '
	 else if(Li%ndim==3)then
	   write(10,*)'    rx        ry       rz        f(r) '
     end if
	 do i=1,Li%nr; write(10,100)Li%r(:,i),Li%f(i); end do
	 write(10,*)'      '
  end do
  close (10)
100 format(1x,5f9.4)

  print*, 'Mform saved to file.'

  call checkMform()
  return
end subroutine setMolecules_

subroutine CheckMForm()
    use workplace
    implicit none

    real*8, dimension (3) :: ri,rj

	integer idim,iLform,io,ino,ia,ina,nbi,i
	integer jdim,jLform,jo,jno,ja,jna,nbj,j

    type (Mformconfig), pointer :: Mi,Mj
	type (Lformconfig), pointer :: Li,Lj

    complex*16 fLi,fLj,fi(ndim)

	do jdim=1,ndim; Mj=>Mformfunc(jdim)

	   fi=0

	   jLform=Mj%Lform; Lj=>Mj%Lformfunc
	   nbj=Lj%nr; ja=Mj%atom; jna=Mj%nbatom
	   jo=Mj%orb; jno=Mj%nborb
	   do j=1,nbj; fLj=Lj%f(j); rj=Lj%r(:,j)
	   
       do idim=1,ndim; Mi=>Mformfunc(idim)
	   iLform=Mi%Lform; Li=>Mi%Lformfunc
	   nbi=Li%nr; ia=Mi%atom; ina=Mi%nbatom
	   io=Mi%orb; ino=Mi%nborb
	   if(ia/=ja.or.ina/=jna)cycle
	   if(io/=jo.or.ino/=jno)cycle

	   do i=1,nbi; fLi=Li%f(i); ri=Li%r(:,i)
       if(sum( abs(ri-rj) )>1.e-5)cycle
  
          fi(idim) = fi(idim) + conjg(fLi)*fLj

	   end do; end do; end do

	   do idim=1,ndim
	      if(idim==jdim.and.abs(fi(idim)-1)>1.e-5)pause 'normalization error @ checkMform'
		  if(idim/=jdim.and.abs(fi(idim))>1.e-5)pause 'orthogonality error @ checkMform'
	   end do
	end do
	print*,'Mform checked to be orthonormal'
    return
end subroutine CheckMForm
































subroutine setMolecules()
  use workplace
  implicit none

  character*3 guess   !pick out molecules of guessed symmetry and with nontrivial Lformfunc
  integer idim,iLform,iOform,i,i1,i2,iatom
  integer findatom
  logical compatible
  type (Mformconfig), pointer :: Mi
  type (Lformconfig), pointer :: Li

  real*8, pointer, dimension (:,:) :: ra

  ra=>model%ra

  ndim=nLform*natom
  
  allocate(Mformfunc(ndim))
 
  idim=0
  iLform=1

  idim=idim+1; Mi=>Mformfunc(idim)
  Mi%atom=1;   Mi%Lform=1;  Mi%form=Lformtable(iLform)
  call formstructure(3,Lformbasis(:,iLform),Mi%form,Mi%Lformfunc,ng,group)

  do iLform=2,nLform
     idim=idim+1; Mi=>Mformfunc(idim)
     Mi%atom=1; Mi%Lform=iLform; Mi%form=Lformtable(iLform)
     call formstructure(3,Lformbasis(:,iLform),Mi%form,Mi%Lformfunc,ng,group)
  end do

  if(natom>2)stop 'error: natom>2'
  do iatom=2,natom
     i1=1+idim*(iatom-1); i2=idim*iatom
     Mformfunc(i1:i2)%atom=iatom
     Mformfunc(i1:i2)%Lform=Mformfunc(1:idim)%Lform
     Mformfunc(i1:i2)%form=Mformfunc(1:idim)%form
	 do i=i1,i2; Mi=>Mformfunc(i)
	    call formstructure(3,-Lformbasis(:,Mi%Lform),Mi%form,Mi%Lformfunc,ng,group)
     end do
  end do

  do idim=1,ndim; Mi=>Mformfunc(idim)
     Mi%nbatom=findatom( ra(:,Mi%atom)+Mi%Lformfunc%r(:,1),model )
  end do

  call checkMform()

  !save details of molecules
  open(10,file='Mform.txt')
  do idim=1,ndim; Mi=>Mformfunc(idim)
     write(10,*)idim,'-th molecule:'
     write(10,*)'atom,nbatom=',Mi%atom,Mi%nbatom
	 Li=>Mi%Lformfunc
	 write(10,*)'L-',Mi%form
	 write(10,*)'nr in Lform=',Li%nr
	 if(Li%ndim==2)then
	   write(10,*)'    rx        ry        f(r) '
	 else if(Li%ndim==3)then
	   write(10,*)'    rx        ry       rz        f(r) '
     end if
	 do i=1,Li%nr; write(10,100)Li%r(:,i),Li%f(i); end do
	 write(10,*)'      '
  end do
  close (10)
100 format(1x,5f9.4)

  print*, 'Mform saved to file.'

  return
end subroutine setMolecules
