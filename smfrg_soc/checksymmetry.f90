subroutine checkHsymmetry()
  use workplace, only : natom,norb,ng,model
  implicit none

  complex*16 hk(norb,norb),hkg(norb,norb)
  complex*16, dimension (norb) :: ug
  real*8 ek(norb),ekg(norb)

  integer ig,i,itest,ntest
  real*8 kv(2),kg(2),pv(3),pg(3)
  complex*16 overlap
  real*8 ran
  integer iseed
  logical symmetric

  print*,'checking symmetry in Hk...'

  symmetric=.true.

  ntest=1000; iseed=29735

  do itest=1,ntest
     kv(1)=ran(iseed)-0.5; kv(2)=ran(iseed)-0.5; kv=kv*4
     call gethk(norb,kv,hk,model); pause '100'
     call ZHEIGEN(norb,hk,ek)
	 do ig=2,ng; pv=0; pv(1:2)=kv
	    call groupaction(ig,pv,pg,model%group); kg=pg(1:2)
		call gethk(norb,kg,hkg,model)
		call ZHEIGEN(norb,hkg,ekg)
		do i=1,norb; if(abs(ekg(i)-ek(i))>1.e-5)stop 'ek not symmetric'
		   call BlochImage(ig,norb,hk(:,i),ug,model)
		   overlap=abs( sum( conjg(ug)*hkg(:,i) ) )
		   if(abs(overlap-1)>1.e-5)then
		     symmetric=.false.
		     print*,'Bloch state not symmetric'
			 print*,'|k>=';   print*,hk(:,i)
			 print*,'|kg>=';  print*,hkg(:,i)
			 print*,'g|k>=';  print*,ug
			 print*,'overlap=',overlap
			 pause 'continue?'
		   end if
		end do
	 end do
  end do
  if(symmetric)then; print*,'Hk satisfies point group symmetry @ checksymmetry'
  else; stop 'Hk does not respect point group symmetry @ checksymmetry'; end if
  return
end subroutine checkHsymmetry


subroutine checkTsymmetry()
  use workplace
  implicit none

  integer itest,ntest,iseed
  integer iband,ia,ispin,iispin,idim,iidim,isign
  complex*16 hk(norb,norb),hkbar(norb,norb),work(norb)
  real*8 ran,vk(2),ev(norb),evbar(norb),error

  ntest=100; iseed=29375
  do itest=1,ntest
     vk(1)=ran(iseed)-0.5; vk(2)=ran(iseed)-0.5; vk=vk*4
     call gethk(norb,vk,hk,model)
	 call ZHEIGEN(norb,hk,ev)
	 call gethk(norb,-vk,hkbar,model)
	 call ZHEIGEN(norb,hkbar,evbar)
     do iband=1,norb
	    if(abs(ev(iband)-evbar(iband))>1.e-5)stop 'eigen value error @ checkTsymmetry'
	    do ispin=1,2; iispin=2; if(ispin==2)iispin=1
           isign=1; if(ispin>iispin)isign=-1
           do ia=1,natom; idim=ia+(ispin-1)*natom; iidim=ia+(iispin-1)*natom
		      work(idim)=conjg(hk(iidim,iband))*isign
		end do; end do
		error=abs( abs(sum( conjg(hkbar(:,iband))*work )) -1 )
		if(error>1.e-5)stop 'Kramers pair error: |-k> /= s_y conjg(|k>) @ checkTsymmetry'
	 end do
  end do
  print*,'T symmetry obeyed @ checkTsymmetry'
  return
end subroutine checkTsymmetry


subroutine checkHsymmetry_()
  use workplace, only : natom,norb,ng,model
  implicit none

  complex*16 hk(norb,norb),hkg(norb,norb),work(norb,norb)

  integer ig,i,itest,ntest
  real*8 kv(2),kg(2),pv(3),pg(3)

  integer ispin,iatom,idim,iispin,iidim
  integer jspin,jatom,jdim,jjspin,jjdim
  integer spinindexg(2)
  complex*16 xph,spinxg(2)

  real*8 ran
  integer iseed
  logical symmetric

  print*,'checking symmetry in Hk...'

  symmetric=.true.
  ntest=1000; iseed=29735

  do itest=1,ntest
     kv(1)=ran(iseed)-0.5; kv(2)=ran(iseed)-0.5; kv=kv*4
     call gethk(norb,kv,hk,model)
	 do ig=2,ng; pv=0; pv(1:2)=kv
	    call groupaction(ig,pv,pg,model%group); kg=pg(1:2)
		call gethk(norb,kg,hkg,model)
	    call spinimage(ig,model%group,SpinIndexg,SpinXg)
        work=0
		do iatom=1,natom; do ispin=1,2; idim=iatom+(ispin-1)*natom
        iispin=spinindexg(ispin); iidim=iatom+(iispin-1)*natom
		do jatom=1,natom; do jspin=1,2; jdim=jatom+(jspin-1)*natom
		jjspin=spinindexg(jspin); jjdim=jatom+(jjspin-1)*natom
	       Xph=conjg(SpinXg(ispin))*SpinXg(jspin)  
	       Xph=conjg(Xph)   !phase factor for ci^dag cj
           work(idim,jdim)=work(idim,jdim)+hk(iidim,jjdim)*Xph
		end do; end do; end do; end do
		symmetric=(sum( abs( work-hkg) )<1.e-5)
		if(.not.symmetric)pause 'error: hk not symmetric under point group'
	 end do
  end do
  if(symmetric)print*,'hkg=g(hk) under point group.'
  return
end subroutine checkHsymmetry_
