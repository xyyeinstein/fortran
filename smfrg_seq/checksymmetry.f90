subroutine checksymmetry(Q)
  use workplace
  implicit none
  real*8 Q

  complex*16 hk(norb,norb),hkg(norb,norb)
  complex*16, dimension (norb) :: ci,cg,cj
  real*8 eval(norb),evalj(norb)

  integer ig,i,itest,ntest,iorb
  real*8 vk(2),pv(3),pg(3)
  complex*16 overlap
  real*8 ran
  integer iseed
  logical symmetric,propergauge

  print*,'checking symmetry in EigenStates...'
  
  propergauge=model%propergauge
  model%propergauge=.false.
  symmetric=.true.
  ntest=100; iseed=29735
  
  do itest=1,ntest
    vk(1)=2*ran(iseed)-1; vk(2)=2*ran(iseed)-1
    
	call gethk(norb,vk,Q,hk,model,MF,.false.)
    !print*,itest; pause
	call ZHEIGEN(norb,hk,evalj)

    do i=1,norb
	   ci=hk(:,i); pv=0; pv(1:2)=vk
	   do ig=2,ng
  	      call BlochImage(ig,norb,ci,cg,model)

		  call groupaction(ig,pv,pg,group)
		  call gethk(norb,pg(1:2),Q,hkg,model,MF,.false.)
		  call ZHEIGEN(norb,hkg,eval); cj=hkg(:,i)
		  !print*,eval; pause 'eval ok?'

		  overlap=sum( conjg(cj)*cg )
		  if(abs(abs(overlap)-1)>1.e-5)then
		    print*,i
			print*,'pv:(',pv,')'
			print*,ig
			print*,'pg:(',pg,')'
			print*,'     ci                    cg                     cj'
            !print*,ci,cg,cj
			do iorb=1,norb
			  write(*,'(2e10.3,2x,2e10.3,2x,2e10.3)')ci(iorb),cg(iorb),cj(iorb)
			end do
			print*,abs(overlap)
			print*,'     eval        evalj'
			do iorb=1,norb
			  write(*,'(1e10.3,2x,1e10.3)')eval(iorb),evalj(iorb)
			end do

			pause 'overlap error'
		  end if
	   end do
	end do
  end do
  if(symmetric)then; print*,'Eigenstates satisfy point group symmetry @ checksymmetry'
  else; stop 'Eigenstate do not respect point group symmetry @ checksymmetry'; end if

  model%propergauge=propergauge
  return
end subroutine checksymmetry

subroutine checkHsymmetry(Q)
  use workplace
  implicit none
  real*8 Q

  complex*16 hk(norb,norb),hkg(norb,norb)
  complex*16, dimension (norb) :: ci,cg,cj
  real*8 eval(norb),evalj(norb)

  integer ig,i,itest,ntest,atomimage
  integer iatom,iorb,idim,jatom,jorb,jdim
  integer iatomg,iorbg,idimg,jatomg,jorbg,jdimg
  integer indexgi(2),indexgj(2)
  real*8 Xgi(2),Xgj(2)
  real*8 vk(2),pv(3),pg(3)
  complex*16 overlap
  real*8 ran
  integer iseed
  logical symmetric,propergauge

  print*,'checking symmetry in Hk...'
  
  propergauge=model%propergauge
  model%propergauge=.false.
  symmetric=.true.
  ntest=100; iseed=29735
  norb1atom=norb/model%natom

  do itest=1,ntest
    vk(1)=2*ran(iseed)-1; vk(2)=2*ran(iseed)-1
    call gethk(norb,vk,Q,hk,model,MF,.false.)
    pv=(/vk(1),vk(2),Q/)
	do ig=2,ng
	  call groupaction(ig,pv,pg,group)
      call gethk(norb,pg(1:2),Q,hkg,model,MF,.false.)
	  do iatom=1,model%natom
	    iatomg=atomimage(ig,iatom,model)
		do iorb=1,norb1atom
          idim=iorb+(iatom-1)*norb1atom
		  call orbitimage(ig,iorb,norb1atom,model%orbits,model%group,indexgi,Xgi)
		  iorbg=indexgi(1)
		  idimg=iorbg+(iatomg-1)*norb1atom
	      do jatom=1,model%natom
	        jatomg=atomimage(ig,jatom,model)
		    do jorb=1,norb1atom
              jdim=jorb+(jatom-1)*norb1atom
		      call orbitimage(ig,jorb,norb1atom,model%orbits(1:5),model%group,indexgj,Xgj)
		      jorbg=indexgj(1)
		      jdimg=jorbg+(jatomg-1)*norb1atom
              
			  if(abs(hkg(idimg,jdimg)*Xgi(1)*Xgj(1)-hk(idim,jdim))>1.e-5)then
			    print*,'       pv       ,   iatom ,iorb  ,idim  ,jatom ,jorb  ,jdim'
				write(*,'(1f5.2,1x,1f5.2,1x,1f5.2,6i7)')pv,iatom,iorb,idim,jatom,jorb,jdim
                print*,'ig=',ig
                print*,'       pg       ,   iatomg,iorbg ,idimg ,jatomg,jorbg ,jdimg'
				write(*,'(1f5.2,1x,1f5.2,1x,1f5.2,6i7)')pg,iatomg,iorbg,idimg,jatomg,jorbg,jdimg
				print*,'                            Xgi                  Xgj     '
				write(*,'(28x,1f4.1,18x,1f4.1)')Xgi(1),Xgj(1)
                print*,' hkg(idim,jdim)    hk(idimg,jdimg)'
				write(*,'(2f8.5,4x,2f8.5)')hkg(idimg,jdimg),hk(idim,jdim)
				pause
			  end if
			end do
		  end do
		end do
	  end do
	end do
  end do

  if(symmetric)then; print*,'Hk satisfies point group symmetry @ checksymmetry'
  else; stop 'Hk does not respect point group symmetry @ checksymmetry'; end if

  model%propergauge=propergauge
  return
end subroutine





