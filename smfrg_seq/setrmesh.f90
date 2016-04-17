subroutine setRmesh()
  use workplace
  implicit none

  integer iR,ix,iy,ig,iRg,iag,jag,iatom,jatom
  real*8 rv(3),rg(3),ra(3,natom)
  integer findR,atomimage
  
  iR=0
  do iy=-Lcontact,Lcontact; do ix=-Lcontact,Lcontact
     rv=ix*model%a+iy*model%b
     iR=iR+1
  end do; end do

  nR=iR; allocate(R(2,nR))

  iR=0
  do iy=-Lcontact,Lcontact; do ix=-Lcontact,Lcontact
     rv=ix*model%a+iy*model%b
     iR=iR+1; R(:,iR)=(/ix,iy/)
  end do; end do

  do iR=1,nR; if(findR(R(:,iR))<0)stop 'R is not searchable @ setRmesh'; end do

  if(.not.usegroup)return
  
  allocate(indexRg(3,ng,nR,natom,natom))
  indexRg=-1; ra=model%ra
  do iR=1,nR; do iatom=1,natom; do jatom=1,natom
     rv=R(1,iR)*model%a+R(2,iR)*model%b+ra(:,jatom)-ra(:,iatom)
	 do ig=1,ng; call groupaction(ig,rv,rg,group)
	    iag=atomimage(ig,iatom,model); jag=atomimage(ig,jatom,model)
		rg=rg-ra(:,jag)+ra(:,iag)
	    ix=nint( sum(rg(1:2)*model%ka) / 2 )
		iy=nint( sum(rg(1:2)*model%kb) / 2 )
		if( sum( abs(rg-ix*model%a-iy*model%b) )>1.e-5 )stop 'Rg not truncatable @ setRmesh'
		iRg=findR( (/ix,iy/) )
		indexRg(:,ig,iR,iatom,jatom)=(/iRg,iag,jag/)
	 end do
  end do; end do; end do

  return
end subroutine setRmesh

function findR(Rc2c)
  use workplace, only : nR,R
  implicit none
  integer findR,Rc2c(2)
  integer iR,iU,iL,info

  findR=-1
  
  info=compare(Rc2c,R(:,1))
  select case (info)
    case (-1); return
	case (0); findR=1; return
  end select

  info=compare(Rc2c,R(:,nR))
  select case (info)
    case (0); findR=nR; return
    case (1); return
  end select

  iU=nR; iL=1
  do while (iU-iL>1)
     iR=(iU+iL)/2
     info=compare(Rc2c,R(:,iR))
	 select case (info)
	   case (-1); iU=iR
	   case (0); findR=iR; return
	   case (1); iL=iR
	 end select
  end do

  return
  contains
    function compare(R1,R2)
	  integer compare,R1(2),R2(2)
	  compare=0
	  if(R1(2)>R2(2))then
	    compare=1
	  else if(R1(2)<R2(2))then
	    compare=-1
	  else
	    if(R1(1)>R2(1))then
		  compare=1
		else if(R1(1)<R2(1))then
		  compare=-1
		else
		  compare=0
		end if
	  end if
	  return
	end function compare
end function findR
