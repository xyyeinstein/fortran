subroutine runfrg()
  use workplace
  implicit none

  real*8 b,w,La,dLa
  real*8 dopinglevel
  complex*16 work(ndim)
  real*8 q(2)
  character*20 explain
  integer icount
  logical isave
  
  integer error

  print*,'filling=',dopinglevel(); 

  if(DataReady.and.(.not.AppendBreakPoint))then
	call GapForm_Band(work,0)              !work is a dummy argument here
	!call BdG_Band(work,0)
	call BdGSlab(model,edgemodel)
	!call realspacetexture('phd',work,q,0)  !work and q are dummy arguments here
	stop
  end if

  b=exp(log(wmax/wmin)/nw)
  open(9,file='flow.dat')
  open(19,file='flow.txt')
  open(29,file='scmax.dat')
  open(21,file='fm.dat')
  w=wmax
  if(AppendBreakPoint)then
    call breakpoint(w)
	do while (.true.); read(9,*,iostat=error)explain; if(error/=0)exit; end do
	do while (.true.); read(19,*,iostat=error)explain; if(error/=0)exit; end do	
	do while (.true.); read(29,*,iostat=error)explain; if(error/=0)exit; end do	
	do while (.true.); read(21,*,iostat=error)explain; if(error/=0)exit; end do
  end if
  icount=0
  do while(w>=wir)
     icount=icount+1
     La=w*(1+1/b)*0.5;	 dLa=w*(1-1/b)
	 call flexiblemesh(La)
	 call singularflow(La,dLa)
	 isave=(mod(icount,50)==0)
  	 call analysis(La,isave)
	 
	 write(*,100)w,output%Vmax 
	 write(*,200)output%q; print*,' '
	 
	 write(9,100)w,output%Vmax 

	 write(19,100)w,output%Vmax 
	 write(19,200)output%q; print*,' '
   
   write(29,*)output%sc
   write(21,*)output%fm
	 
	 w=w/b
	 if(output%Vmax(1)>diverge)exit
  end do
  close(9); close(19); close(29); close(21)
100 format(1x,5e15.6)
200 format(1x,26(' '),2('(',2f7.3,')  ') )

  if(QuickSearch)call FinalProjection() 

  call analysis(La,.true.)
  call outputV(w)
  return
end subroutine runfrg

subroutine OutputV(w)
  use workplace
  real*8 w
  integer idim,iq
  real*8 q(2)

  open(10,file='Pair.dat')
  open(11,file='phDW.dat')
  do idim=1,ndim; do iq=1,nq
     write(10,100)qv(1:2,iq),P(idim,idim,iq)
	 write(11,100)qv(1:2,iq),D(idim,idim,iq)-C(idim,idim,iq)
  end do; end do
100 format(1x,4f15.6)

  close(10); close(11)

  open(10,file='PD.dat')
  do iy=-40,40; do ix=-40,40; q=(/ix,iy/)*1.4/40
     call search(q,iq,qmesh)
	 do idim=1,ndim
        write(10,200)q,P(idim,idim,iq),D(idim,idim,iq)-C(idim,idim,iq)
	 end do
  end do; end do
  close(10)
200 format(1x,6f12.5)


  open(10,file='Vout.dat')
  write(10,*)w
  write(10,*)P,C,D
  close(10)  
  return
end subroutine OutputV


subroutine breakpoint(w)
  use workplace
  real*8 w

  open(10,file='Vout.dat')
  read(10,*)w
  read(10,*)P,C,D
  close(10)  
  return
end subroutine breakpoint

subroutine FinalProjection()
  use workplace
  implicit none

  complex*16, dimension (ndim,ndim,nq) :: dP,dC,dD,Pwork,Cwork,Dwork,work

  Pwork=P; Cwork=C; Dwork=D

  dP=0; dC=Cwork; dD=Dwork
  QuickSearch=.true.;  call contact(dP,dC,dD); work=P-dP
  dP=0; dC=Cwork; dD=Dwork
  QuickSearch=.false.; call contact(dP,dC,dD); P=work+dP


  dP=Pwork; dC=0; dD=Dwork
  QuickSearch=.true.;  call contact(dP,dC,dD); work=C-dC
  dP=Pwork; dC=0; dD=Dwork
  QuickSearch=.false.;  call contact(dP,dC,dD); C=work+dC

  dP=Pwork; dC=Cwork; dD=0
  QuickSearch=.true.;  call contact(dP,dC,dD); work=D-dD
  dP=Pwork; dC=Cwork; dD=0
  QuickSearch=.false.;  call contact(dP,dC,dD); D=work+dD

  return
end subroutine FinalProjection


