subroutine runfrg()
  use workplace
  implicit none

  integer error
  real*8 b,w,La,dLa
  real*8 dopinglevel
  complex*16 work(ndim)
  real*8 q(2)
  integer icount,imode
  logical isave,iXq
  character*20 explain

  type (leadingchannel), pointer :: pair

  if(DataReady.and.(.not.AppendBreakPoint))then
    call GapForm_Band(work,0)               !work is a dummy argument here
    call Tmatrix()
    !call realspacetexture('sdw',work,q,0)  !work and q are dummy arguments here
    !call realspacetexture('cdw',work,q,0)  !work and q are dummy arguments here
    !call RPAXq(norb*2,5.d-1,5.d-1,1.d0,.true.)	
    !call bubblesusceptibility()
    return
  end if
  
  b=exp(log(wmax/wmin)/nw)
  open(30,file='flow.dat')
  open(31,file='flow.txt')
  open(22,file='Spp.dat')
  open(23,file='Spp.txt')

  w=wmax; icount=0
  if(AppendBreakPoint)then
    call breakpoint(w)
    do while (.true.); read(30,*,iostat=error)explain; if(error/=0)exit; end do
    do while (.true.); read(31,*,iostat=error)explain; if(error/=0)exit; end do
    do while (.true.); read(22,*,iostat=error)explain; if(error/=0)exit; end do
    do while (.true.); read(23,*,iostat=error)explain; if(error/=0)exit; end do
  end if

  La=w*(1+1/b)*0.5

  do while(w>=wir)
    icount=icount+1
    La=w*(1+1/b)*0.5
    dLa=w*(1-1/b)
    if(.not.model%ladder)call flexiblemesh(La)
    if(w<wirx0)useX0=.false.
    call singularflow(La,dLa,icount)
      
    if(useX0.and.(.not.AppendX0))then
      w=w/b
      cycle
    end if
    
    isave=(mod(icount,40)==0)
    iXq=.false.
    if(QuickSearch)isave=.false.
    call analysis(La,isave,iXq)
    !if(isave)call outputV(w)     !this is used if breakpoint in the intermediate stages is necessary	 
    !if(isave.and.(.not.BCSflow))call runBCS_k(w,.false.,.true.)

    write(*,100)w,output%Vmax 
    write(*,200)output%q; print*,' '
    write(31,100)w,output%Vmax 
    write(31,200)output%q
    write(30,100)w,output%Vmax 

    if(ndim>1)then
      do imode=1,10
        pair=>output%leadingpair(imode)
        write(22,220)pair%s,pair%q
      end do
      write(23,230)w,output%leadingpair(1)%s,output%leadingpair(1)%form,output%leadingpair(2)%s,output%leadingpair(2)%form
    end if

    !if(output%Vmax(1)>diverge)BCSflow=.true.
    if(output%Vmax(1)>diverge)exit
    w=w/b
  end do
  close(30); close(31); close(22); close(23)
100 format(1x,5e15.6)
200 format(1x,26(' '),3('(',2f7.3,')  ') )

220 format(1x,3e20.6)
230 format(1x,f15.6,2(f15.6,a6))

  if(useX0.and.(.not.AppendX0))return

  if(QuickSearch)call FinalProjection()

  call analysis(La,.true.,.false.)
  call outputV(w)
  !call EigenBCS_R(3)
  !call runBCS_k(w,.false.,.true.)

  return
end subroutine runfrg





subroutine OutputV(w)
  use workplace
  real*8 w
  integer idim,iq
  real*8 q(2),scale
  logical hexagonal
  integer findq
  complex*16 work(ndim,ndim)
  real*8 eval(ndim)

  open(10,file='Pair.dat')
  open(11,file='SDW.dat')
  open(12,file='CDW.dat')
  do iq=1,nq
     work=P(:,:,iq); call ZHEIGEN(ndim,work,eval); write(10,*)eval(ndim)
     work=-C(:,:,iq); call ZHEIGEN(ndim,work,eval); write(11,*)eval(ndim)
     work=2*D(:,:,iq)-C(:,:,iq); call ZHEIGEN(ndim,work,eval); write(12,*)eval(ndim)
  end do
100 format(1x,4f15.6)
  close(10); close(11); close(12)

  open(10,file='Vout.dat',form='UNFORMATTED')
  write(10)w,P,C,D
  close(10)  
  return
end subroutine OutputV

subroutine breakpoint(w)
  use workplace
  real*8 w

  open(10,file='Vout.dat',form='UNFORMATTED')
  read(10)w,P,C,D
  close(10)  
  return
end subroutine breakpoint

subroutine FinalProjection()
  use workplace
  implicit none

  complex*16, dimension (ndim,ndim,nq) :: dP,dC,dD,Pwork,Cwork,Dwork,work

  Pwork=P; Cwork=C; Dwork=D

  dP=0; dC=Cwork; dD=Dwork
  QuickSearch=.true.; call contact(dP,dC,dD); work=Pwork-dP
  dP=0; dC=Cwork; dD=Dwork
  QuickSearch=.false.; call contact(dP,dC,dD); P=work+dP

  dP=Pwork; dC=0; dD=Dwork
  QuickSearch=.true.; call contact(dP,dC,dD); work=Cwork-dC
  dP=Pwork; dC=0; dD=Dwork
  QuickSearch=.false.; call contact(dP,dC,dD); C=work+dC

  dP=Pwork; dC=Cwork; dD=0
  QuickSearch=.true.; call contact(dP,dC,dD); work=Dwork-dD
  dP=Pwork; dC=Cwork; dD=0
  QuickSearch=.false.; call contact(dP,dC,dD); D=work+dD

  return
end subroutine FinalProjection
