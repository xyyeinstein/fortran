subroutine setOformfunc()
  use workplace

  if( sum( model%a*model%b ) ==0 ) then
    call setOformfunc_square()
  else
    call setOformfunc_LowSymmetry()
  end if

  return
end subroutine setOformfunc

subroutine setOformfunc_square()
  use workplace
  implicit none

  integer icount,iorb,jorb,iform,np
  integer iE(2)
  logical Eorbit
  character*3 combine
  type (Oformconfig), pointer :: O,Oi
  character*3, pointer, dimension (:) :: orbits

  integer morb,ig
  integer indexOg(4)
  real*8 XOg(4)

  morb=norb/natom
  
  if(model%sample=='rice')then; call setOformfunc_rice(); return; end if
  if(group=='C2v'.or.group=='D2h')then; call setOformfunc_LowSymmetry(); return; end if
  
  if(.true.)then
    nOform=morb*morb
	  allocate(Oformfunc(nOform))
    icount=0
	  do iorb=1,morb
	    do jorb=1,morb
	      icount=icount+1
		    O=>Oformfunc(icount)
		    O%parity='E'
		    O%npair=1; O%form='A1g'
	      allocate(O%pair(2,1),O%f(1))
	      O%pair(:,1)=(/iorb,jorb/); 
		    O%f=1
	    end do
	  end do
	  return
  end if

  orbits=>model%orbits

  iE=0;  icount=0
  do iorb=1,morb 
     if( Eorbit(orbits(iorb)) )then
	   icount=icount+1
	   iE(icount)=iorb
	 end if
  end do

  icount=0
  do jorb=1,morb; if(Eorbit( orbits(jorb) ) )cycle
  do iorb=1,jorb; if(Eorbit( orbits(iorb) ) )cycle
     icount=icount+1
  end do; end do
  nOform=icount

  if(iE(1)*iE(2)/=0)then
     nOform=icount+3
	 if(.not.skipOddOrbPair)nOform=icount+4
  end if

  allocate(Oformfunc(nOform))

  icount=0
  do jorb=1,morb; if(Eorbit( orbits(jorb) ) )cycle
  do iorb=1,jorb; if(Eorbit( orbits(iorb) ) )cycle
     icount=icount+1
     O=>Oformfunc(icount); O%parity='E'
	 if(iorb==jorb)then
	   O%npair=1; O%form='A1g'
	   allocate(O%pair(2,1),O%f(1))
	   O%pair=iorb; O%f=1
     else
	   O%npair=2; O%form=combine(orbits(iorb),orbits(jorb),group)
	   allocate(O%pair(2,2),O%f(2))
	   O%pair(:,1)=(/iorb,jorb/); O%pair(:,2)=(/jorb,iorb/)
	   O%f=1/sqrt(2.d0)
	 end if   
  end do; end do
  
  if(iE(1)*iE(2)/=0)then
    icount=icount+1
	O=>Oformfunc(icount); O%parity='E'
    O%npair=2
	allocate(O%pair(2,2),O%f(2))
	O%form='A1g'
	O%pair(:,1)=iE(1); O%pair(:,2)=iE(2)
	O%f=(/1,1/)/sqrt(2.d0)

    icount=icount+1
	O=>Oformfunc(icount); O%parity='E'
    O%npair=2
	allocate(O%pair(2,2),O%f(2))
	O%form='B1g'
	O%pair(:,1)=iE(1); O%pair(:,2)=iE(2)
	O%f=(/1,-1/)/sqrt(2.d0)

    icount=icount+1
	O=>Oformfunc(icount); O%parity='E'
    O%npair=2
	allocate(O%pair(2,2),O%f(2))
	O%form='B2g'
	O%pair(:,1)=iE; O%pair(:,2)=(/iE(2),iE(1)/)
	O%f=(/1,1/)/sqrt(2.d0)

    if(.not.skipOddOrbPair)then
      icount=icount+1
	  O=>Oformfunc(icount); O%parity='O'
      O%npair=2
	  allocate(O%pair(2,2),O%f(2))
	  O%form='A2g'
	  O%pair(:,1)=iE; O%pair(:,2)=(/iE(2),iE(1)/)
	  O%f=(/1,-1/)/sqrt(2.d0)
	end if
  end if

  return
end subroutine setOformfunc_square


function Eorbit(form)
  implicit none
  logical Eorbit
  character*3, intent (in) :: form
  Eorbit=(form=='E1u'.or.form=='E2u'.or.form=='E1g'.or.form=='E2g'.or.form=='E3u'.or.form=='A1u')
  return
end function Eorbit


subroutine setOformfunc_LowSymmetry()
  use workplace
  implicit none

  integer icount,iorb,jorb,iform,ig
  integer morb
  type (Oformconfig), pointer :: O,Oi
  character*3, pointer, dimension (:) :: orbits
  integer indexOg(4)
  real*8 XOg(4)
  character*3 formi,formj
  logical Eorbit
  character*1 parity

  !notice: In low symmetry groups C2v and D2h, all irreducible representations are one-dimensional: A1g, B2g, E1u, E2u. 

  orbits=>model%orbits
  morb=norb/natom; if(model%sample=='CuO2')morb=norb

  icount=0
  do jorb=1,morb; do iorb=1,morb
     formi=orbits(iorb); formj=orbits(jorb); parity='E'
	 if(Eorbit(formi).and.(.not.Eorbit(formj)))parity='O'
	 if(Eorbit(formj).and.(.not.Eorbit(formi)))parity='O'
	 if(skipOddOrbPair.and.parity=='O')cycle
     icount=icount+1
  end do; end do

  nOform=icount

  allocate(Oformfunc(nOform))

  icount=0
  do jorb=1,morb; do iorb=1,morb
     formi=orbits(iorb); formj=orbits(jorb); parity='E'
	 if(Eorbit(formi).and.(.not.Eorbit(formj)))parity='O'
	 if(Eorbit(formj).and.(.not.Eorbit(formi)))parity='O'
	 if(skipOddOrbPair.and.parity=='O')cycle
     icount=icount+1;  O=>Oformfunc(icount)
	 O%parity=parity;  O%npair=1; O%form='***'
	 allocate(O%pair(2,1),O%f(1))
	 O%pair(:,1)=(/iorb,jorb/); O%f=1
  end do; end do
  
  return
end subroutine setOformfunc_LowSymmetry

subroutine setOformfunc_rice()
  use workplace
  implicit none

  integer icount,iorb,jorb,iform,np,ig
  integer morb
  logical Eorbit
  type (Oformconfig), pointer :: O,Oi
  character*3, pointer, dimension (:) :: orbits
  integer OformImage

  nOform=1

  allocate(Oformfunc(nOform))

  O=>Oformfunc(1); O%parity='*'
  O%npair=1; O%form='A1g'
  allocate(O%pair(2,1),O%f(1))
  O%pair(:,1)=1; O%f=1
  
  return
end subroutine setOformfunc_rice
