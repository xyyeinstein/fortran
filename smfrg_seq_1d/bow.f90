program main
implicit none

integer phase,stat1
real*8 U,V,mu,la

open(10,file='phd.dat')
open(11,file='bow.dat')
do while(.true.)
	read(10,'(4f10.5,1i4)',iostat=stat1)U,V,mu,la,phase
	if(phase==5)then
		write(11,'(2f10.5)')U,V
	endif
    if(stat1/=0)exit
end do
close(10)
close(11)
stop
end program
