subroutine definejob()
  use workplace
  implicit none

  !call Tmatrix_(); stop

  one=cmplx(0,1)
  pi=asin(1.d0)*2
  twopi=pi*2

  !call checkla(); stop
  !call ChenWei(); stop

  !default parameters
  model%n4logqmesh=1
  model%U=0; model%JH=0; model%Uab=0
  model%Vnn=0; model%Jnn=0
  model%periodiclayers=.false.
  model%nest=0
  model%a=(/1,0,0/); model%b=(/0,1,0/); model%c=(/0,0,1/)
  appendfs=.false.  
  
  !call compass()
  !call cuprate2L()
  !call Dagotto()
  !call KFeSE1L10d()
  !call KFeSe1L5d()
  !call KFeSe1Lav()
  !call FeSe1Fe()
  !call FS1FLee()
  !call FS2F()
  !call KFeSe2L20d()
  !call KFeSe2vL()
  !call KFeSePrimitive()
  !call Kuroki()
  !call Raghu()
  call ehm1d()
  !call openladder()
  !call rice()
  !call square1atom1orb()
  !call square2atom2orb()
  !call WangFa5d()
  !call BiS2()
  !call BiS2_1band()
  !call SrRuO1orb()
  !call SrRuO2orb()
  !call SrRuO3orb()
  !call CuO2()
  !call kagome()
  !call LiFeAs1L10d()
  !call LiFeAs1L5d()
  !call LiFeAs2L5d()
  !call hujp()
  !call hujp2()
  !call lhq()
  !call BFA1L10d()
  !call BaFeAs1L5d()
  
  call writeinputs()

  return
end subroutine definejob



subroutine WriteInputs()
  use workplace
  implicit none

  open(10,file='inputed.txt')

  write(10,*)'Sample = ',model%sample
  write(10,*)'DataReady=',DataReady
  write(10,*)'norb,natom=',model%norb,model%natom
  write(10,*)'atomic positions:'
  write(10,*)model%ra
  write(10,*)'sample=',model%sample
  write(10,*)'orbits=',model%orbits
  write(10,*)'zerotz=',model%zerotz
  write(10,*)'U,JH,Uab,Vnn,mu='
  write(10,*)model%U,model%JH,model%Uab,model%Vnn,model%mu
  if(model%sample=='CuO2')then
    write(10,*)'Udd,Upp,Vdd,Vpp,Vpd='
	write(10,*)model%Udd,model%Upp,model%Vdd,model%Vpp,model%Vpd
  end if
  if(model%ephcoupled)write(10,*)'Veph=',model%Veph
  write(10,*)'Appendkmesh=',appendkmesh
  write(10,*)'kmeshfile=',kmeshfile
  write(10,*)'usegroup,group,ng=',usegroup,group,ng
  write(10,*)'nkx,nky=',nkx,nky
  write(10,*)'nLform,Lcontact=',nLform,Lcontact
  write(10,*)'Lformtable=',Lformtable
  write(10,*)'Lformbasis='
  write(10,*)Lformbasis
  write(10,*)'SkipInterOrbitPair=',SkipInterOrbitPair
  write(10,*)'guessed pairing symmetry=',guess
  write(10,*)'nw=',nw
  write(10,*)'wmax,wmin=',wmax,wmin
  write(10,*)'wir,diverge=',wir,diverge

  close(10)

  return
end subroutine WriteInputs



  
	     
