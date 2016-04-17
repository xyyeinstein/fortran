subroutine definejob()
  use workplace
  implicit none

  integer imodel
  
  bcsflow=.false.
  dataready=.false.
  quicksearch=.false.
  appendbreakpoint=.false.
  appendkmesh=.false.
  appendprojectors=.false.
  useprojectors=.false.
  nw=400
  wmax=100
  wmin=1.d-4
  wir=1.d-5
  diverge=40
  model%u=2.0
  model%vnn=0.0
  model%t1=0
  model%mu=0
  model%vl2l=0
  model%vbonding=0
  model%jnn=0
  model%la=0.01
  imodel=1
  select case (imodel)
    case(1)
      call cuprate()
    case(2)
      call cuprate2L()
    case default
  end select

  call writeinputs()
  return
end subroutine definejob


subroutine WriteInputs()
  use workplace
  implicit none

  open(10,file='inputed.txt')

  write(10,*)model%sample
  write(10,*)'DataReady=',DataReady
  write(10,*)'norb,natom=',model%norb,model%natom
  write(10,*)'atomic positions:'
  write(10,*)model%ra
  write(10,*)'U,mu='
  write(10,*)model%U,model%mu
  write(10,*)'Appendkmesh=',appendkmesh
  write(10,*)'kmeshfile=',kmeshfile
  write(10,*)'usegroup,group,ng=',usegroup,group,ng
  write(10,*)'nLform,Lcontact=',nLform,Lcontact
  write(10,*)'Lformtable=',Lformtable
  write(10,*)'Lformbasis='
  write(10,*)Lformbasis
  write(10,*)'nw=',nw
  write(10,*)'wmax,wmin=',wmax,wmin
  write(10,*)'wir,diverge=',wir,diverge

  close(10)

  return
end subroutine WriteInputs








  
	     
