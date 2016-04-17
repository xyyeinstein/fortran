subroutine compass()
  use workplace
  implicit none

  DataReady=.false.
  AppendBreakPoint=.false.
  AppendKmesh=.true.; kmeshfile='cmpsskmesh.dat'

  norb=4; natom=4
  model%norb=4; model%natom=4; model%sample='cmpss'
  print*,'The sample chosen is ',model%sample;

  allocate(model%ra(3,4),model%orbits(4))
  model%ra(:,1)=(/0.25,0.,0./); model%ra(:,2)=(/0.,0.25,0./)
  model%ra(:,3)=(/-0.25,0.,0./); model%ra(:,4)=(/0.,-0.25,0./)
  model%a=(/1,0,0/); model%b=(/0,1,0/); model%c=(/0,0,1/)
  model%orbits='A1g'
  model%periodiclayers=.false.
   
  usegroup=.true.;  group='C4v';  ng=8
  model%group=group; model%ng=ng

  model%U=4; model%Vnn=0; model%mu=0.

  model%propergauge=.false.
  if(.not.model%propergauge)call checksymmetry(0.d0)
  nkx=40; nky=40   
  nw=400; wmax=1.e2
  wmin=1.e-4; wir=1.e-4
  diverge=30

  nLform=1; Lcontact=0
  allocate(Lformtable(nLform))
  allocate(Lformbasis(3,nLform))

  Lformtable='A1g'
  Lformbasis=0

  call setLformfunc(ng)
  call setOformfunc()
  call setmolecules()

  return
end subroutine compass


subroutine compasshk(norb,kv,hk,model,mf,slvbsn)
  use standard_derived_types
  implicit none

  integer norb
  real*8 kv(2)
  complex*16 hk(norb,norb)
  type (modelblock) :: model
  type (mfblock) :: mf
	logical slvbsn

  real*8 pi,bond(2),dr(2),t1
  complex*16 one
  integer i,j

  pi=asin(1.d0)*2; one=cmplx(0,1)
  t1=1

  hk=0
 
  do i=1,norb; hk(i,i)=-model%mu; end do
  
  !1-2 NN-bond
  bond=(/-0.25,0.25/)
  dr=0; if(model%propergauge)dr=model%ra(1:2,2)-model%ra(1:2,1)  
  hk(1,2)=hk(1,2)-t1*exp(one*pi*sum((bond-dr)*kv))

  !1-3 NN-bond
  bond=(/0.5,0./)
  dr=0; if(model%propergauge)dr=model%ra(1:2,3)-model%ra(1:2,1)  
  hk(1,3)=hk(1,3)-t1*exp(one*pi*sum((bond-dr)*kv))

  !1-4 NN-bond
  bond=(/-0.25,-0.25/)
  dr=0; if(model%propergauge)dr=model%ra(1:2,4)-model%ra(1:2,1)  
  hk(1,4)=hk(1,4)-t1*exp(one*pi*sum((bond-dr)*kv))

  !2-3 NN-bond
  bond=(/-0.25,-0.25/)
  dr=0; if(model%propergauge)dr=model%ra(1:2,3)-model%ra(1:2,2)  
  hk(2,3)=hk(2,3)-t1*exp(one*pi*sum((bond-dr)*kv))

  !2-4 NN-bond
  bond=(/0.,0.5/)
  dr=0; if(model%propergauge)dr=model%ra(1:2,4)-model%ra(1:2,2)  
  hk(2,4)=hk(2,4)-t1*exp(one*pi*sum((bond-dr)*kv))

  !3-4 NN-bond
  bond=(/0.25,-0.25/)
  dr=0; if(model%propergauge)dr=model%ra(1:2,4)-model%ra(1:2,3)  
  hk(3,4)=hk(3,4)-t1*exp(one*pi*sum((bond-dr)*kv))

  do j=1,4; do i=1,j-1; hk(j,i)=conjg(hk(i,j)); end do; end do

  return
end subroutine compasshk
