subroutine rpa(cutoff,appendChi0,chi0file)
  use workplace
  implicit none
  
  real*8 cutoff
  logical appendChi0
  character*10 chi0file
  
  integer iq,idim,imq,mq
  real*8 q(2)

  type (Mformconfig), pointer :: Mi
  complex*16, dimension (ndim,ndim,nq) :: dP,dC,dD
  complex*16, allocatable, dimension (:,:,:) :: rpaXc,rpaXs
  complex*16, dimension (ndim,ndim) :: Xpp,Xph,work,unity


  integer iqmax
  logical reducible
  real*8 eval(ndim),evmax

  integer findq
  real*8, allocatable :: qc(:,:)
  logical qcut

  open(10,file=chi0file,form='UNFORMATTED')
  unity=0; do idim=1,ndim; unity(idim,idim)=1; end do
  
  iqmax=nq-model%nest-1
  if(uniformqmesh)iqmax=nq

  qcut=.false.
  if(qcut)then
    mq=21
	  allocate(rpaXc(ndim,ndim,mq),rpaXs(ndim,ndim,mq),qc(2,mq))
    do imq=1,mq
	  	q=qc(:,imq)
	    print*,'rpa ...',imq
	  	print*,q
	    iq=findq(q,.false.)
      call bareSuscept(q,Xpp,Xph,cutoff)
	  	work = unity + matmul( C(:,:,iq), Xph ); call ZINVERT(ndim,work)
	    dC(:,:,iq) = matmul( work, C(:,:,iq) ) - C(:,:,iq)
	    rpaXs(:,:,imq) = matmul( Xph, work )
		end do 
  else 
    mq=nq
    allocate(rpaXc(ndim,ndim,nq),rpaXs(ndim,ndim,nq))
    do iq=1,nq; q=qv(1:2,iq); print*,'rpa ...',iq
      if(iq<=iqmax.and.usegroup.and.reducible(q,group))cycle
	    !q(1)=q(1)+2*pi
      if(.not.appendChi0)then
	    call bareSuscept(q,Xpp,Xph,cutoff)
	    write(10)Xpp,Xph
	  else; read(10)Xpp,Xph; end if
	 
	  work = unity - matmul( P(:,:,iq), Xpp ); call ZINVERT(ndim,work)    !P/(1-pX)
	  dP(:,:,iq) = matmul( work, P(:,:,iq) ) - P(:,:,iq)

	  work = unity + matmul( C(:,:,iq), Xph ); call ZINVERT(ndim,work)
	  dC(:,:,iq) = matmul( work, C(:,:,iq) ) - C(:,:,iq)
	  rpaXs(:,:,iq) = matmul( Xph, work )   

	  work = unity + matmul( C(:,:,iq) - 2* D(:,:,iq), Xph ); call ZINVERT(ndim,work)
	  dD(:,:,iq) = matmul( work, C(:,:,iq) - 2* D(:,:,iq) )  - ( C(:,:,iq) - 2* D(:,:,iq) )
	  rpaXc(:,:,iq) = matmul( Xph, work )    

	  dD(:,:,iq) = -0.5* ( dD(:,:,iq) - dC(:,:,iq) )
     
	  if((.not.usegroup).or.iq>iqmax )cycle
			  call symmetrizePCD_q(iq,dP,dC,dD)
	  call symmetrizePCD_q(iq,rpaXc,rpaXs,dD)  !dD is used to fix the empty position, not for any other purpose
    end do; close(10)
  end if
   
  open(11,file='rpaXc.dat')
  open(12,file='rpaXs.dat')
  do iq=1,mq
     work=rpaXc(:,:,iq)
	 call ZHEIGEN(ndim,work,eval)
	 write(11,*)eval(ndim)
	 if(eval(ndim)<0)print*,'cdw instability!'
	 work=rpaXs(:,:,iq)
	 call ZHEIGEN(ndim,work,eval)
	 write(12,*)eval(ndim)
	 if(eval(ndim)<0)print*,'sdw instability!'
  end do
  close(11)
  close(12)
    
  P=P+dP; C=C+dC; D=D+dD
  
  print*,'RPA sus OK!'

  call runBCS_k(cutoff,.false.,.true.)

  return
end subroutine rpa

subroutine bareSuscept(q,Xpp,Xph,cutoff)
  use workplace
  implicit none

  real*8 q(2),cutoff
  complex*16, dimension (ndim,ndim) :: Xpp,Xph

  integer ik,ib,jb,iorb1,jorb1,iorb2,jorb2,i,j,idim,iLform,iOform,npi,ip
  real*8 vk(2),sk,fpp,fph,fOi,Qc,lindhard
  complex*16 fki,formk
  complex*16, dimension (norb,norb) :: hkq,hk,hkk
  real*8, dimension (norb) :: evk,evkk,evkq

  type (Mformconfig), pointer :: Mi
  type (Oformconfig), pointer :: Oi
  complex*16, dimension (norb*norb,norb*norb) :: X,Y
  complex*16 fdag(ndim,norb*norb),f(norb*norb,ndim)

  Qc=0  !this refers to z-momentum, but is not used for 2d systems

  Xpp=0; Xph=0
  
  do ik=1,nk; vk=kv(1:2,ik); sk=kv(3,ik)
  
     fdag=0
     do idim=1,ndim; Mi=>Mformfunc(idim)	 
     iOform=Mi%Oform; iLform=Mi%Lform; fki=formk(-vk,iLform)
     Oi=>Oformfunc(iOform); npi=Oi%npair
     do ip=1,npi; fOi=Oi%f(ip); iorb1=Oi%pair(1,ip); iorb2=Oi%pair(2,ip)
     if(natom>1)call modifyindex(iorb1,Mi%atom,iorb2,Mi%nbatom)
	    i=iorb1+(iorb2-1)*norb
		fdag(idim,i)=fdag(idim,i)+fOi*fki
	 end do; end do

	 f=conjg(transpose(fdag))  

     call gethk(norb,vk+q,Qc,hkq,model,MF,.false.); call ZHEIGEN(norb,hkq,evkq)
     call gethk(norb,vk,Qc,hk,model,MF,.false.); call ZHEIGEN(norb,hk,evk)
	 call gethk(norb,-vk,Qc,hkk,model,MF,.false.); call ZHEIGEN(norb,hkk,evkk)   !for time-invariant system, h(-k)=h(k)^*

     X =0; Y=0
	 do ib=1,norb; do jb=1,norb
	    
		fpp = lindhard( evkq(ib), -evkk(jb), cutoff )
		fph = lindhard( evkq(ib), evk(jb), cutoff )
             
        do jorb1=1,norb; do jorb2=1,norb; j=jorb1+(jorb2-1)*norb
        do iorb1=1,norb; do iorb2=1,norb; i=iorb1+(iorb2-1)*norb
           X(i,j) =  X(i,j) + hkq(iorb1,ib)*conjg(hkq(jorb1,ib))*hkk(iorb2,jb)*conjg(hkk(jorb2,jb))*fpp
	       Y(i,j) =  Y(i,j) + hkq(iorb1,ib)*conjg(hkq(jorb1,ib))* hk(jorb2,jb)*conjg( hk(iorb2,jb))*fph
	    end do; end do; end do; end do
	 end do; end do

	 X = X*sk;	 Xpp = Xpp + matmul( fdag, matmul( X, f ) )	  
	 Y = Y*sk;   Xph = Xph + matmul( fdag, matmul( Y, f ) )
  
  end do

  return
end subroutine bareSuscept

function Lindhard(ei,ej,cutoff)
  implicit none
  real*8 Lindhard,ei,ej,cutoff
  real*8 a,b,fa,fb,temp,pi
  complex*16 z
  character*20 scheme

  scheme='temperature'
  a=ei; b=ej; if(abs(a-b)<1.e-8)b=a+1.e-8

  select case (scheme)
    case ('temperature')
      fa=1/(1+exp(-abs(a)/cutoff)); if(a>0)fa=1-fa
      fb=1/(1+exp(-abs(b)/cutoff)); if(b>0)fb=1-fb
      Lindhard=(fa-fb)/(b-a)
	case ('energy')
	  Lindhard=0
	  if(abs(ei)>cutoff.or.abs(ej)>cutoff)return
	  temp=1.e-6
      fa=1/(1+exp(-abs(a)/temp)); if(a>0)fa=1-fa
      fb=1/(1+exp(-abs(b)/temp)); if(b>0)fb=1-fb
      Lindhard=(fa-fb)/(b-a) 
	case ('frequency')
	  z=cmplx(0,cutoff); pi=3.1415926
	  Lindhard=aimag( log( (z-a)/(z-b) ) ) /(a-b) /pi
	case default
	  stop 'error: scheme undefined @ Lindhard'
  end select

  return
end function Lindhard
