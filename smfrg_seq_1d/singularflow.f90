
subroutine singularflow(La,dLa,icount)
  use workplace
  implicit none
  
  real*8 La,dLa
  integer icount
  
  integer iq
  character*3 fileid,digitize
  real*8 q(2)

  complex*16, dimension (ndim,ndim,nq) :: dP,dC,dD,Peph,Ceph
  complex*16, dimension (ndim,ndim) :: Xpp,Xph

  integer ig,iqg,iqmax
  real*8 pv(3),pvg(3)
  logical reducible
  integer idim,jdim,iidim,jjdim
  integer isign,jsign,gsign
  integer sgn  
  logical binary; data binary/.true./

  if(useX0)then
    fileid=digitize(icount)
    !if(binary)open(10,file='./kuroki_KFS/mu11.30a0.8b1.0/'//'X0.'//fileid,form='binary')
    !if(binary)open(10,file='./kfs1lav1/'//'X0.'//fileid,form='binary')
    !if(binary)open(10,file='./lhq/lhqt0.3/'//'X0.'//fileid,form='binary')
    !if(binary)open(10,file='./111/lifeas_new/'//'X0.'//fileid,form='binary')
    if(binary)open(10,file='X0.'//fileid,form='UNFORMATTED')
    if(.not.binary)open(10,file='X0.'//fileid)
    !if(icount==321)then; useX0=.false.; appendX0=.false.; end if
  end if

  dP=0; dC=0; dD=0; Peph=0; Ceph=0
  if(model%ephcoupled)call ElectronPhonon(La,Peph,Ceph)

  iqmax=nq-model%nest-1
  do iq=1,nq; q=qv(1:2,iq)
  
     if(iq<=iqmax.and.usegroup.and.reducible(q,group))cycle
     
	 !-------------------------------------
	 if(useX0)then
       if(appendX0)then
	     if(binary)read(10)Xpp,Xph
		 if(.not.binary)read(10,*)Xpp,Xph
		 if(BCSflow.and.sum(abs(q))>1.e-6)cycle
	   else
	     call bubble(La,q,Xpp,Xph)
	     Xpp=Xpp*dLa/twopi; Xph=Xph*dLa/twopi
	     if(binary)write(10)Xpp,Xph
		 if(.not.binary)write(10,*)Xpp,Xph
		 cycle
	   end if
	 else
       if(BCSflow.and.sum(abs(q))>1.e-6)cycle
	   call bubble(La,q,Xpp,Xph)
	   Xpp=Xpp*dLa/twopi; Xph=Xph*dLa/twopi       
	 end if
	 !------------------------------------------
	 
	 dP(:,:,iq)=matmul( P(:,:,iq)+Peph(:,:,iq), matmul(Xpp,P(:,:,iq)+Peph(:,:,iq)) )
	 
	 if(.not.BCSflow)then
	   dC(:,:,iq)=matmul( C(:,:,iq)+Ceph(:,:,iq), matmul(Xph,C(:,:,iq)+Ceph(:,:,iq)) )
	   dD(:,:,iq)=matmul( C(:,:,iq)+Ceph(:,:,iq)-D(:,:,iq), matmul(Xph,D(:,:,iq)) ) &
	                       + matmul( D(:,:,iq), matmul(Xph,C(:,:,iq)+Ceph(:,:,iq)-D(:,:,iq)) )
	 end if
	 
	 if((.not.usegroup).or.iq>iqmax )cycle

	 call symmetrizePCD_q(iq,dP,dC,dD)
  end do

  if(useX0)then
    close(10)
	if(.not.AppendX0)then
	  print*,'X0.'//fileid,' for w,dw=',La,dLa; return
	end if
  end if
    
  if(.not.BCSflow)then
    if(useprojectors)then; call contact_clever(dP,dC,dD); else; call contact(dP,dC,dD); end if
  else
    dC=0; dD=0
  end if

  P=P+dP; C=C+dC; D=D+dD
  return
end subroutine singularflow

function digitize(icount)
  integer icount
  character*3 digitize

  character*1 string(3),stringtable(0:9)

  stringtable=(/'0','1','2','3','4','5','6','7','8','9'/)
  string(1)=stringtable( icount/100 )
  string(2)=stringtable( mod(icount,100)/10 )
  string(3)=stringtable( mod(icount,10) )
  digitize = (string(1)//string(2)//string(3) )

  return
end function digitize
