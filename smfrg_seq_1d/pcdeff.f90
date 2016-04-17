subroutine symmetrize(nkf,Peff,model)
  use standard_derived_types
  implicit none

  integer nkf
  complex*16 Peff(nkf,nkf),Pwork(nkf,nkf),doswork(nkf)
  type (modelblock) :: model
 
  integer ig,ng,band(nkf)
  character*3 group
  real*8 k(3),kg(3),kf(2,nkf)
  integer ik,jk,indexg(nkf)

  Pwork=Peff; doswork=model%pdos; kf=model%kf; band=model%band
  ng=model%ng; group=model%group

  k(3)=0
  do ig=2,ng
	 do ik=1,nkf
	    k(1:2)=kf(:,ik); call groupaction(ig,k,kg,group)
		indexg(ik)=findindex(band(ik),kg(1:2))
	 end do     
	 do jk=1,nkf; model%pdos(jk)=model%pdos(jk)+doswork(indexg(jk))
	 do ik=1,nkf; Peff(ik,jk)=Peff(ik,jk)+pwork(indexg(ik),indexg(jk))
	 end do; end do
  end do

  Peff=Peff/ng; model%pdos=model%pdos/ng

  return
  contains
    function findindex(ib,k)
	  integer findindex
	  integer ib,i,icount
	  real*8 k(2)
	  findindex=0; icount=0
	  do i=1,nkf
	     if(band(i)==ib.and.sum(abs(kf(1:2,i)-k))<1.e-3)then
		   findindex=i; icount=icount+1
		 end if
	  end do
	  if(findindex==0)stop 'k not identifiable @ findindex @ symmetrize'
	  if(icount>1)stop 'two many indices found @ findindex @ symmetrize'
      return
	end function findindex
end subroutine symmetrize 


subroutine getPeff(nkf,kf,band,Peff)
  use workplace
  implicit none

  integer nkf
  real*8 kf(2,nkf)
  integer band(nkf)
  complex*16, dimension(nkf,nkf) :: Peff

  integer iq,iP,ik,jk,iD,iC,idim,iatom,iOform,iorb,iiorb,leg1,leg2,jdim,jatom,jOform,jorb,jjorb,leg3,leg4
  real*8 G(2),ekk(norb),ep(norb)
  complex*16 hk(norb,norb),hp(norb,norb),uk(norb),up(norb)
  complex*16,dimension(ndim,ndim,nq) :: dP,dC,dD
  complex*16,dimension(ndim,ndim) :: CtoP,DtoP
  
  do iq=1,nq
    dP(:,:,iq)=P(:,:,iq)
	dC(:,:,iq)=C(:,:,iq)-Cloc
	dD(:,:,iq)=D(:,:,iq)-Dloc
  end do

  !contact terms to be eliminated to avoid overcounting
  do jdim=1,ndim; do idim=1,ndim
     CtoP(idim,jdim)=sum(dC(idim,jdim,:)*qv(3,:))
	 DtoP(idim,jdim)=sum(dD(idim,jdim,:)*qv(3,:))
  end do; end do

  iP=getindex( (/0d0,0d0/), G )
  
  Peff=0
  
  do ik=1,nkf; call gethk(norb,kf(:,ik),0d0,hk,model,MF,.false.); call ZHEIGEN(norb,hk,ekk); uk=hk(:,band(ik))
  do jk=1,ik; call gethk(norb,kf(:,jk),0d0,hp,model,MF,.false.); call ZHEIGEN(norb,hp,ep); up=hp(:,band(jk))

    iD=getindex( kf(:,jk)-kf(:,ik), G );  iC=getindex( kf(:,jk)+kf(:,ik), G )

    !only holds for on-site lattice forms and free orbital forms
	do idim=1,ndim
	  iatom=Mformfunc(idim)%atom
	  iOform=Mformfunc(idim)%Oform
	  if(Oformfunc(iOform)%npair/=1)stop 'Oform error @ getPeff' 
	  iorb=Oformfunc(iOform)%pair(1,1)
	  iiorb=Oformfunc(iOform)%pair(2,1)
	  if(natom>1)call modifyindex(iorb,iatom,iiorb,iatom)
      leg1=iorb
	  leg2=iiorb
	  do jdim=1,ndim
	    jatom=Mformfunc(jdim)%atom
	    jOform=Mformfunc(jdim)%Oform
	    if(Oformfunc(jOform)%npair/=1)stop 'Oform error @ getPeff' 
	    jorb=Oformfunc(jOform)%pair(1,1)
	    jjorb=Oformfunc(jOform)%pair(2,1)
	    if(natom>1)call modifyindex(jorb,jatom,jjorb,jatom)
        leg3=jorb
	    leg4=jjorb 

	    !Pairing channel: legs = ( k^,-k^,p,-p )
		Peff(ik,jk) = Peff(ik,jk) + dP(idim,jdim,iP) * conjg(uk(leg1)) * uk(leg2) * up(leg3) * conjg(up(leg4))
		                
		!Crossing channel: legs = ( k^,-p,p,-k^ )
		Peff(ik,jk) = Peff(ik,jk) + dC(idim,jdim,iC) * conjg(uk(leg1)) *conjg(up(leg2)) * up(leg3) * uk(leg4)
       
	    !Direct channel:  legs = ( k^,p,-p,-k^ )
		Peff(ik,jk) = Peff(ik,jk) + dD(idim,jdim,iD) * conjg(uk(leg1)) * up(leg2) * conjg(up(leg3)) * uk(leg4)

		if(iatom==jatom)then
		  !subtract overcounted contact term from Crossing channel: legs = ( k^,-p,p,-k^ )
		  !Peff(ik,jk) = Peff(ik,jk) - CtoP(idim,jdim) * conjg(uk(leg1)) *conjg(up(leg2)) * up(leg3) * uk(leg4)
       
	      !subtract overcounted contact term from Direct channel:  legs = ( k^,p,-p,-k^ )
		  !Peff(ik,jk) = Peff(ik,jk) - DtoP(idim,jdim) * conjg(uk(leg1)) * up(leg2) * conjg(up(leg3)) * uk(leg4)
		end if  
     end do; end do
	 Peff(jk,ik)=conjg(Peff(ik,jk))
  end do; end do

  return
  contains
    function getindex(q,G)
      integer getindex
	  real*8 q(2),p(2),g(2)
	  real*8 error,dq
	  integer iq
	  call rbzkv(q,model%a,model%b,p); G=q-p
	  getindex=1; error=sum(abs(p-qv(1:2,1)))
	  do iq=2,nq; dq=sum(abs(p-qv(1:2,iq))); if(dq<error)then; error=dq; getindex=iq; end if; end do
	  return
	end function getindex
end subroutine

subroutine getPeff_k(nkf,kf,band,Peff,modecoupled)
  use workplace
  implicit none

  integer nkf
  integer band(nkf)  !native band indices for each selected momenta on the fermi surface
  real*8 kf(2,nkf)
  complex*16, dimension (nkf,nkf) :: Peff
  logical modecoupled

  integer ik,jk,iq
  integer iqP,iqC,iqD
  complex*16, dimension (norb,norb) :: hk,hp
  complex*16, dimension (norb) :: uk,up,Tuk,Tup
  real*8 ework(norb)

  real*8 fOi,fOj,Q
  complex*16 fki,fkj,fatom,adag,a,formk
  complex*16, dimension (ndim) :: fPdag,fP,fCdag,fC,fDdag,fD

  integer idim,iOform,iLform,npi,ip,jdim
  integer iorb1,iorb2,leg1,leg2

  type (Mformconfig), pointer :: Mi
  type (Oformconfig), pointer :: Oi

  complex*16, dimension (ndim,ndim,nq) :: dP,dC,dD

  real*8 dr(2),GC(2),GD(2)

  !prepare matrices for states at k and p, assuming |-k> = |k>^* by time-reversal symmetry
 
  if(modecoupled)then
    dP=0; dC=C; dD=D; call contact(dP,dC,dD); dP=P-dP
  else
    dP=P
  end if
    
  iqP=getindex( (/0.d0,0.d0/), GD )

  Peff=0
  do jk=1,nkf; call gethk(norb,kf(:,jk),Q,hp,model,MF,.false.); call ZHEIGEN(norb,hp,ework); up=hp(:,band(jk)); Tup=conjg(up)
  do ik=1,jk;  call gethk(norb,kf(:,ik),Q,hk,model,MF,.false.); call ZHEIGEN(norb,hk,ework); uk=hk(:,band(ik)); Tuk=conjg(uk)

     iqD=getindex( kf(:,ik)-kf(:,jk), GD );  iqC=getindex( kf(:,jk)+kf(:,ik), GC )

     fPdag=0; fP=0; fCdag=0; fC=0; fDdag=0; fD=0

     do idim=1,ndim; Mi=>Mformfunc(idim)	 
     iOform=Mi%Oform; iLform=Mi%Lform; fki=formk(-kf(:,ik),iLform); fkj=formk(-kf(:,jk),iLform)
     Oi=>Oformfunc(iOform); npi=Oi%npair
     do ip=1,npi; fOi=Oi%f(ip); iorb1=Oi%pair(1,ip); iorb2=Oi%pair(2,ip)
     if(natom>1)call modifyindex(iorb1,Mi%atom,iorb2,Mi%nbatom)

        leg1=iorb1; leg2=iorb2
	   
	    !Pairing channel: legs = ( k^,-k^,p,-p )

        fPdag(idim) = fPdag(idim) + conjg(fki*uk(leg1)*Tuk(leg2))*fOi
		fP(idim) = fP(idim) + fkj*up(leg1)*Tup(leg2)*fOi
				                
		if(.not.modecoupled)cycle

		!Crossing channel: legs = ( k^,-p,p,-k^ )
		adag = fkj*conjg(uk(leg1)) * Tup(leg2) *fOi
		a = conjg(fki) * up(leg1) * conjg(Tuk(leg2)) * fOi
		if(.not.model%propergauge)then
		  fatom = exp(one*pi*sum(model%ra(1:2,Mi%atom)*GC) )
		  adag = adag * conjg(fatom); a = a*fatom
		end if
		fCdag(idim) = fCdag(idim) + adag; fC(idim) = fC(idim) + a
       
	    !Direct channel:  legs = ( k^,p,-p,-k^ )
	    adag = conjg( fkj * uk(leg1) ) * up(leg2) * fOi
		a = conjg(fki) * Tup(leg1) * conjg( Tuk(leg2) ) * fOi
		if(.not.model%propergauge)then
		  fatom = exp(one*pi*sum(model%ra(1:2,Mi%atom)*GD) )
		  adag = adag * conjg(fatom); a = a*fatom
		end if
		fddag(idim) = fDdag(idim) + adag; fD(idim) = fD(idim) + a
     end do; end do

     Peff(ik,jk) = Peff(ik,jk) + sum( fPdag * matmul( dP(:,:,iqP),fP) )
	 Peff(ik,jk) = Peff(ik,jk) + sum( fCdag * matmul( dC(:,:,iqC),fC) ) + sum( fDdag * matmul(dD(:,:,iqD),fD) ) 
	 Peff(jk,ik)=conjg(Peff(ik,jk))
	 if(ik==jk)Peff(ik,jk)=real(Peff(ik,jk))
  end do; end do
   
  return
  contains
    function getindex(q,G)
      integer getindex
	  real*8 q(2),p(2),g(2)
	  real*8 error,dq
	  integer iq
	  call rbzkv(q,model%a,model%b,p); G=q-p
	  getindex=1; error=sum(abs(p-qv(1:2,1)))
	  do iq=2,nq; dq=sum(abs(p-qv(1:2,iq))); if(dq<error)then; error=dq; getindex=iq; end if; end do
	  return
	end function getindex
end subroutine getPeff_k


subroutine getDeff_k(nkf,kf,band,Deff,modecoupled)
  use workplace
  implicit none

  integer nkf
  integer band(nkf)  !native band indices for each selected momenta on the fermi surface
  real*8 kf(2,nkf)
  complex*16, dimension (nkf,nkf) :: Deff
  logical modecoupled

  integer ik,jk
  integer iqP,iqC,iqD
  complex*16, dimension (norb,norb) :: hk,hp
  complex*16, dimension (norb) :: uk,up
  real*8 ework(norb)

  real*8 fOi,fOj,Q
  complex*16 fki,fkj,fatom,adag,a,formk

  integer idim,iOform,iLform,npi,ip
  integer iorb1,iorb2,leg1,leg2

  type (Mformconfig), pointer :: Mi
  type (Oformconfig), pointer :: Oi

  complex*16, dimension (ndim,ndim,nq) :: dP,dC,dD
  complex*16, dimension (ndim) :: fDdag,fD,fCdag,fC,fPdag,fP
  real*8, dimension (2) :: dr,GC,GP

  !prepare matrices for states at k and p, assuming |-k> = |k>^* by time-reversal symmetry
  if(modecoupled)then   
    !subtracting contact terms
    dD=0; dP=P; dC=C
	call contact(dP,dC,dD)
	dD=D-dD
  else
    dD=C-2*D
  end if

  iqD=getindex( (/0.d0,0.d0/),GP )

  Deff=0
  do jk=1,nkf; call gethk(norb,kf(:,jk),Q,hp,model,MF,.false.); call ZHEIGEN(norb,hp,ework); up=hp(:,band(jk))
  do ik=1,jk; call gethk(norb,kf(:,ik),Q,hk,model,MF,.false.); call ZHEIGEN(norb,hk,ework); uk=hk(:,band(ik))

	 iqP=getindex( kf(:,ik)+kf(:,jk),GP);  iqC=getindex( kf(:,ik)-kf(:,jk),GC )

     fDdag=0; fD=0; fPdag=0; fP=0; fCdag=0; fC=0

     do idim=1,ndim; Mi=>Mformfunc(idim)	 
     iOform=Mi%Oform; iLform=Mi%Lform; fki=formk(-kf(:,ik),iLform); fkj=formk(-kf(:,ik),iLform)
     Oi=>Oformfunc(iOform); npi=Oi%npair
     do ip=1,npi; fOi=Oi%f(ip); iorb1=Oi%pair(1,ip); iorb2=Oi%pair(2,ip)
     if(natom>1)call modifyindex(iorb1,Mi%atom,iorb2,Mi%nbatom)

        leg1=iorb1; leg2=iorb2
	   	   
        fDdag(idim) = fDdag(idim) + conjg(fki * uk(leg1) ) * uk(leg2) * fOi
		fD(idim) = fD(idim) + fkj * up(leg1) * conjg(up(leg2)) * fOi

        if(.not.modecoupled)cycle

		adag = conjg( fkj * uk(leg1) ) * up(leg2) * fOi
		a = fkj * uk(leg1) * conjg( up(leg2) ) * fOi
		if(.not.model%propergauge)then
		  fatom = exp( one*pi*sum(model%ra(1:2,Mi%atom)*GC) )
		  adag = adag * conjg(fatom); a = a * fatom
		end if
		fCdag(idim) = fCdag(idim) + adag; fC(idim) = fC(idim) + a

        adag = fkj * conjg( uk(leg1) * up(leg2) ) * fOi
		a = conjg(fki) *  up(leg1) * uk(leg2) * fOi - 2*conjg(fkj) * uk(leg1) * up(leg2) * fOi
		if(.not.model%propergauge)then
		  fatom = exp( one*pi*sum(model%ra(1:2,Mi%atom)*GC) )
		  adag = adag * conjg(fatom); a = a * fatom
		end if
		fPdag(idim) = fPdag(idim) + adag; fP(idim) = fP(idim) + a
     end do; end do
     
	 if(modecoupled)then
       Deff(ik,jk) = Deff(ik,jk) + sum( fDdag * matmul( C(:,:,iqD) - 2*dD(:,:,iqD), fD) )  &
	                             + sum( fCdag * matmul( dD(:,:,iqC) - 2*C(:,:,iqC), fC) ) &
								 + sum( fPdag * matmul( P(:,:,iqP), fP) )
	 else
	   Deff(ik,jk) = Deff(ik,jk) + sum( fDdag * matmul( dD(:,:,iqD), fD) )
	 end if
		
	 Deff(jk,ik)=conjg(Deff(ik,jk))
	 if(ik==jk)Deff(ik,jk)=real(Deff(ik,jk))
  end do; end do
   
  return
  contains
    function getindex(q,g)
      integer getindex
	  real*8 q(2),p(2),g(2)
	  real*8 error,dq
	  integer iq
	  call rbzkv(q,model%a,model%b,p); g=q-p
	  getindex=1; error=sum(abs(p-qv(1:2,1)))
	  do iq=2,nq; dq=sum(abs(p-qv(1:2,iq))); if(dq<error)then; error=dq; getindex=iq; end if; end do
	  return
	end function getindex
end subroutine getDeff_k


subroutine getCeff_k(nkf,kf,band,Ceff,modecoupled)
  use workplace
  implicit none

  integer nkf
  integer band(nkf)  !native band indices for each selected momenta on the fermi surface
  real*8 kf(2,nkf)
  complex*16, dimension (nkf,nkf) :: Ceff
  logical modecoupled

  integer ik,jk
  integer iqP,iqC,iqD
  complex*16, dimension (norb,norb) :: hk,hp
  complex*16, dimension (norb) :: uk,up
  real*8 ework(norb)

  real*8 fOi,fOj,Q
  complex*16 fki,fkj,adag,a,fatom,formk

  integer idim,iOform,iLform,npi,ip
  integer iorb1,iorb2,leg1,leg2

  type (Mformconfig), pointer :: Mi
  type (Oformconfig), pointer :: Oi

  complex*16, dimension (ndim,ndim,nq) :: dP,dC,dD
  complex*16, dimension (ndim) :: fCdag,fC,fPdag,fP,fDdag,fD

  real*8, dimension (2) :: dr,GP,GD

  if(modecoupled)then   
    !subtracting contact terms
    dP=P; dC=0; dD=D; call contact(dP,dC,dD); dC=C-dC
  else
    dC=C
  end if

  iqC=getindex( (/0.d0,0.d0/),GD )

  Ceff=0
  do jk=1,nkf; call gethk(norb,kf(:,jk),Q,hp,model,MF,.false.); call ZHEIGEN(norb,hp,ework); up=hp(:,band(jk))
  do ik=1,jk; call gethk(norb,kf(:,ik),Q,hk,model,MF,.false.); call ZHEIGEN(norb,hk,ework); uk=hk(:,band(ik))

     iqP=getindex( kf(:,jk)+kf(:,ik),GP );  iqD=getindex( kf(:,ik)-kf(:,jk),GD )

     fCdag=0; fC=0; fPdag=0; fP=0; fDdag=0; fD=0

     do idim=1,ndim; Mi=>Mformfunc(idim)	 
     iOform=Mi%Oform; iLform=Mi%Lform; fki=formk(-kf(:,ik),iLform); fkj=formk(-kf(:,jk),iLform)
     Oi=>Oformfunc(iOform); npi=Oi%npair
     do ip=1,npi; fOi=Oi%f(ip); iorb1=Oi%pair(1,ip); iorb2=Oi%pair(2,ip)
     if(natom>1)call modifyindex(iorb1,Mi%atom,iorb2,Mi%nbatom)

        leg1=iorb1; leg2=iorb2
	    
		fCdag(idim) = fCdag(idim) + conjg( fki*uk(leg1) ) * uk(leg2) *fOi
		fC(idim) = fC(idim) + fkj * up(leg1) * conjg(up(leg2)) * fOi

        if(.not.modecoupled)cycle
	
	    adag = fkj * conjg( uk(leg1) * up(leg2) ) * fOi
		a = conjg(fki) * up(leg1) * uk(leg2) * fOi
		if(.not.model%propergauge)then
		  fatom = exp(one*pi*sum(model%ra(1:2,Mi%atom)*GP) )
		  adag = adag * conjg(fatom); a = a * fatom
		end if
        fPdag(idim) = fPdag(idim) + adag; fP(idim) = fP(idim) + a


        adag = conjg( fkj * uk(leg1) ) * up(leg2) * fOi
		a = fkj * uk(leg1) * conjg(up(leg2)) * fOi	
		if(.not.model%propergauge)then
		  fatom = exp(one*pi*sum(model%ra(1:2,Mi%atom)*GD) )
		  adag = adag * conjg(fatom); a = a * fatom
		end if
		fDdag(idim) = fDdag(idim) + adag; fD(idim) = fD(idim) + a	
     end do; end do

	 Ceff(ik,jk) = Ceff(ik,jk) + sum( fCdag * matmul( dC(:,:,iqC), fC) )
	 if(modecoupled) Ceff(ik,jk) = Ceff(ik,jk) + sum( fPdag * matmul( P(:,:,iqP), fP) ) + sum( fDdag * matmul( D(:,:,iqD), fD) )
	 Ceff(jk,ik)=conjg(Ceff(ik,jk))
	 if(ik==jk)Ceff(ik,jk)=real(Ceff(ik,jk))
  end do; end do
  
  return
  contains
    function getindex(q,G)
      integer getindex
	  real*8 q(2),p(2),G(2)
	  real*8 error,dq
	  integer iq
	  call rbzkv(q,model%a,model%b,p); g=q-p
	  getindex=1; error=sum(abs(p-qv(1:2,1)))
	  do iq=2,nq; dq=sum(abs(p-qv(1:2,iq))); if(dq<error)then; error=dq; getindex=iq; end if; end do
	  return
	end function getindex
end subroutine getCeff_k
