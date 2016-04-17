subroutine contact(dP,dC,dD)
  use workplace
  implicit none

  complex*16, dimension (ndim,ndim,nq) :: dP,dC,dD

  complex*16, dimension (ndim,ndim,nR) :: dPR,dCR,dDR
  integer iq,iR,iqmax
  integer idim,jdim,iatom,jatom

  real*8, dimension (3) :: Sc2c,S,ri,rj
	real*8, pointer, dimension (:,:) :: ra

  complex*16 factor,fq(nq)
	logical reducible

  ra=>model%ra

  !FT dP,dC and dD to dPR,dCR,dDR
  dPR=0; dCR=0; dDR=0
  do iR=1,nR; Sc2c=R(1,iR)*model%a+R(2,iR)*model%b
  if(QuickSearch.and.sum(Sc2c(1:2)**2)>1.e-5)cycle
  do iatom=1,natom; do jatom=1,natom
    S=Sc2c+ra(:,jatom)-ra(:,iatom)
    if( usegroup.and.reducible(S(1:2),group) )cycle

    do iq=1,nq; fq(iq) = exp( -one*pi*sum( qv(1:2,iq)*S(1:2) )  ) * qv(3,iq); end do

    do jdim=1,ndim; if(Mformfunc(jdim)%atom/=jatom)cycle
    do idim=1,ndim; if(Mformfunc(idim)%atom/=iatom)cycle
      dPR(idim,jdim,iR) =  sum( fq*dP(idim,jdim,:) ) 
      dCR(idim,jdim,iR) =  sum( fq*dC(idim,jdim,:) )
      dDR(idim,jdim,iR) =  sum( fq*dD(idim,jdim,:) )
    end do; end do
    if(usegroup)call symmetrizePCD_R(iR,iatom,jatom,dPR,dCR,dDR)
  end do; end do; end do

  !get contact contributions
  call contactkernel(dPR,dCR,dDR)

  !FT dPR,dCR and dDR and add to dP,dC,dD
  iqmax=nq-model%nest-1
  do iq=1,nq
    if(iq<=iqmax.and.usegroup.and.reducible(qv(1:2,iq),group) )cycle
    do iR=1,nR; Sc2c=R(1,iR)*model%a+R(2,iR)*model%b
      if(QuickSearch.and.sum(Sc2c(1:2)**2)>1.e-5)cycle
      do jdim=1,ndim
        rj=ra(:,Mformfunc(jdim)%atom )
        do idim=1,ndim
          ri=ra(:,Mformfunc(idim)%atom )
          S=Sc2c+rj-ri
          factor=exp( one*pi*sum( qv(1:2,iq)*S(1:2) )  )  
          dP(idim,jdim,iq)=dP(idim,jdim,iq) +  factor*dPR(idim,jdim,iR)
          dC(idim,jdim,iq)=dC(idim,jdim,iq) +  factor*dCR(idim,jdim,iR)
          dD(idim,jdim,iq)=dD(idim,jdim,iq) +  factor*dDR(idim,jdim,iR)
        end do
      end do
    end do
    if((.not.usegroup).or.iq>iqmax)cycle
    call symmetrizePCD_q(iq,dP,dC,dD)
  end do

  return
end subroutine contact

subroutine contactkernel(dPR,dCR,dDR)
  use workplace
  implicit none

  complex*16, dimension (ndim,ndim,nR) :: dPR,dCR,dDR

  !purpose: get mutual overlaps in dP,dC,dD
  !notice:  on entry, dPR,dCR,dDR are fourier transform of dP,dC,dD from singular flow
  !         on exit,  dPR,dCR,dDR are contact contributions, e.g., dPR=Proj(dCR,dDR), ...
  
  complex*16, dimension (ndim,ndim,nR) :: dPc,dCc,dDc
  integer iR,ia,ja
  integer idim,iOform,iLform,npi,nbi,ip,i
  integer jdim,jOform,jLform,npj,nbj,jp,j
  integer iorb1,iorb2,jorb1,jorb2,orb(4)
  integer iatom,iiatom,jatom,jjatom
  integer findatom

  type (Mformconfig), pointer :: Mi,Mj
  type (Lformconfig), pointer :: Li,Lj
  type (Oformconfig), pointer :: Oi,Oj

  real*8, dimension (3) :: R1,R2,R3,R4,ri,rj,Sc2c,S,SS,Rc2c
  real*8 fOi,fOj,fLi,fLj
	real*8, pointer, dimension (:,:) :: ra

  complex*16 add,factor
	logical near,Pcoupled,reducible

  Pcoupled=.true.

  ra=>model%ra

  dPc=dPR; dCc=dCR; dDc=dDR
  dPR=0; dCR=0; dDR=0

  do iR=1,nR; S=R(1,iR)*model%a+R(2,iR)*model%b; if(QuickSearch.and.sum(S*S)>1.e-5)cycle
  do ja=1,natom; do ia=1,natom; SS=S+ra(:,ja)-ra(:,ia)
    if(usegroup.and.reducible(SS(1:2),group))cycle
    do jdim=1,ndim
      if(Mformfunc(jdim)%atom/=ja)cycle
      do idim=1,ndim; if(Mformfunc(idim)%atom/=ia)cycle
        dPR(idim,jdim,iR)=0; dCR(idim,jdim,iR)=0; dDR(idim,jdim,iR)=0
      end do
    end do

    do jdim=1,ndim
      Mj=>Mformfunc(jdim)
      jatom=Mj%atom
      if(jatom/=ja)cycle
      jOform=Mj%Oform
      jLform=Mj%Lform
      if(QuickSearch.and.jLform/=1)cycle
      Oj=>Oformfunc(jOform)
      Lj=>Lformfunc(jLform)
      npj=Oj%npair
      do jp=1,npj
        fOj=Oj%f(jp)
        jorb1=Oj%pair(1,jp)
        jorb2=Oj%pair(2,jp)
        if(natom>1)call modifyindex(jorb1,Mj%atom,jorb2,Mj%nbatom)
        nbj=Lj%nr
        do j=1,nbj
          fLj=Lj%f(j)
          rj=Lj%r(:,j)
          do idim=1,ndim; Mi=>Mformfunc(idim); iatom=Mi%atom; if(iatom/=ia)cycle
            iOform=Mi%Oform; iLform=Mi%Lform; if(QuickSearch.and.iLform/=1)cycle
	   Oi=>Oformfunc(iOform); Li=>Lformfunc(iLform)
	   npi=Oi%npair; do ip=1,npi; fOi=Oi%f(ip)
	   iorb1=Oi%pair(1,ip); iorb2=Oi%pair(2,ip)
	   if(natom>1)call modifyindex(iorb1,Mi%atom,iorb2,Mi%nbatom)
	   nbi=Li%nr; do i=1,nbi; fLi=Li%f(i); ri=Li%r(:,i)

       factor=fLi*fOi*fLj*fOj

	   !get contributions to dP from dC and dD
	   if(Pcoupled)then
	     add=0; R1=ra(:,iatom); R2=R1+ri; R4=S+ra(:,jatom); R3=R4+rj
 	     if(near(R3-R1).and.near(R2-R4))then
	       Rc2c=S; orb=(/iorb1,jorb2,jorb1,iorb2/)
	       call addfrom(dCc,orb,R3-R1,Rc2c,R2-R4,add)
         end if
	     if(near(R4-R1).and.near(R2-R3))then
	       iiatom=iatom; jjatom=findatom(R3)
	       Rc2c=R3-ra(:,jjatom)-(R1-ra(:,iiatom))
		   orb=(/iorb1,jorb1,jorb2,iorb2/)
	       call addfrom(dDc,orb,R4-R1,Rc2c,R2-R3,add)	    
	     end if
		 dPR(idim,jdim,iR)=dPR(idim,jdim,iR)+factor*add
       end if

       !get contributions to dC from dP and dD
	   add=0;  R1=ra(:,iatom); R3=R1+ri; R4=S+ra(:,jatom); R2=R4+rj
	   if(Pcoupled.and.near(R2-R1).and.near(R3-R4))then
	     Rc2c=S; orb=(/iorb1,jorb2,jorb1,iorb2/)
	     call addfrom(dPc,orb,R2-R1,Rc2c,R3-R4,add)
	   end if
	   if(near(R4-R1).and.near(R2-R3))then
	     iiatom=iatom; jjatom=findatom(R3)
	     Rc2c= R3-ra(:,jjatom)-(R1-ra(:,iiatom)) 
         orb=(/iorb1,jorb1,iorb2,jorb2/)
	     call addfrom(dDc,orb,R4-R1,Rc2c,R2-R3,add)
	   end if
	   dCR(idim,jdim,iR)=dCR(idim,jdim,iR)+factor*add

	   !get contributions to dD from dP and dC
	   add=0; R1=ra(:,iatom); R4=R1+ri; R3=S+ra(:,jatom); R2=R3+rj
	   if(Pcoupled.and.near(R2-R1).and.near(R3-R4))then
	     iiatom=iatom; jjatom=findatom(R4)
	     Rc2c= R4-ra(:,jjatom)-(R1-ra(:,iiatom)) 
         orb=(/iorb1,jorb2,iorb2,jorb1/)
	     call addfrom(dPc,orb,R2-R1,Rc2c,R3-R4,add)
	   end if
	   if(near(R3-R1).and.near(R2-R4))then
	     iiatom=iatom; jjatom=findatom(R4)
	     Rc2c= R4-ra(:,jjatom)-(R1-ra(:,iiatom)) 
         orb=(/iorb1,jorb1,iorb2,jorb2/)
	     call addfrom(dCc,orb,R3-R1,Rc2c,R2-R4,add)
	   end if
	   dDR(idim,jdim,iR)=dDR(idim,jdim,iR)+factor*add

       end do; end do; end do
	   end do; end do; end do

	   if(usegroup)call symmetrizePCD_R(iR,ia,ja,dPR,dCR,dDR)
    end do; end do; end do

    return
end subroutine contactkernel

subroutine addfrom(dV,orb,ri,Sc2c,rj,add)
    use workplace
    implicit none

    integer orb(4)
    real*8, dimension (3) :: Sc2c,ri,rj
    complex*16 dV(ndim,ndim,nR),add

    real*8, dimension (3) :: rii,rjj

	integer idim,iOform,iLform,npi,nbi,ip,i
	integer jdim,jOform,jLform,npj,nbj,jp,j
    integer iorb1,iorb2,jorb1,jorb2

    type (Mformconfig), pointer :: Mi,Mj
	type (Lformconfig), pointer :: Li,Lj
	type (Oformconfig), pointer :: Oi,Oj

    integer iR,Rc2c(2)
    real*8 fOi,fLi,fOj,fLj,fi(ndim),fj(ndim)
	logical mismatch
	integer findR
	integer ni,nj,infoi(ndim),infoj(ndim)
  
    Rc2c(1)=nint( sum(Sc2c(1:2)*model%ka)/2 )
	Rc2c(2)=nint( sum(Sc2c(1:2)*model%kb)/2 )
	if(sum(abs( Sc2c-Rc2c(1)*model%a-Rc2c(2)*model%b ) )>1.e-5)stop 'Sc2c error @ addfrom'

    !do i=1,2; if(abs(Rc2c(i))>Lcontact)pause 'Rc2c error @ addto'; end do
    iR=findR(Rc2c); if(iR<=0)stop 'Rc2c error @ addfrom'

    fj=0; nj=0; infoj=0
	do jdim=1,ndim; Mj=>Mformfunc(jdim)
	   jOform=Mj%Oform; jLform=Mj%Lform; if(QuickSearch.and.jLform/=1)cycle
	   Oj=>Oformfunc(jOform); Lj=>Lformfunc(jLform)
	   npj=Oj%npair; do jp=1,npj; fOj=Oj%f(jp)
       jorb1=Oj%pair(1,jp); jorb2=Oj%pair(2,jp)
	   if(natom>1)call modifyindex(jorb1,Mj%atom,jorb2,Mj%nbatom)
	   if( mismatch(orb(3:4),(/jorb1,jorb2/)) )cycle
	   nbj=Lj%nr; do j=1,nbj; fLj=Lj%f(j); rjj=Lj%r(:,j)
       if(sum( abs(rj-rjj) )>1.e-5)cycle
	   fj(jdim) = fj(jdim) + fOj*fLj
       end do; end do 
	   if(abs(fj(jdim))<1.e-5)cycle
	   nj=nj+1; infoj(nj)=jdim
	end do

    fi=0; ni=0; infoi=0
    do idim=1,ndim; Mi=>Mformfunc(idim)
	   iOform=Mi%Oform; iLform=Mi%Lform; if(QuickSearch.and.iLform/=1)cycle
	   Oi=>Oformfunc(iOform); Li=>Lformfunc(iLform)
	   npi=Oi%npair; do ip=1,npi; fOi=Oi%f(ip)
       iorb1=Oi%pair(1,ip); iorb2=Oi%pair(2,ip)
	   if(natom>1)call modifyindex(iorb1,Mi%atom,iorb2,Mi%nbatom)
	   if( mismatch(orb(1:2),(/iorb1,iorb2/) ) )cycle
	   nbi=Li%nr; do i=1,nbi; fLi=Li%f(i); rii=Li%r(:,i)
	   if(sum( abs(ri-rii) )>1.e-5)cycle
       fi(idim) = fi(idim) + fOi*fLi
	   end do; end do
	   if(abs(fi(idim))<1.e-5)cycle
	   ni=ni+1; infoi(ni)=idim
	end do

	do j=1,nj; jdim=infoj(j); do i=1,ni; idim=infoi(i)		
       add = add + dV(idim,jdim,iR)*fi(idim)*fj(jdim)    
    end do; end do
 
    return
end subroutine addfrom	        

function mismatch(p1,p2)
  implicit none
  integer p1(2),p2(2)
  logical mismatch
  mismatch=( sum( abs(p1-p2) )>0 )
  return
end function mismatch

function near(r)
    use workplace, only : Lcontact
	logical near
	real*8 r(3),radius
	integer i
	radius=Lcontact/2.+1.e-5
	near=.true.
	do i=1,3; if(abs(r(i))>radius)near=.false.; end do
	return
end function near

































subroutine contact_old(dP,dC,dD)
    use workplace
    implicit none

    complex*16, dimension (ndim,ndim,nq) :: dP,dC,dD
    complex*16, dimension (ndim,ndim,nR) :: dPR,dCR,dDR,dPc,dCc,dDc

    integer iq,iR,iqmax
	integer idim,iOform,iLform,npi,nbi,ip,i
	integer jdim,jOform,jLform,npj,nbj,jp,j
    integer iorb1,iorb2,jorb1,jorb2,orb(4)
    integer iatom,iiatom,jatom,jjatom
	integer findatom

    type (Mformconfig), pointer :: Mi,Mj
	type (Lformconfig), pointer :: Li,Lj
	type (Oformconfig), pointer :: Oi,Oj

    real*8, dimension (3) :: R1,R2,R3,R4,ri,rj,Sc2c,S,Rc2c
    real*8 fOi,fOj,fLi,fLj
	real*8, pointer, dimension (:,:) :: ra

    complex*16 add,factor,fq(nq)
	logical near,Pcoupled,reducible

    Pcoupled=.true.

    ra=>model%ra

    dPR=0; dCR=0; dDR=0

    do iR=1,nR; Sc2c=R(1,iR)*model%a+R(2,iR)*model%b
	if(QuickSearch.and.sum(Sc2c(1:2)**2)>1.e-5)cycle
	do iatom=1,natom; do jatom=1,natom	 
	   S=Sc2c+ra(:,jatom)-ra(:,iatom)
	   if( usegroup.and.reducible(S(1:2),group) )cycle
	   	
	   do iq=1,nq; fq(iq) = exp( -one*pi*sum( qv(1:2,iq)*S(1:2) )  ) * qv(3,iq); end do

       do jdim=1,ndim; if(Mformfunc(jdim)%atom/=jatom)cycle
	   do idim=1,ndim; if(Mformfunc(idim)%atom/=iatom)cycle
	      dPR(idim,jdim,iR) =  sum( fq*dP(idim,jdim,:) ) 
	      dCR(idim,jdim,iR) =  sum( fq*dC(idim,jdim,:) )
	      dDR(idim,jdim,iR) =  sum( fq*dD(idim,jdim,:) )
       end do; end do
	   if(usegroup)call symmetrizePCD_R(iR,iatom,jatom,dPR,dCR,dDR)
	end do; end do; end do

    dPc=0; dCc=0; dDc=0
    do iR=1,nR; S=R(1,iR)*model%a+R(2,iR)*model%b; if(QuickSearch.and.sum(S*S)>1.e-5)cycle
       do jdim=1,ndim; Mj=>Mformfunc(jdim); jatom=Mj%atom
	   jOform=Mj%Oform; jLform=Mj%Lform; if(QuickSearch.and.jLform/=1)cycle
	   Oj=>Oformfunc(jOform); Lj=>Lformfunc(jLform)
	   npj=Oj%npair; do jp=1,npj; fOj=Oj%f(jp)
	   jorb1=Oj%pair(1,jp); jorb2=Oj%pair(2,jp)
	   if(natom>1)call modifyindex(jorb1,Mj%atom,jorb2,Mj%nbatom)
	   nbj=Lj%nr; do j=1,nbj; fLj=Lj%f(j); rj=Lj%r(:,j)

       do idim=1,ndim; Mi=>Mformfunc(idim); iatom=Mi%atom
	   iOform=Mi%Oform; iLform=Mi%Lform; if(QuickSearch.and.iLform/=1)cycle
	   Oi=>Oformfunc(iOform); Li=>Lformfunc(iLform)
	   npi=Oi%npair; do ip=1,npi; fOi=Oi%f(ip)
	   iorb1=Oi%pair(1,ip); iorb2=Oi%pair(2,ip)
	   if(natom>1)call modifyindex(iorb1,Mi%atom,iorb2,Mi%nbatom)
	   nbi=Li%nr; do i=1,nbi; fLi=Li%f(i); ri=Li%r(:,i)

       factor=fLi*fOi*fLj*fOj

	   !project dP to dC and dD
	   if(Pcoupled)then
	   add=factor*dPR(idim,jdim,iR)
	   R1=ra(:,iatom); R2=R1+ri; R4=S+ra(:,jatom); R3=R4+rj
 	   if(near(R3-R1).and.near(R2-R4))then
	     Rc2c=S; orb=(/iorb1,jorb2,jorb1,iorb2/)
	     call addto(dCc,orb,R3-R1,Rc2c,R2-R4,add)
       end if
	   if(near(R4-R1).and.near(R2-R3))then
	     iiatom=iatom; jjatom=findatom(R3)
	     Rc2c= R3-ra(:,jjatom)-(R1-ra(:,iiatom)) 
		 orb=(/iorb1,jorb1,jorb2,iorb2/)
	     call addto(dDc,orb,R4-R1,Rc2c,R2-R3,add)	    
	   end if; end if

       !project dC to dP and dD
	   add=factor*dCR(idim,jdim,iR)
       R1=ra(:,iatom); R3=R1+ri; R4=S+ra(:,jatom); R2=R4+rj
	   if(near(R2-R1).and.near(R3-R4))then
	     Rc2c=S; orb=(/iorb1,jorb2,jorb1,iorb2/)
	     if(Pcoupled)call addto(dPc,orb,R2-R1,Rc2c,R3-R4,add)
	   end if
	   if(near(R4-R1).and.near(R2-R3))then
	     iiatom=iatom; jjatom=findatom(R3)
	     Rc2c= R3-ra(:,jjatom)-(R1-ra(:,iiatom)) 
         orb=(/iorb1,jorb1,iorb2,jorb2/)
	     call addto(dDc,orb,R4-R1,Rc2c,R2-R3,add)
	   end if

	   !project dD to dP and dC
	   add=factor*dDR(idim,jdim,iR)
	   R1=ra(:,iatom); R4=R1+ri; R3=S+ra(:,jatom); R2=R3+rj
	   if(near(R2-R1).and.near(R3-R4))then
	     iiatom=iatom; jjatom=findatom(R4)
	     Rc2c= R4-ra(:,jjatom)-(R1-ra(:,iiatom)) 
         orb=(/iorb1,jorb2,iorb2,jorb1/)
	     if(Pcoupled)call addto(dPc,orb,R2-R1,Rc2c,R3-R4,add)
	   end if
	   if(near(R3-R1).and.near(R2-R4))then
	     iiatom=iatom; jjatom=findatom(R4)
	     Rc2c= R4-ra(:,jjatom)-(R1-ra(:,iiatom)) 
         orb=(/iorb1,jorb1,iorb2,jorb2/)
	     call addto(dCc,orb,R3-R1,Rc2c,R2-R4,add)
	   end if

       end do; end do; end do
	   end do; end do; end do
    end do

    !FT dPc,dCc and dDc and add to dP,dC,dD
	iqmax=nq-model%nest-1
    do iq=1,nq; if(iq<=iqmax.and.usegroup.and.reducible(qv(1:2,iq),group) )cycle
	   do iR=1,nR; Sc2c=R(1,iR)*model%a+R(2,iR)*model%b
	   if(QuickSearch.and.sum(Sc2c(1:2)**2)>1.e-5)cycle
       do jdim=1,ndim; rj=ra(:,Mformfunc(jdim)%atom )
	   do idim=1,ndim; ri=ra(:,Mformfunc(idim)%atom )
	      S=Sc2c+rj-ri
	  	  factor=exp( one*pi*sum( qv(1:2,iq)*S(1:2) )  )  
	      dP(idim,jdim,iq)=dP(idim,jdim,iq) +  factor*dPc(idim,jdim,iR)
	      dC(idim,jdim,iq)=dC(idim,jdim,iq) +  factor*dCc(idim,jdim,iR)
	      dD(idim,jdim,iq)=dD(idim,jdim,iq) +  factor*dDc(idim,jdim,iR)
       end do; end do; end do
       if((.not.usegroup).or.iq>iqmax)cycle
	   call symmetrizePCD_q(iq,dP,dC,dD)
    end do

    return
end subroutine contact_old

subroutine addto(dV,orb,ri,Sc2c,rj,add)
    use workplace
    implicit none

    integer orb(4)
    real*8, dimension (3) :: ri,rj,Sc2c
    complex*16 dV(ndim,ndim,nR),add

    
    real*8, dimension (3) :: rii,rjj

    integer Rc2c(2)
	integer idim,iOform,iLform,npi,nbi,ip,i
	integer jdim,jOform,jLform,npj,nbj,jp,j
    integer iorb1,iorb2,jorb1,jorb2

    type (Mformconfig), pointer :: Mi,Mj
	type (Lformconfig), pointer :: Li,Lj
	type (Oformconfig), pointer :: Oi,Oj

    integer iR
    real*8 fOi,fLi,fOj,fLj,fi(ndim),fj(ndim)
	logical mismatch
	integer findR
  
    Rc2c(1)=nint( sum(Sc2c(1:2)*model%ka)/2 )
	Rc2c(2)=nint( sum(Sc2c(1:2)*model%kb)/2 )
	if(sum(abs( Sc2c-Rc2c(1)*model%a-Rc2c(2)*model%b ) )>1.e-5)stop 'Sc2c error @ addfrom'

    !do i=1,2; if(abs(Rc2c(i))>Lcontact)pause 'Rc2c error @ addto'; end do
    iR=findR(Rc2c); if(iR<=0)stop 'Rc2c error @ addfrom'

    fj=0
	do jdim=1,ndim; Mj=>Mformfunc(jdim)
	   jOform=Mj%Oform; jLform=Mj%Lform; if(QuickSearch.and.jLform/=1)cycle
	   Oj=>Oformfunc(jOform); Lj=>Lformfunc(jLform)
	   npj=Oj%npair; do jp=1,npj; fOj=Oj%f(jp)
       jorb1=Oj%pair(1,jp); jorb2=Oj%pair(2,jp)
	   if(natom>1)call modifyindex(jorb1,Mj%atom,jorb2,Mj%nbatom)
	   if( mismatch(orb(3:4),(/jorb1,jorb2/)) )cycle
	   nbj=Lj%nr; do j=1,nbj; fLj=Lj%f(j); rjj=Lj%r(:,j)
       if(sum( abs(rj-rjj) )>1.e-5)cycle
	   fj(jdim) = fj(jdim) + fOj*fLj
    end do; end do; end do

    fi=0
    do idim=1,ndim; Mi=>Mformfunc(idim)
	   iOform=Mi%Oform; iLform=Mi%Lform; if(QuickSearch.and.iLform/=1)cycle
	   Oi=>Oformfunc(iOform); Li=>Lformfunc(iLform)
	   npi=Oi%npair; do ip=1,npi; fOi=Oi%f(ip)
       iorb1=Oi%pair(1,ip); iorb2=Oi%pair(2,ip)
	   if(natom>1)call modifyindex(iorb1,Mi%atom,iorb2,Mi%nbatom)
	   if( mismatch(orb(1:2),(/iorb1,iorb2/) ) )cycle
	   nbi=Li%nr; do i=1,nbi; fLi=Li%f(i); rii=Li%r(:,i)
	   if(sum( abs(ri-rii) )>1.e-5)cycle
       fi(idim) = fi(idim) + fOi*fLi
	end do; end do; end do

	do jdim=1,ndim; do idim=1,ndim		
       dV(idim,jdim,iR)=dV(idim,jdim,iR)+add*fi(idim)*fj(jdim)    
    end do; end do
 
    return
end subroutine addto

subroutine projection(dPR,dCR,dDR,dPc,dCc,dDc,Vsource,iR,idim,jdim)
    use workplace
    implicit none

    complex*16, dimension (ndim,ndim,nR) :: dPR,dCR,dDR,dPc,dCc,dDc

    character*1 Vsource
    integer iR,idim,jdim
	integer iOform,iLform,npi,nbi,ip,i
	integer jOform,jLform,npj,nbj,jp,j
    integer iorb1,iorb2,jorb1,jorb2,orb(4)
    integer iatom,iiatom,jatom,jjatom
	integer findatom

    type (Mformconfig), pointer :: Mi,Mj
	type (Lformconfig), pointer :: Li,Lj
	type (Oformconfig), pointer :: Oi,Oj

    real*8, dimension (3) :: R1,R2,R3,R4,ri,rj,Rc2c,S
    real*8 fOi,fOj,fLi,fLj
	real*8, pointer, dimension (:,:) :: ra

    complex*16 add,factor,factorc
	logical near,Pcoupled

    Pcoupled=.true.

    ra=>model%ra

    dPc=0; dCc=0; dDc=0
    S=R(1,iR)*model%a+R(2,iR)*model%b
       Mj=>Mformfunc(jdim); jatom=Mj%atom
	   jOform=Mj%Oform; jLform=Mj%Lform
	   Oj=>Oformfunc(jOform); Lj=>Lformfunc(jLform)
	   npj=Oj%npair; do jp=1,npj; fOj=Oj%f(jp)
	   jorb1=Oj%pair(1,jp); jorb2=Oj%pair(2,jp)
	   if(natom>1)call modifyindex(jorb1,Mj%atom,jorb2,Mj%nbatom)
	   nbj=Lj%nr; do j=1,nbj; fLj=Lj%f(j); rj=Lj%r(:,j)

       Mi=>Mformfunc(idim); iatom=Mi%atom
	   iOform=Mi%Oform; iLform=Mi%Lform
	   Oi=>Oformfunc(iOform); Li=>Lformfunc(iLform)
	   npi=Oi%npair; do ip=1,npi; fOi=Oi%f(ip)
	   iorb1=Oi%pair(1,ip); iorb2=Oi%pair(2,ip)
	   if(natom>1)call modifyindex(iorb1,Mi%atom,iorb2,Mi%nbatom)
	   nbi=Li%nr; do i=1,nbi; fLi=Li%f(i); ri=Li%r(:,i)

       factor=fLi*fOi*fLj*fOj

       if(Vsource=='P')then
	   !project dP to dC and dD
	   if(Pcoupled)then
	   add=factor*dPR(idim,jdim,iR)
	   R1=ra(:,iatom); R2=R1+ri; R4=S+ra(:,jatom); R3=R4+rj
 	   if(near(R3-R1).and.near(R2-R4))then
	     Rc2c=S; orb=(/iorb1,jorb2,jorb1,iorb2/)
	     call addto(dCc,orb,R3-R1,Rc2c,R2-R4,add)
       end if
	   if(near(R4-R1).and.near(R2-R3))then
	     iiatom=iatom; jjatom=findatom(R3)
	     Rc2c= R3-ra(:,jjatom)-(R1-ra(:,iiatom)) 
		 orb=(/iorb1,jorb1,jorb2,iorb2/)
	     call addto(dDc,orb,R4-R1,Rc2c,R2-R3,add)	    
	   end if; end if
	   end if
       
	   if(Vsource=='C')then
       !project dC to dP and dD
	   add=factor*dCR(idim,jdim,iR)
       R1=ra(:,iatom); R3=R1+ri; R4=S+ra(:,jatom); R2=R4+rj
	   if(near(R2-R1).and.near(R3-R4))then
	     Rc2c=S; orb=(/iorb1,jorb2,jorb1,iorb2/)
	     if(Pcoupled)call addto(dPc,orb,R2-R1,Rc2c,R3-R4,add)
	   end if
	   if(near(R4-R1).and.near(R2-R3))then
	     iiatom=iatom; jjatom=findatom(R3)
	     Rc2c= R3-ra(:,jjatom)-(R1-ra(:,iiatom)) 
         orb=(/iorb1,jorb1,iorb2,jorb2/)
	     call addto(dDc,orb,R4-R1,Rc2c,R2-R3,add)
	   end if
	   end if

       if(Vsource=='D')then
	   !project dD to dP and dC
	   add=factor*dDR(idim,jdim,iR)
	   R1=ra(:,iatom); R4=R1+ri; R3=S+ra(:,jatom); R2=R3+rj
	   if(near(R2-R1).and.near(R3-R4))then
	     iiatom=iatom; jjatom=findatom(R4)
	     Rc2c=R4-ra(:,jjatom)-(R1-ra(:,iiatom))
         orb=(/iorb1,jorb2,iorb2,jorb1/)
	     if(Pcoupled)call addto(dPc,orb,R2-R1,Rc2c,R3-R4,add)
	   end if
	   if(near(R3-R1).and.near(R2-R4))then
	     iiatom=iatom; jjatom=findatom(R4)
	     Rc2c=R4-ra(:,jjatom)-(R1-ra(:,iiatom))
         orb=(/iorb1,jorb1,iorb2,jorb2/)
	     call addto(dCc,orb,R3-R1,Rc2c,R2-R4,add)
	   end if
	   end if

       end do; end do; end do; end do
    !finished for specified source
    return
end subroutine projection





 subroutine contact_clever(dP,dC,dD)
    use workplace
    implicit none

    complex*16, dimension (ndim,ndim,nq) :: dP,dC,dD
    complex*16, dimension (ndim,ndim,nR) :: dPR,dCR,dDR,dPc,dCc,dDc

    integer iq,iR
	integer idim,iLform,nbi,i,iorb1,iorb2
	integer jdim,jLform,nbj,j,jorb1,jorb2
    integer iatom,iiatom,jatom,jjatom
	integer orb(4)
	integer findatom

    type (Mformconfig), pointer :: Mi,Mj
	type (Lformconfig), pointer :: Li,Lj

    real*8, dimension (3) :: R1,R2,R3,R4,ri,rj,S,Sc2c
	integer Rc2c(2)
    real*8 fLi,fLj
	real*8, pointer, dimension (:,:) :: ra

    complex*16 add,factor,factorc
	logical near,Pcoupled

    Pcoupled=.true.

    ra=>model%ra

    dPR=0; dCR=0; dDR=0

    do iR=1,nR; S=R(1,iR)*model%a+R(2,iR)*model%b
	   do iq=1,nq; factor=exp(-one*pi*sum(qv(1:2,iq)*S(1:2)))*qv(3,iq)
       do jdim=1,ndim; rj=ra(:,Mformfunc(jdim)%atom )
	   do idim=1,ndim; ri=ra(:,Mformfunc(idim)%atom )
	  	  factorc = exp( -one*pi*sum( qv(1:2,iq)*(rj(1:2)-ri(1:2)) )  )  
	      dPR(idim,jdim,iR)=dPR(idim,jdim,iR)+factor*dP(idim,jdim,iq)*factorc
	      dCR(idim,jdim,iR)=dCR(idim,jdim,iR)+factor*dC(idim,jdim,iq)*factorc
	      dDR(idim,jdim,iR)=dDR(idim,jdim,iR)+factor*dD(idim,jdim,iq)*factorc
       end do; end do
	end do; end do

    dPc=0; dCc=0; dDc=0
		
	call addto_clever(dPR,nP2C,P2C,infoP2C,dCc)
    call addto_clever(dPR,nP2D,P2D,infoP2D,dDc)
	
	call addto_clever(dCR,nC2P,C2P,infoC2P,dPc)
	call addto_clever(dCR,nC2D,C2D,infoC2D,dDc)

	call addto_clever(dDR,nD2C,D2C,infoD2C,dCc)
	call addto_clever(dDR,nD2P,D2P,infoD2P,dPc)

    !FT dPc,dCc and dDc and add to dP,dC,dD
    do iq=1,nq; do iR=1,nR; S=R(1,iR)*model%a+R(2,iR)*model%b
	   factor=exp( one*pi*sum(qv(1:2,iq)*S(1:2)))
       do jdim=1,ndim; rj=ra(:,Mformfunc(jdim)%atom )
	   do idim=1,ndim; ri=ra(:,Mformfunc(idim)%atom )
	  	  factorc=exp( one*pi*sum( qv(1:2,iq)*(rj(1:2)-ri(1:2)) )  )  
	      dP(idim,jdim,iq)=dP(idim,jdim,iq) +  factor*dPc(idim,jdim,iR)*factorc
	      dC(idim,jdim,iq)=dC(idim,jdim,iq) +  factor*dCc(idim,jdim,iR)*factorc
	      dD(idim,jdim,iq)=dD(idim,jdim,iq) +  factor*dDc(idim,jdim,iR)*factorc
       end do; end do
    end do; end do

    return
end subroutine contact_clever

subroutine addto_clever(A,nA2B,A2B,infoA2B,B)
    use workplace
    implicit none

    complex*16, dimension (ndim,ndim,nR) :: A,B
	integer nA2B,infoA2B(6,nA2B)
	complex*16 A2B(nA2B)

    integer i
	integer idim,jdim,iR
	integer iidim,jjdim,iRR

	do i=1,nA2B
	   iidim=infoA2B(1,i); jjdim=infoA2B(2,i); iRR=infoA2B(3,i)
	   idim=infoA2B(4,i); jdim=infoA2B(5,i); iR=infoA2B(6,i)
	   B(iidim,jjdim,iRR)=B(iidim,jjdim,iRR)+A(idim,jdim,iR)*A2B(i)
	end do
    return
end subroutine addto_clever

subroutine projectors()
    use workplace
    implicit none

    complex*16, dimension (ndim,ndim,nR) :: dPR,dCR,dDR,dPc,dCc,dDc
    
	integer i
	integer idim,jdim,iR
	integer iidim,jjdim,iRR
	integer iP2C,iP2D,iC2P,iC2D,iD2C,iD2P
    real*8 eps
	logical isave

    if(AppendProjectors)then
      open(10,file=projectorfile)
	  print*,'reading projectors...'
	  read(10,*)nP2C,nP2D,nC2P,nC2D,nD2C,nD2P
      allocate(P2C(nP2C),infoP2C(6,nP2C),P2D(nP2D),infoP2D(6,nP2D))
      allocate(C2P(nC2P),infoC2P(6,nC2P),C2D(nC2D),infoC2D(6,nC2D))
      allocate(D2C(nD2C),infoD2C(6,nD2C),D2P(nD2P),infoD2P(6,nD2P))
	  read(10,*)infoP2C,infoP2D,infoC2P,infoC2D,infoD2C,infoD2P
	  read(10,*)P2C,P2D,C2P,C2D,D2C,D2P
	  close(10)
	  return
    end if


    open(11,file='work1.dat');	open(12,file='work2.dat')
    eps=1.e-6
	print*,'setting projectors...'
    !Projector from P to C and D
    iP2C=0; iP2D=0 
    do iR=1,nR; do jdim=1,ndim; do idim=1,ndim
       dPR=0; dCR=0; dDR=0; dPR(idim,jdim,iR)=1
       call projection(dPR,dCR,dDR,dPc,dCc,dDc,'P',iR,idim,jdim)
       do iRR=1,nR; do jjdim=1,ndim; do iidim=1,ndim
          if(abs(dCc(iidim,jjdim,iRR))>eps)then
		    iP2C=iP2C+1
			write(11,*)dCc(iidim,jjdim,iRR)
			write(11,*)iidim,jjdim,iRR,idim,jdim,iR
		   end if
	       if(abs(dDc(iidim,jjdim,iRR))>eps)then
	         iP2D=iP2D+1
			 write(12,*)dDc(iidim,jjdim,iRR)
			 write(12,*)iidim,jjdim,iRR,idim,jdim,iR
		   end if
       end do; end do; end do
	   !print*,iP2C,iP2D
	end do; end do; end do
    nP2C=iP2C; nP2D=iP2D
	allocate(P2C(nP2C),infoP2C(6,nP2C),P2D(nP2D),infoP2D(6,nP2D))
    rewind(11); rewind(12)
    do iP2C=1,nP2C; read(11,*)P2C(iP2C); read(11,*)infoP2C(:,iP2C); end do
	do iP2D=1,nP2D; read(12,*)P2D(iP2D); read(12,*)infoP2D(:,iP2D); end do
	print*,'nP2C,nP2D=',nP2C,nP2D

    !Projector from C to P and D
    rewind(11); rewind(12)
    iC2P=0; iC2D=0 
    do iR=1,nR; do jdim=1,ndim; do idim=1,ndim
       dPR=0; dCR=0; dDR=0; dCR(idim,jdim,iR)=1
       call projection(dPR,dCR,dDR,dPc,dCc,dDc,'C',iR,idim,jdim)
       do iRR=1,nR; do jjdim=1,ndim; do iidim=1,ndim
          if(abs(dPc(iidim,jjdim,iRR))>eps)then
 	        iC2P=iC2P+1
			write(11,*)dPc(iidim,jjdim,iRR)
			write(11,*)iidim,jjdim,iRR,idim,jdim,iR
	       end if
	       if(abs(dDc(iidim,jjdim,iRR))>eps)then
  	         iC2D=iC2D+1
			 write(12,*)dDc(iidim,jjdim,iRR)
			 write(12,*)iidim,jjdim,iRR,idim,jdim,iR
	       end if
       end do; end do; end do
	end do; end do; end do
    nC2P=iC2P; nC2D=iC2D
	allocate(C2P(nC2P),infoC2P(6,nC2P),C2D(nC2D),infoC2D(6,nC2D))
    rewind(11); rewind(12)
    do iC2P=1,nC2P; read(11,*)C2P(iC2P); read(11,*)infoC2P(:,iC2P); end do
	do iC2D=1,nC2D; read(12,*)C2D(iC2D); read(12,*)infoC2D(:,iC2D); end do
	print*,'nC2P,nC2D=',nC2P,nC2D
     

    !Projector from D to C and P
    rewind(11); rewind(12)
    iD2C=0; iD2P=0 
    do iR=1,nR; do jdim=1,ndim; do idim=1,ndim
       dPR=0; dCR=0; dDR=0; dDR(idim,jdim,iR)=1
       call projection(dPR,dCR,dDR,dPc,dCc,dDc,'D',iR,idim,jdim)
       do iRR=1,nR; do jjdim=1,ndim; do iidim=1,ndim
          if(abs(dCc(iidim,jjdim,iRR))>eps)then
  	        iD2C=iD2C+1
			write(11,*)dCc(iidim,jjdim,iRR)
			write(11,*)iidim,jjdim,iRR,idim,jdim,iR
		  end if
	      if(abs(dPc(iidim,jjdim,iRR))>eps)then
		    iD2P=iD2P+1
			write(12,*)dPc(iidim,jjdim,iRR)
 	        write(12,*)iidim,jjdim,iRR,idim,jdim,iR
 	      end if
       end do; end do; end do
	end do; end do; end do
    nD2C=iD2C; nD2P=iD2P
    allocate(D2C(nD2C),infoD2C(6,nD2C),D2P(nD2P),infoD2P(6,nD2P))
    rewind(11); rewind(12)
    do iD2C=1,nD2C; read(11,*)D2C(iD2C); read(11,*)infoD2C(:,iD2C); end do
	do iD2P=1,nD2P; read(12,*)D2P(iD2P); read(12,*)infoD2P(:,iD2P); end do
	print*,'nD2C,nD2P=',nD2C,nD2P

    rewind(11); rewind(12)
	write(11,*)'trashed file'; write(12,*)'trashed file'
	close(11); close(12)

    open(10,file=projectorfile)
	write(10,*)nP2C,nP2D,nC2P,nC2D,nD2C,nD2P
	write(10,*)infoP2C,infoP2D,infoC2P,infoC2D,infoD2C,infoD2P
	write(10,*)P2C,P2D,C2P,C2D,D2C,D2P
	close(10)
	print*,'projectors saved to file.'
 
    return
end subroutine projectors

