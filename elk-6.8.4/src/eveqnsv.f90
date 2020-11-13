
! Copyright (C) 2002-2010 J. K. Dewhurst, S. Sharma, C. Ambrosch-Draxl,
! F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine eveqnsv(ngp,igpig,vgpc,apwalm,evalfv,evecfv,evalsvp,evecsv)
use modmain
use moddftu
use modomp
implicit none
! arguments
integer, intent(in) :: ngp,igpig(ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
real(8), intent(in) :: evalfv(nstfv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
real(8), intent(out) :: evalsvp(nstsv)
complex(8), intent(out) :: evecsv(nstsv,nstsv)
! local variables
logical socz
integer nsc,nsd,ld,ist,jst
integer ispn,jspn,is,ias
integer nrc,nrci,nrco,irco,irc
integer l,lm,nm,npc,npci,ipco
integer igp,i,j,k,nthd
real(8) a1,a2,a3,t1
real(8) ts0,ts1
complex(8) z1
! automatic arrays
complex(8) wfmt2(npcmtmax),wfmt3(npcmtmax)
complex(8) wfmt4(npcmtmax,3),wfmt5(npcmtmax,3)
complex(8) wfgp(ngkmax,3)
! allocatable arrays
complex(8), allocatable :: wfmt1(:,:),gwfmt(:,:,:)
complex(8), allocatable :: wfir1(:),wfir2(:),gwfgp(:,:)
! external functions
complex(8), external :: zdotc
! no calculation of second-variational eigenvectors
if (.not.tevecsv) then
  do i=1,nstsv
    evalsvp(i)=evalfv(i)
  end do
  evecsv(:,:)=0.d0
  do i=1,nstsv
    evecsv(i,i)=1.d0
  end do
  return
end if
call timesec(ts0)
if (tafield) then
! coupling constant of the external A-field (-1/c)
  t1=-1.d0/solsc
  a1=t1*afieldc(1); a2=t1*afieldc(2); a3=t1*afieldc(3)
end if
! number of spin combinations after application of Hamiltonian
if (spinpol) then
  if (ncmag.or.spinorb) then
    nsc=3
  else
    nsc=2
  end if
  nsd=2
else
  nsc=1
  nsd=1
end if
! special case of spin-orbit coupling and collinear magnetism
if (spinorb.and.cmagz) then
  socz=.true.
else
  socz=.false.
end if
ld=lmmaxdm*nspinor
! zero the second-variational Hamiltonian (stored in the eigenvector array)
evecsv(:,:)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
allocate(wfmt1(npcmtmax,nstfv))
if (xcgrad.eq.4) allocate(gwfmt(npcmtmax,3,nstfv))
call holdthd(nstfv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfmt2,wfmt3,wfmt4,wfmt5) &
!$OMP PRIVATE(ias,is,nrc,nrci,nrco,irco) &
!$OMP PRIVATE(npc,npci,ipco,ist,jst,i,j,k) &
!$OMP PRIVATE(irc,t1,l,nm,lm,ispn,jspn,z1) &
!$OMP NUM_THREADS(nthd)
! begin loop over atoms
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  nrco=nrc-nrci
  irco=nrci+1
  npc=npcmt(is)
  npci=npcmti(is)
  ipco=npci+1
! compute the first-variational wavefunctions
!$OMP DO
  do ist=1,nstfv
    call wavefmt(lradstp,ias,ngp,apwalm(:,:,:,ias),evecfv(:,ist),wfmt1(:,ist))
  end do
!$OMP END DO
! begin loop over states
!$OMP DO
  do jst=1,nstfv
    if (spinpol) then
! convert wavefunction to spherical coordinates
      call zbsht(nrc,nrci,wfmt1(:,jst),wfmt2)
! apply Kohn-Sham effective magnetic field
      wfmt3(1:npc)=bsmt(1:npc,ias,ndmag)*wfmt2(1:npc)
! convert to spherical harmonics and store in wfmt4
      call zfsht(nrc,nrci,wfmt3,wfmt4)
      wfmt4(1:npc,2)=-wfmt4(1:npc,1)
! non-collinear magnetic field
      if (socz) then
        wfmt4(1:npc,3)=0.d0
      else if (ncmag) then
        wfmt3(1:npc)=cmplx(bsmt(1:npc,ias,1),-bsmt(1:npc,ias,2),8)*wfmt2(1:npc)
        call zfsht(nrc,nrci,wfmt3,wfmt4(:,3))
      end if
! apply spin-orbit coupling if required
      if (spinorb) then
        call lopnzflm(lmaxi,nrci,lmmaxi,wfmt1(1,jst),wfmt5,wfmt5(1,2), &
         wfmt5(1,3))
        call lopnzflm(lmaxo,nrco,lmmaxo,wfmt1(ipco,jst),wfmt5(ipco,1), &
         wfmt5(ipco,2),wfmt5(ipco,3))
        i=1
! inner part of muffin-tin
        do irc=1,nrci
          t1=socfr(irc,ias)
          do lm=1,lmmaxi
            wfmt4(i,1)=wfmt4(i,1)+t1*wfmt5(i,3)
            wfmt4(i,2)=wfmt4(i,2)-t1*wfmt5(i,3)
            wfmt4(i,3)=wfmt4(i,3)+t1*(wfmt5(i,1) &
             +cmplx(aimag(wfmt5(i,2)),-dble(wfmt5(i,2)),8))
            i=i+1
          end do
        end do
! outer part of muffin-tin
        do irc=irco,nrc
          t1=socfr(irc,ias)
          do lm=1,lmmaxo
            wfmt4(i,1)=wfmt4(i,1)+t1*wfmt5(i,3)
            wfmt4(i,2)=wfmt4(i,2)-t1*wfmt5(i,3)
            wfmt4(i,3)=wfmt4(i,3)+t1*(wfmt5(i,1) &
             +cmplx(aimag(wfmt5(i,2)),-dble(wfmt5(i,2)),8))
            i=i+1
          end do
        end do
      end if
    else
      do k=1,nsc
        wfmt4(1:npc,k)=0.d0
      end do
    end if
! apply muffin-tin potential matrix if required
    if (tvmatmt) then
      do l=0,lmaxdm
        if (tvmmt(l,ias)) then
          nm=2*l+1
          lm=idxlm(l,-l)
          do k=1,nsc
            if (k.eq.1) then
              ispn=1
              jspn=1
            else if (k.eq.2) then
              ispn=2
              jspn=2
            else
              ispn=1
              jspn=2
            end if
            if (l.le.lmaxi) then
              call zgemm('N','N',nm,nrci,nm,zone,vmatmt(lm,ispn,lm,jspn,ias), &
               ld,wfmt1(lm,jst),lmmaxi,zone,wfmt4(lm,k),lmmaxi)
            end if
            i=npci+lm
            call zgemm('N','N',nm,nrco,nm,zone,vmatmt(lm,ispn,lm,jspn,ias),ld, &
             wfmt1(i,jst),lmmaxo,zone,wfmt4(i,k),lmmaxo)
          end do
        end if
      end do
    end if
! apply vector potential if required
    if (tafield) then
      call gradzfmt(nrc,nrci,rlcmt(:,-1,is),wcrcmt(:,:,is),wfmt1(:,jst), &
       npcmtmax,wfmt5)
      do i=1,npc
        z1=a1*wfmt5(i,1)+a2*wfmt5(i,2)+a3*wfmt5(i,3)
        z1=cmplx(aimag(z1),-dble(z1),8)
        wfmt4(i,1:nsd)=wfmt4(i,1:nsd)+z1
      end do
    end if
! apply the radial integral weights
    do k=1,nsc
      call zfmtwr(nrc,nrci,wrcmt(:,is),wfmt4(:,k))
    end do
! add to second-variational Hamiltonian matrix
! upper diagonal block
    do ist=1,jst
      evecsv(ist,jst)=evecsv(ist,jst)+zdotc(npc,wfmt1(:,ist),1,wfmt4,1)
    end do
! lower diagonal block
    if (nsc.ge.2) then
      j=jst+nstfv
      do ist=1,jst
        i=ist+nstfv
        evecsv(i,j)=evecsv(i,j)+zdotc(npc,wfmt1(:,ist),1,wfmt4(:,2),1)
      end do
    end if
! off-diagonal block
    if (nsc.eq.3) then
      do ist=1,nstfv
        evecsv(ist,j)=evecsv(ist,j)+zdotc(npc,wfmt1(:,ist),1,wfmt4(:,3),1)
      end do
    end if
! end loop over states
  end do
!$OMP END DO
! apply meta-GGA non-multiplicative potential if required
  if (xcgrad.eq.4) then
!$OMP DO
    do ist=1,nstfv
      call gradzfmt(nrc,nrci,rlcmt(:,-1,is),wcrcmt(:,:,is),wfmt1(:,ist), &
       npcmtmax,gwfmt(:,:,ist))
    end do
!$OMP END DO
!$OMP DO
    do jst=1,nstfv
      do k=1,3
        call zbsht(nrc,nrci,gwfmt(:,k,jst),wfmt2)
        wfmt2(1:npc)=wsmt(1:npc,ias)*wfmt2(1:npc)
        call zfsht(nrc,nrci,wfmt2,wfmt5(:,k))
        call zfmtwr(nrc,nrci,wrcmt(:,is),wfmt5(:,k))
      end do
      do ist=1,nstfv
        z1=0.d0
        do k=1,3
          z1=z1+zdotc(npc,gwfmt(:,k,ist),1,wfmt5(:,k),1)
        end do
        z1=0.5d0*z1
        evecsv(ist,jst)=evecsv(ist,jst)+z1
        if (spinpol) then
          i=ist+nstfv
          j=jst+nstfv
          evecsv(i,j)=evecsv(i,j)+z1
        end if
      end do
    end do
!$OMP END DO
  end if
! end loop over atoms
end do
!$OMP END PARALLEL
call freethd(nthd)
deallocate(wfmt1)
if (xcgrad.eq.4) deallocate(gwfmt)
!---------------------------!
!     interstitial part     !
!---------------------------!
if (spinpol.or.tafield.or.(xcgrad.eq.4)) then
  if (socz) nsc=2
  if (xcgrad.eq.4) allocate(gwfgp(ngkmax,nstfv))
  call holdthd(nstfv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfir1,wfir2,wfgp) &
!$OMP PRIVATE(ist,jst,igp,t1,z1,i,j,k) &
!$OMP NUM_THREADS(nthd)
  allocate(wfir1(ngtot),wfir2(ngtot))
! begin loop over states
!$OMP DO
  do jst=1,nstfv
    wfir1(:)=0.d0
    do igp=1,ngp
      wfir1(igfft(igpig(igp)))=evecfv(igp,jst)
    end do
! Fourier transform wavefunction to real-space
    call zfftifc(3,ngridg,1,wfir1)
! multiply with magnetic field and transform to G-space
    if (spinpol) then
      wfir2(:)=bsir(:,ndmag)*wfir1(:)
      call zfftifc(3,ngridg,-1,wfir2)
      do igp=1,ngp
        wfgp(igp,1)=wfir2(igfft(igpig(igp)))
      end do
      wfgp(1:ngp,2)=-wfgp(1:ngp,1)
      if (ncmag) then
        wfir2(:)=cmplx(bsir(:,1),-bsir(:,2),8)*wfir1(:)
        call zfftifc(3,ngridg,-1,wfir2)
        do igp=1,ngp
          wfgp(igp,3)=wfir2(igfft(igpig(igp)))
        end do
      end if
    else
      wfgp(1:ngp,1:nsd)=0.d0
    end if
! apply vector potential if required
    if (tafield) then
      wfir1(:)=0.d0
      do igp=1,ngp
        t1=a1*vgpc(1,igp)+a2*vgpc(2,igp)+a3*vgpc(3,igp)
        wfir1(igfft(igpig(igp)))=t1*evecfv(igp,jst)
      end do
      call zfftifc(3,ngridg,1,wfir1)
      wfir1(:)=wfir1(:)*cfunir(:)
      call zfftifc(3,ngridg,-1,wfir1)
      do igp=1,ngp
        z1=wfir1(igfft(igpig(igp)))
        wfgp(igp,1:nsd)=wfgp(igp,1:nsd)+z1
      end do
    end if
! add to second-variational Hamiltonian matrix
! upper diagonal block
    do ist=1,jst
      evecsv(ist,jst)=evecsv(ist,jst)+zdotc(ngp,evecfv(:,ist),1,wfgp,1)
    end do
! lower diagonal block
    if (nsc.ge.2) then
      j=jst+nstfv
      do ist=1,jst
        i=ist+nstfv
        evecsv(i,j)=evecsv(i,j)+zdotc(ngp,evecfv(:,ist),1,wfgp(:,2),1)
      end do
    end if
! off-diagonal block
    if (nsc.eq.3) then
      do ist=1,nstfv
        evecsv(ist,j)=evecsv(ist,j)+zdotc(ngp,evecfv(:,ist),1,wfgp(:,3),1)
      end do
    end if
! end loop over states
  end do
!$OMP END DO
! apply meta-GGA non-multiplicative potential if required
  if (xcgrad.eq.4) then
    do k=1,3
! determine the gradient of the wavefunctions
!$OMP DO
      do ist=1,nstfv
        do igp=1,ngp
          z1=evecfv(igp,ist)
          gwfgp(igp,ist)=vgpc(k,igp)*cmplx(-aimag(z1),dble(z1),8)
        end do
      end do
!$OMP END DO
!$OMP DO
      do jst=1,nstfv
        wfir1(:)=0.d0
        do igp=1,ngp
          wfir1(igfft(igpig(igp)))=gwfgp(igp,jst)
        end do
        call zfftifc(3,ngridg,1,wfir1)
        wfir1(:)=wsir(:)*wfir1(:)
        call zfftifc(3,ngridg,-1,wfir1)
        do igp=1,ngp
          wfgp(igp,1)=wfir1(igfft(igpig(igp)))
        end do
        do ist=1,nstfv
          z1=0.5d0*zdotc(ngp,gwfgp(:,ist),1,wfgp,1)
          evecsv(ist,jst)=evecsv(ist,jst)+z1
          if (spinpol) then
            i=ist+nstfv
            j=jst+nstfv
            evecsv(i,j)=evecsv(i,j)+z1
          end if
        end do
      end do
!$OMP END DO
    end do
  end if
  deallocate(wfir1,wfir2)
!$OMP END PARALLEL
  call freethd(nthd)
  if (xcgrad.eq.4) deallocate(gwfgp)
end if
! add the diagonal first-variational part
i=0
do ispn=1,nspinor
  do ist=1,nstfv
    i=i+1
    evecsv(i,i)=evecsv(i,i)+evalfv(ist)
  end do
end do
if (spcpl.or.(.not.spinpol)) then
! spins are coupled; or spin-unpolarised: full diagonalisation
  call eveqnzh(nstsv,nstsv,evecsv,evalsvp)
else
! spins not coupled: block diagonalise H
  call eveqnzh(nstfv,nstsv,evecsv,evalsvp)
  i=nstfv+1
  call eveqnzh(nstfv,nstsv,evecsv(i,i),evalsvp(i))
  do i=1,nstfv
    do j=1,nstfv
      evecsv(i,j+nstfv)=0.d0
      evecsv(i+nstfv,j)=0.d0
    end do
  end do
end if
call timesec(ts1)
!$OMP ATOMIC
timesv=timesv+ts1-ts0
end subroutine

