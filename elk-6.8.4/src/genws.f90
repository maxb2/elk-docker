
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genws
use modmain
use modomp
implicit none
! local variables
integer is,ias,nthd
integer nrc,nrci
! automatic arrays
real(8) rfmt(npcmtmax)
if (xcgrad.ne.4) return
! muffin-tin effective meta-GGA potential
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(rfmt,is,nrc,nrci) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
! convert to coarse radial mesh and spherical coordinates
  call rfmtftoc(nrc,nrci,wxcmt(:,ias),rfmt)
  call rbsht(nrc,nrci,rfmt,wsmt(:,ias))
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! interstitial meta-GGA potential
wsir(1:ngtot)=wxcir(1:ngtot)*cfunir(1:ngtot)
end subroutine

