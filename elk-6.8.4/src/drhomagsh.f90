
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine drhomagsh
use modmain
use modphonon
use modomp
implicit none
! local variables
integer idm,is,ias,nthd
integer nrc,nrci,npc
! allocatable arrays
complex(8), allocatable :: zfmt(:)
call holdthd(natmtot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zfmt,is,nrc,nrci,npc,idm) &
!$OMP NUM_THREADS(nthd)
allocate(zfmt(npcmtmax))
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
! convert the density derivative to spherical harmonics
  zfmt(1:npc)=drhomt(1:npc,ias)
  call zfsht(nrc,nrci,zfmt,drhomt(:,ias))
end do
!$OMP END DO NOWAIT
do idm=1,ndmag
!$OMP DO
  do ias=1,natmtot
    is=idxis(ias)
    nrc=nrcmt(is)
    nrci=nrcmti(is)
    npc=npcmt(is)
! convert the magnetisation derivative to spherical harmonics
    zfmt(1:npc)=dmagmt(1:npc,ias,idm)
    call zfsht(nrc,nrci,zfmt,dmagmt(:,ias,idm))
  end do
!$OMP END DO
end do
deallocate(zfmt)
!$OMP END PARALLEL
call freethd(nthd)
end subroutine

