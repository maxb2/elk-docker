
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rhomagq
use modmain
use modulr
use modomp
implicit none
! local variables
integer ifq,idm,is,ias,nthd
! automatic arrays
complex(8) zfmt(npcmtmax)
! partial Fourier transform of density to Q-space
call rfzfftq(-1,1,ngtc,rhormt,rhorir,rhoqmt,rhoqir)
! convert density to spherical harmonics
call holdthd(nfqrz,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(zfmt,ias,is) &
!$OMP NUM_THREADS(nthd)
do ifq=1,nfqrz
  do ias=1,natmtot
    is=idxis(ias)
    call zcopy(npcmt(is),rhoqmt(:,ias,ifq),1,zfmt,1)
    call zfsht(nrcmt(is),nrcmti(is),zfmt,rhoqmt(:,ias,ifq))
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
if (.not.spinpol) return
! partial Fourier transform of magnetisation to Q-space
call rfzfftq(-1,ndmag,ngtc,magrmt,magrir,magqmt,magqir)
! convert magnetisation to spherical harmonics
call holdthd(nfqrz,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(zfmt,idm,ias,is) &
!$OMP NUM_THREADS(nthd)
do ifq=1,nfqrz
  do idm=1,ndmag
    do ias=1,natmtot
      is=idxis(ias)
      call zcopy(npcmt(is),magqmt(:,ias,idm,ifq),1,zfmt,1)
      call zfsht(nrcmt(is),nrcmti(is),zfmt,magqmt(:,ias,idm,ifq))
    end do
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
end subroutine

