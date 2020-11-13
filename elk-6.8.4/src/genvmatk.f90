
! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genvmatk(vmt,vir,ngp,igpig,wfmt,ld,wfgp,vmat)
use modmain
use modomp
implicit none
! arguments
! the potential is multiplied by the radial integration weights in the
! muffin-tin and by the characteristic function in the interstitial region
real(8), intent(in) :: vmt(npcmtmax,natmtot),vir(ngtot)
integer, intent(in) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
complex(8), intent(in) :: wfmt(npcmtmax,natmtot,nspinor,nstsv)
integer, intent(in) :: ld
complex(8), intent(in) :: wfgp(ld,nspinor,nstsv)
complex(8), intent(out) :: vmat(nstsv,nstsv)
! local variables
integer ist,jst,ispn,jspn
integer is,ias,npc,igp,nthd
! automatic arrays
complex(8) wfmt1(npcmtmax),z(ngkmax)
! allocatable arrays
complex(8), allocatable :: wfir(:)
! external functions
complex(8), external :: zdotc
! zero the matrix elements
vmat(:,:)=0.d0
!-------------------------!
!     muffin-tin part     !
!-------------------------!
call holdthd(nstsv,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(wfmt1,ispn,ias) &
!$OMP PRIVATE(is,npc,ist) &
!$OMP NUM_THREADS(nthd)
do jst=1,nstsv
  do ispn=1,nspinor
    do ias=1,natmtot
      is=idxis(ias)
      npc=npcmt(is)
! apply potential to wavefunction
      wfmt1(1:npc)=vmt(1:npc,ias)*wfmt(1:npc,ias,ispn,jst)
! compute the inner products
      do ist=1,jst
        vmat(ist,jst)=vmat(ist,jst)+zdotc(npc,wfmt(:,ias,ispn,ist),1,wfmt1,1)
      end do
    end do
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
!---------------------------!
!     interstitial part     !
!---------------------------!
call holdthd(nstsv,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfir,z,ispn,jspn) &
!$OMP PRIVATE(igp,ist) &
!$OMP NUM_THREADS(nthd)
allocate(wfir(ngtot))
!$OMP DO
do jst=1,nstsv
  do ispn=1,nspinor
    jspn=jspnfv(ispn)
! Fourier transform wavefunction to real-space
    wfir(:)=0.d0
    do igp=1,ngp(jspn)
      wfir(igfft(igpig(igp,jspn)))=wfgp(igp,ispn,jst)
    end do
    call zfftifc(3,ngridg,1,wfir)
! apply potential to wavefunction
    wfir(1:ngtot)=vir(1:ngtot)*wfir(1:ngtot)
! Fourier transform to G+p-space
    call zfftifc(3,ngridg,-1,wfir)
    do igp=1,ngp(jspn)
      z(igp)=wfir(igfft(igpig(igp,jspn)))
    end do
    do ist=1,jst
! compute inner product
      vmat(ist,jst)=vmat(ist,jst)+zdotc(ngp(jspn),wfgp(:,ispn,ist),1,z,1)
    end do
  end do
end do
!$OMP END DO
deallocate(wfir)
!$OMP END PARALLEL
call freethd(nthd)
! lower triangular part
do ist=1,nstsv
  do jst=1,ist-1
    vmat(ist,jst)=conjg(vmat(jst,ist))
  end do
end do
end subroutine

