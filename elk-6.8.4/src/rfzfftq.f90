
! Copyright (C) 2019 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rfzfftq(sgn,nf,ngt,rfmt,rfir,zfmt,zfir)
use modmain
use modmpi
use modomp
implicit none
! arguments
integer, intent(in) :: sgn,nf,ngt
real(8), intent(inout) :: rfmt(npcmtmax,natmtot,nf,nqpt)
real(8), intent(inout) :: rfir(ngt,nf,nqpt)
complex(8), intent(inout) :: zfmt(npcmtmax,natmtot,nf,nfqrz)
complex(8), intent(inout) :: zfir(ngt,nf,nfqrz)
! local variables
integer jf,is,ias,ir
integer npc,i,n,nthd
! automatic arrays
real(8) r(nqpt)
complex(8) z(nfqrz)
if (sgn.eq.-1) then
! loop over the number of functions
  do jf=1,nf
! Fourier transform the muffin-tin function
    do ias=1,natmtot
      is=idxis(ias)
      npc=npcmt(is)
      if (np_mpi.gt.1) zfmt(1:npc,ias,jf,:)=0.d0
      call holdthd(npc/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(r,z) &
!$OMP NUM_THREADS(nthd)
      do i=1,npc
! distribute among MPI processes
        if (mod(i-1,np_mpi).ne.lp_mpi) cycle
        r(1:nqpt)=rfmt(i,ias,jf,1:nqpt)
        call rzfftifc(3,ngridq,-1,r,z)
        zfmt(i,ias,jf,1:nfqrz)=z(1:nfqrz)
      end do
!$OMP END PARALLEL DO
      call freethd(nthd)
    end do
! Fourier transform the interstitial function
    if (np_mpi.gt.1) zfir(:,jf,:)=0.d0
    call holdthd(ngt/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(r,z) &
!$OMP NUM_THREADS(nthd)
    do ir=1,ngt
! distribute among MPI processes
      if (mod(ir-1,np_mpi).ne.lp_mpi) cycle
      r(1:nqpt)=rfir(ir,jf,1:nqpt)
      call rzfftifc(3,ngridq,-1,r,z)
      zfir(ir,jf,1:nfqrz)=z(1:nfqrz)
    end do
!$OMP END PARALLEL DO
    call freethd(nthd)
! end loop over number of functions
  end do
  if (np_mpi.gt.1) then
    n=npcmtmax*natmtot*nf*nfqrz
    call mpi_allreduce(mpi_in_place,zfmt,n,mpi_double_complex,mpi_sum,mpicom, &
     ierror)
    n=ngt*nf*nfqrz
    call mpi_allreduce(mpi_in_place,zfir,n,mpi_double_complex,mpi_sum,mpicom, &
     ierror)
  end if
else
! loop over the number of functions
  do jf=1,nf
    do ias=1,natmtot
      is=idxis(ias)
      npc=npcmt(is)
      if (np_mpi.gt.1) rfmt(1:npc,ias,jf,:)=0.d0
      call holdthd(npc/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(z,r) &
!$OMP NUM_THREADS(nthd)
      do i=1,npc
! distribute among MPI processes
        if (mod(i-1,np_mpi).ne.lp_mpi) cycle
        z(1:nfqrz)=zfmt(i,ias,jf,1:nfqrz)
        call rzfftifc(3,ngridq,1,r,z)
        rfmt(i,ias,jf,1:nqpt)=r(1:nqpt)
      end do
!$OMP END PARALLEL DO
      call freethd(nthd)
    end do
    if (np_mpi.gt.1) rfir(:,jf,:)=0.d0
    call holdthd(ngt/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(z,r) &
!$OMP NUM_THREADS(nthd)
    do ir=1,ngt
! distribute among MPI processes
      if (mod(ir-1,np_mpi).ne.lp_mpi) cycle
      z(1:nfqrz)=zfir(ir,jf,1:nfqrz)
      call rzfftifc(3,ngridq,1,r,z)
      rfir(ir,jf,1:nqpt)=r(1:nqpt)
    end do
!$OMP END PARALLEL DO
    call freethd(nthd)
! end loop over number of functions
  end do
  if (np_mpi.gt.1) then
    n=npcmtmax*natmtot*nf*nqpt
    call mpi_allreduce(mpi_in_place,rfmt,n,mpi_double_precision,mpi_sum, &
     mpicom,ierror)
    n=ngt*nf*nqpt
    call mpi_allreduce(mpi_in_place,rfir,n,mpi_double_precision,mpi_sum, &
     mpicom,ierror)
  end if
end if
end subroutine

