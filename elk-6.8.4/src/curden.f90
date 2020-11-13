
! Copyright (C) 2018 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine curden(afield)
use modmain
use modmpi
use modomp
implicit none
! arguments
real(8), intent(in) :: afield(3)
! local variables
integer ik,is,ias,ispn
integer nr,nri,iro,np
integer nrc,nrci,npc
integer ir,n,i,nthd
real(8) ca,t1
! automatic arrays
integer(8) lock(natmtot)
real(8) rfmt(npmtmax)
! external functions
real(8), external :: rfint
! coupling constant of the external A-field (-1/c)
ca=-1.d0/solsc
! set the current density to zero
do i=1,3
  do ias=1,natmtot
    is=idxis(ias)
    cdmt(1:npcmt(is),ias,i)=0.d0
  end do
end do
cdir(1:ngtc,:)=0.d0
! current density cannot be computed if wavefunctions do not exist
if (iscl.le.0) return
! initialise the OpenMP locks
do ias=1,natmtot
  call omp_init_lock(lock(ias))
end do
call holdthd(nkpt/np_mpi,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
do ik=1,nkpt
! distribute among MPI processes
  if (mod(ik-1,np_mpi).ne.lp_mpi) cycle
  call curdenk(ik,lock)
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! destroy the OpenMP locks
do ias=1,natmtot
  call omp_destroy_lock(lock(ias))
end do
! convert muffin-tin current density to spherical harmonics
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(rfmt,is,nrc,nrci,npc,i) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
  do i=1,3
    rfmt(1:npc)=cdmt(1:npc,ias,i)
    call rfsht(nrc,nrci,rfmt,cdmt(:,ias,i))
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
! symmetrise the current density
call symrvf(.false.,.true.,nrcmt,nrcmti,npcmt,ngdgc,ngtc,ngvc,igfc,npmtmax, &
 cdmt,ngtot,cdir)
! convert muffin-tin and interstitial current density from coarse to fine grids
call holdthd(3,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do i=1,3
  call rfmtctof(cdmt(:,:,i))
end do
!$OMP END DO NOWAIT
!$OMP DO
do i=1,3
  call rfirctof(cdir(:,i),cdir(:,i))
end do
!$OMP END DO
!$OMP END PARALLEL
call freethd(nthd)
! add current densities from each process and redistribute
if (np_mpi.gt.1) then
  n=npmtmax*natmtot*3
  call mpi_allreduce(mpi_in_place,cdmt,n,mpi_double_precision,mpi_sum,mpicom, &
   ierror)
  n=ngtot*3
  call mpi_allreduce(mpi_in_place,cdir,n,mpi_double_precision,mpi_sum,mpicom, &
   ierror)
end if
! add vector potential contribution to make current gauge invariant
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(rfmt,is,nr,nri,np) &
!$OMP PRIVATE(iro,ispn,i,ir,t1) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  np=npmt(is)
  iro=nri+1
! remove the core density from the muffin-tin density
  rfmt(1:np)=rhomt(1:np,ias)
  do ispn=1,nspncr
    i=1
    do ir=1,nri
      rfmt(i)=rfmt(i)-rhocr(ir,ias,ispn)
      i=i+lmmaxi
    end do
    do ir=iro,nr
      rfmt(i)=rfmt(i)-rhocr(ir,ias,ispn)
      i=i+lmmaxo
    end do
  end do
  do i=1,3
    t1=ca*afield(i)
    call daxpy(np,t1,rfmt,1,cdmt(:,ias,i),1)
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
do i=1,3
  t1=ca*afield(i)
  call daxpy(ngtot,t1,rhoir,1,cdir(:,i),1)
end do
! compute the total current in the unit cell
do i=1,3
  curtot(i)=rfint(cdmt(:,:,i),cdir(:,i))
end do
! total current magnitude
curtotm=sqrt(curtot(1)**2+curtot(2)**2+curtot(3)**2)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
end subroutine

