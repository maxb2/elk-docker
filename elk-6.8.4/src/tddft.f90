
! Copyright (C) 2014 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine tddft
use modmain
use modtddft
use moddftu
use modmpi
use modomp
use modstore
use modtest
implicit none
! local variables
logical exist
integer ik,itimes0
! allocatable arrays
real(8), allocatable :: evalfv(:,:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
if (tshift) then
  write(*,*)
  write(*,'("Error(tddft): use tshift = .false. for the ground-state run")')
  write(*,*)
  stop
end if
! average force can be non-zero (allow for translation of atomic basis)
tfav0_=tfav0
tfav0=.false.
tforce_=tforce
! Ehrenfest dynamics
if (any(task.eq.[462,463])) then
! forces should not be calculated
  tforce=.false.
! enable small amplitude displacements
  tdatpos=.true.
! zero the displacements and velocities
  datposc(:,:,:,:)=0.d0
end if
! initialise global variables
call init0
call init1
! read the charge density and potentials from file
call readstate
! generate the first- and second-variational eigenvectors and eigenvalues for
! the k-point set reduced with the symmetries which leave A(t) invariant for all
! time steps
call genvsig
call gencore
call energykncr
call readfermi
call linengy
call genapwlofr
call gensocfr
call genevfsv
call occupy
! DFT+U
if (dftu.ne.0) then
  call gendmatmt
  call genvmatmt
end if
! generate the kinetic matrix elements in the second-variational basis
call genkmat(.false.,.false.)
! write the momentum matrix elements in the second-variational basis
call genpmat
! write the power density to file
if (mp_mpi) call writeafpdt
! copy EVALFV.OUT, EVECFV.OUT, OCCSV.OUT and EVECSV.OUT to _TD.OUT extension
if (mp_mpi.and.((task.eq.460).or.(task.eq.462))) then
  allocate(evalfv(nstfv,nspnfv),evecfv(nmatmax,nstfv,nspnfv))
  allocate(evecsv(nstsv,nstsv))
  do ik=1,nkpt
    call getevalfv('.OUT',ik,vkl(:,ik),evalfv)
    call putevalfv('_TD.OUT',ik,evalfv)
    call getevecfv('.OUT',ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
    call putevecfv('_TD.OUT',ik,evecfv)
    call putoccsv('_TD.OUT',ik,occsv(:,ik))
    call getevecsv('.OUT',ik,vkl(:,ik),evecsv)
! randomise eigenvectors at t=0 if required
    call rndevsv(rndevt0,evecsv)
    call putevecsv('_TD.OUT',ik,evecsv)
  end do
  deallocate(evalfv,evecfv,evecsv)
end if
! set global file extension
filext='_TD.OUT'
! output the new k-point set to file
if (mp_mpi) call writekpts
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
if (any(task.eq.[461,463])) then
! restart if required
  call tdrestart(itimes0)
else
! start from t=0
  itimes0=1
end if
! read the forces calculated during the previous TDDFT run
if (tdatpos) call readforcet
! read the atomic displacements and velocities
if (trddatpos) call readdatposc
! set the stop signal to .false.
tstop=.false.
!---------------------------------!
!    main loop over time steps    !
!---------------------------------!
if (mp_mpi) write(*,*)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
do itimes=itimes0,ntimes-1
  if (mp_mpi) then
    write(*,'("Info(tddft): time step ",I8," of ",I8,",   t = ",G18.10)') &
     itimes,ntimes-1,times(itimes)
  end if
! reset the OpenMP thread variables
  call omp_reset
! evolve the wavefunctions across a single time step
  call timestep
! generate the density and magnetisation at current time step
  call rhomag
! compute the total current
  call curden(afieldt(:,itimes))
! time step the induced A-field
  if (tafindt) call afindtstep
! calculate the electric field
  call genefieldt
! compute the time-dependent Kohn-Sham potentials and magnetic fields
  call potkst
! DFT+U
  if (dftu.ne.0) then
    call gendmatmt
    call genvmatmt
  end if
! compute the total energy
  call energytd
! write muffin-tin L, S and J if required
  if (tdlsj) call writetdlsj
! calculate the atomic forces if required
  if (tforce) then
    if ((itimes.eq.itimes0).or.(mod(itimes-1,ntsforce).eq.0)) then
      call force
    end if
  end if
! time step the atomic positions for Ehrenfest dynamics using forces calculated
! during the previous TDDFT run
  if (tdatpos) call atptstep
  if (mp_mpi) then
! write TDDFT output
    call writetddft
! check for STOP file
    inquire(file='STOP',exist=exist)
    if (exist) then
      open(50,file='STOP')
      close(50,status='DELETE')
      tstop=.true.
    end if
  end if
! broadcast tstop from master process to all other processes
  call mpi_bcast(tstop,1,mpi_logical,0,mpicom,ierror)
  if (tstop) exit
end do
filext='.OUT'
tfav0=tfav0_
tforce=tforce_
tdatpos=.false.
! write the total current of the last step to test file
call writetest(460,'total current of last time step',nv=3,tol=1.d-4,rva=curtot)
end subroutine

