
! Copyright (C) 2019 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine puthdbulr(ik0)
use modmain
use modulr
implicit none
! arguments
integer, intent(in) :: ik0
! local variables
integer ik,ikk,ikpa
integer ist,recl
real(8) t1
! automatic arrays
integer idx(nstsv)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:),evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: wfmtk(:,:,:,:),wfgkk(:,:,:)
complex(8), allocatable :: hdb(:,:,:),ok(:,:),b(:,:)
! central k-point
ik=(ik0-1)*nkpa+1
! get the ground-state eigenvectors from file for central k-point
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
call getevecfv('.OUT',ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
call getevecsv('.OUT',ik,vkl(:,ik),evecsv)
! find the matching coefficients
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
call match(ngk(1,ik),vgkc(:,:,1,ik),gkc(:,1,ik),sfacgk(:,:,1,ik),apwalm)
! index to all states
do ist=1,nstsv
  idx(ist)=ist
end do
! calculate the wavefunctions for all states of the central k-point
allocate(wfmt(npcmtmax,natmtot,nspinor,nstsv),wfir(ngtot,nspinor,nstsv))
call genwfsv(.false.,.false.,nstsv,idx,ngridg,igfft,ngk(1,ik),igkig(:,1,ik), &
 apwalm,evecfv,evecsv,wfmt,ngtot,wfir)
! compute the diagonal blocks of the ultra long-range Hamiltonian in the basis
! of the states at the central k-point
allocate(wfmtk(npcmtmax,natmtot,nspinor,nstsv),wfgkk(ngkmax,nspinor,nstsv))
allocate(hdb(nstsv,nstsv,2:nkpa),ok(nstsv,nstsv),b(nstsv,nstsv))
do ikpa=2,nkpa
  ikk=(ik0-1)*nkpa+ikpa
  call getevecfv('.OUT',ikk,vkl(:,ikk),vgkl(:,:,:,ikk),evecfv)
  call getevecsv('.OUT',ikk,vkl(:,ikk),evecsv)
  call match(ngk(1,ikk),vgkc(:,:,1,ikk),gkc(:,1,ikk),sfacgk(:,:,1,ikk),apwalm)
  call genwfsv(.false.,.true.,nstsv,idx,ngridg,igfft,ngk(:,ikk),igkig(:,:,ikk),&
   apwalm,evecfv,evecsv,wfmtk,ngkmax,wfgkk)
! compute the overlap matrix between the states at k and k+kappa
  call genolpq(nstsv,expqmt(:,:,ikpa),ngk(:,ikk),igkig(:,:,ikk),wfmt,wfir, &
   wfmtk,wfgkk,ok)
! use singular value decompostion to make the matrix strictly unitary
  call unitary(nstsv,ok)
! apply the overlap matrix from the right to the eigenvalues at k+kappa
  do ist=1,nstsv
    t1=evalsv(ist,ikk)
    b(ist,:)=t1*ok(ist,:)
  end do
! apply the conjugate transpose of the overlap matrix from the left to form the
! Hamiltonian matrix in the basis of states at k
  call zgemm('C','N',nstsv,nstsv,nstsv,zone,ok,nstsv,b,nstsv,zzero, &
   hdb(:,:,ikpa),nstsv)
end do
! determine the record length
inquire(iolength=recl) vkl(:,1),nstsv,nkpa,hdb
!$OMP CRITICAL(u190)
open(190,file='HDBULR.OUT',form='UNFORMATTED',access='DIRECT',action='WRITE', &
 recl=recl)
write(190,rec=ik0) vkl(:,ik),nstsv,nkpa,hdb
close(190)
!$OMP END CRITICAL(u190)
deallocate(apwalm,evecfv,evecsv)
deallocate(wfmt,wfir,wfmtk,wfgkk)
deallocate(hdb,ok,b)
end subroutine

