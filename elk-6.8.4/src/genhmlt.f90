
! Copyright (C) 2014 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genhmlt(ik,vmt,vir,bmt,bir,kmat,pmat,h)
use modmain
use modtddft
use modmpi
implicit none
! arguments
integer, intent(in) :: ik
real(8), intent(in) :: vmt(npcmtmax,natmtot),vir(ngtot)
real(8), intent(in) :: bmt(npcmtmax,natmtot,ndmag),bir(ngtot,ndmag)
complex(8), intent(in) :: kmat(nstsv,nstsv),pmat(nstsv,nstsv,3)
complex(8), intent(out) :: h(nstsv,nstsv)
! local variables
integer ist,i
real(8) ca,t1
! automatic arrays
integer idx(nstsv)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:),evecfv(:,:),evecsv(:,:)
complex(8), allocatable :: wfmt(:,:,:,:),wfgk(:,:,:)
! coupling constant of the external A-field (-1/c)
ca=-1.d0/solsc
allocate(evecfv(nmatmax,nstfv),evecsv(nstsv,nstsv))
! get the ground-state eigenvectors from file for input k-point
call getevecfv('.OUT',ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
call getevecsv('.OUT',ik,vkl(:,ik),evecsv)
! find the matching coefficients
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
call match(ngk(1,ik),vgkc(:,:,1,ik),gkc(:,1,ik),sfacgk(:,:,1,ik),apwalm)
! index to all states
do ist=1,nstsv
  idx(ist)=ist
end do
! calculate the wavefunctions for all states of the input k-point
allocate(wfmt(npcmtmax,natmtot,nspinor,nstsv),wfgk(ngkmax,nspinor,nstsv))
call genwfsv(.false.,.true.,nstsv,idx,ngridg,igfft,ngk(:,ik),igkig(:,:,ik), &
 apwalm,evecfv,evecsv,wfmt,ngkmax,wfgk)
deallocate(apwalm,evecfv,evecsv)
! Kohn-Sham potential and magnetic field matrix elements
if (spinpol) then
  call genvbmatk(vmt,vir,bmt,bir,ngk(:,ik),igkig(:,:,ik),wfmt,ngkmax,wfgk,h)
else
  call genvmatk(vmt,vir,ngk(:,ik),igkig(:,:,ik),wfmt,ngkmax,wfgk,h)
end if
! add the DFT+U matrix elements if required
call genvmmtsv(wfmt,h)
deallocate(wfmt,wfgk)
! add the kinetic matrix elements in the second-variational basis
h(:,:)=h(:,:)+kmat(:,:)
! add the A-field matrix elements in the second-variational basis
do i=1,3
  t1=ca*afieldt(i,itimes)
  h(:,:)=h(:,:)+t1*pmat(:,:,i)
end do
end subroutine

