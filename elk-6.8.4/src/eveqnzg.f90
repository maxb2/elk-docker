
! Copyright (C) 2019 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine eveqnzg(n,ld,a,w)
use modomp
implicit none
! arguments
integer, intent(in) :: n,ld
complex(8), intent(inout) :: a(ld,n)
complex(8), intent(out) :: w(n)
! local variables
integer lwork,info,nthd
! allocatable arrays
real(8), allocatable :: rwork(:)
complex(8), allocatable :: vr(:,:),work(:)
lwork=2*n
allocate(rwork(2*n),vr(n,n),work(lwork))
! enable MKL parallelism
call holdthd(maxthdmkl,nthd)
call mkl_set_num_threads(nthd)
! determine the eigenvalues and right eigenvectors
call zgeev('N','V',n,a,ld,w,vr,1,vr,n,work,lwork,rwork,info)
call freethd(nthd)
call mkl_set_num_threads(1)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(eveqnzg): diagonalisation failed")')
  write(*,'(" ZGEEV returned INFO = ",I8)') info
  write(*,*)
  stop
end if
! copy right eigenvectors to output array
a(1:n,1:n)=vr(1:n,1:n)
deallocate(rwork,vr,work)
end subroutine

