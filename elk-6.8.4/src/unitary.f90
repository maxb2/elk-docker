
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine unitary(n,a)
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(inout) :: a(n,n)
! local variables
integer lwork,info
complex(8), parameter :: zzero=(0.d0,0.d0)
complex(8), parameter :: zone=(1.d0,0.d0)
! allocatable arrays
real(8), allocatable :: s(:),rwork(:)
complex(8), allocatable :: u(:,:),vt(:,:),work(:)
! perform singular value decomposition on matrix
lwork=3*n
allocate(s(n),rwork(5*n))
allocate(u(n,n),vt(n,n),work(lwork))
call zgesvd('A','A',n,n,a,n,s,u,n,vt,n,work,lwork,rwork,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(unitary): singular value decomposition failed")')
  write(*,'(" ZGESVD returned INFO = ",I8)') info
  write(*,*)
  stop
end if
! multiply the two unitary matrices together and store in the input matrix
call zgemm('N','N',n,n,n,zone,u,n,vt,n,zzero,a,n)
deallocate(s,rwork,u,vt,work)
end subroutine
