
! Copyright (C) 2019 Chung-Yu Wang, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine eveqnuv(n,au,bv,w,un)
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(inout) :: au(n,n),bv(n,n)
real(8), intent(out) :: w(n),un(n)
! local variables
integer n2,i,j
! allocatable arrays
real(8), allocatable :: w2(:)
complex(8), allocatable :: h(:,:)
! external functions
real(8), external :: dznrm2
n2=2*n
! setup the Bogoliubov Hamiltonian
allocate(w2(n2),h(n2,n2))
do j=1,n
  do i=1,n
    h(i,j)=au(i,j)
    h(n+i,n+j)=-conjg(au(i,j))
    h(i,n+j)=bv(i,j)
  end do
end do
! find the eigenvalues and eigenvectors
call eveqnzh(n2,n2,h,w2)
! copy to output arrays
do i=1,n
! choose the positive eigenvalues
  j=n+i
  w(i)=w2(j)
  call zcopy(n,h(1,j),1,au(1,i),1)
  call zcopy(n,h(n+1,j),1,bv(1,i),1)
! store the U-norm
  un(i)=dznrm2(n,au(1,i),1)**2
end do
deallocate(w2,h)
end subroutine

