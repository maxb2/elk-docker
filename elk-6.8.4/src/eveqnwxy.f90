
! Copyright (C) 2019 Chung-Yu Wang, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine eveqnwxy(n,dw,ex,fy,w)
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(inout) :: dw(n,n),ex(n,n),fy(n)
complex(8), intent(out) :: w(n)
! local variables
integer n2,i,j
real(8) t1
complex(8), parameter :: zzero=(0.d0,0.d0),zone=(1.d0,0.d0)
! allocatable arrays
integer, allocatable :: idx(:)
real(8), allocatable :: r(:)
complex(8), allocatable :: w2(:),h(:,:)
complex(8), allocatable :: x(:),a(:,:)
n2=2*n
! setup the Bogoliubov Hamiltonian
allocate(w2(n2),h(n2,n2))
do j=1,n
  do i=1,n
    h(i,j)=dw(i,j)
    h(n+i,n+j)=-conjg(dw(i,j))
    h(i,n+j)=-ex(i,j)
    h(n+i,j)=conjg(ex(i,j))
  end do
end do
! find the eigenvalues and right eigenvectors
call eveqnzg(n2,n2,h,w2)
! select the positive eigenvalues and corresponding eigenvectors
allocate(idx(n2),r(n2))
do i=1,n2
  r(i)=dble(w2(i))
end do
call sortidx(n2,r,idx)
do i=1,n
  j=idx(n+i)
  w(i)=w2(j)
  call zcopy(n,h(1,j),1,dw(1,i),1)
  call zcopy(n,h(n+1,j),1,ex(1,i),1)
end do
deallocate(idx,r,w2,h)
! solve for the vector y
allocate(x(n),a(n,n))
a(:,:)=dw(:,:)-ex(:,:)
x(:)=fy(:)
call zgemv('T',n,n,zone,a,n,x,1,zzero,fy,1)
do i=1,n
  t1=abs(dble(w(i)))+abs(aimag(w(i)))
  if (t1.gt.1.d-6) then
    fy(i)=fy(i)/w(i)
  else
    fy(i)=0.d0
  end if
end do
deallocate(x,a)
end subroutine

