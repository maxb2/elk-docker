
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

complex(8) function zmdet(n,a)
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(in) :: a(n,n)
! local variables
integer i,m,info
! automatic arrays
integer ipiv(n)
call zgetrf(n,n,a,n,ipiv,info)
! multiply diagonal elements of U together
zmdet=a(1,1)
do i=2,n
  zmdet=zmdet*a(i,i)
end do
! determine the sign from the number of row interchanges
m=1
do i=1,n
  if (ipiv(i).ne.i) m=-m
end do
if (m.eq.-1) zmdet=-zmdet
end function

