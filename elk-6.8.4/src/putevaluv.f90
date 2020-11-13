
! Copyright (C) 2019 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putevaluv(ik,evaluvp)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
real(8), intent(in) :: evaluvp(nstsv)
! local variables
integer recl
! find the record length
inquire(iolength=recl) vkl(:,ik),nstsv,evaluvp
!$OMP CRITICAL(u274)
open(274,file='EVALUV.OUT',form='UNFORMATTED',access='DIRECT',action='WRITE', &
 recl=recl)
write(274,rec=ik) vkl(:,ik),nstsv,evaluvp
close(274)
!$OMP END CRITICAL(u274)
end subroutine

