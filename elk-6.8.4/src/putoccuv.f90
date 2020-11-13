
! Copyright (C) 2019 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putoccuv(ik,occuvp)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
real(8), intent(in) :: occuvp(nstsv)
! local variables
integer recl
! find the record length
inquire(iolength=recl) vkl(:,ik),nstsv,occuvp
!$OMP CRITICAL(u278)
open(278,file='OCCUV.OUT',form='UNFORMATTED',access='DIRECT',action='WRITE', &
 recl=recl)
write(278,rec=ik) vkl(:,ik),nstsv,occuvp
close(278)
!$OMP END CRITICAL(u278)
end subroutine

