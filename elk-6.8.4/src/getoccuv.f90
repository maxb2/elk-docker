
! Copyright (C) 2019 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine getoccuv(ik,occuvp)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
real(8), intent(out) :: occuvp(nstsv)
! local variables
integer recl,nstsv_
real(8) vkl_(3),t1
! find the record length
inquire(iolength=recl) vkl_,nstsv_,occuvp
!$OMP CRITICAL(u278)
open(278,file='OCCUV.OUT',form='UNFORMATTED',access='DIRECT',action='READ', &
 recl=recl)
read(278,rec=ik) vkl_,nstsv_,occuvp
close(278)
!$OMP END CRITICAL(u278)
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getoccuv): differing vectors for k-point ",I8)') ik
  write(*,'(" current   : ",3G18.10)') vkl(:,ik)
  write(*,'(" OCCUV.OUT : ",3G18.10)') vkl_
  write(*,*)
  stop
end if
if (nstsv.ne.nstsv_) then
  write(*,*)
  write(*,'("Error(getoccuv): differing nstsv for k-point ",I8)') ik
  write(*,'(" current   : ",I8)') nstsv
  write(*,'(" OCCUV.OUT : ",I8)') nstsv_
  write(*,*)
  stop
end if
end subroutine
