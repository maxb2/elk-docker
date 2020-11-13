
! Copyright (C) 2019 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putevecuv(ik,evecu,evecv)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: evecu(nstsv,nstsv),evecv(nstsv,nstsv)
! local variables
integer recl
! find the record length
inquire(iolength=recl) vkl(:,ik),nstsv,evecu,evecv
!$OMP CRITICAL(u276)
open(276,file='EVECUV.OUT',form='UNFORMATTED',access='DIRECT',action='WRITE',&
 recl=recl)
write(276,rec=ik) vkl(:,ik),nstsv,evecu,evecv
close(276)
!$OMP END CRITICAL(u276)
end subroutine

