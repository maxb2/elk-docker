
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine putepsinv(iq,epsi)
use modmain
implicit none
! arguments
integer, intent(in) :: iq
complex(8), intent(in) :: epsi(ngrf,ngrf,nwrf)
! local variables
integer recl
! determine the record length for EPSINV.OUT
inquire(iolength=recl) vql(:,iq),ngrf,nwrf,epsi
!$OMP CRITICAL(u180)
open(180,file='EPSINV.OUT',form='UNFORMATTED',access='DIRECT',action='WRITE',&
 recl=recl)
write(180,rec=iq) vql(:,iq),ngrf,nwrf,epsi
close(180)
!$OMP END CRITICAL(u180)
end subroutine

