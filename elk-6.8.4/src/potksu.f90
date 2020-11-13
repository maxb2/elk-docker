
! Copyright (C) 2018 T. Mueller, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine potksu
use modmain
use modulr
implicit none
! compute the ultra long-range Coulomb potential
call potcoulu
! compute the ultra long-range exchange-correlation potential and fields
call potxcu
! reduce the external magnetic field if required
if (reducebf.lt.1.d0) then
  bfcq(:,:)=bfcq(:,:)*reducebf
  bfcmtq(:,:,:)=bfcmtq(:,:,:)*reducebf
end if
end subroutine

