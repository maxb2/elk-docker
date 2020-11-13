
! Copyright (C) 2019 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine bandstruv
use modmain
use modbog
implicit none
! local variables
integer ip
! generate k-points along a path for band structure plots
call plotpt1d(bvec,nvp1d,npp1d,vvlp1d,vplp1d,dvp1d,dpp1d)
! loop over plot points along path
do ip=1,npp1d


end do
end subroutine

