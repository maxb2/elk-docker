
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine grzvclncr(ias,gzfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: ias
complex(8), intent(out) :: gzfmt(npmtmax,3)
! local variables
integer is,nr,nri,ir,i
! automatic arrays
complex(8) zrhomt(npmtmax),zvclmt(npmtmax)
is=idxis(ias)
nr=nrmt(is)
nri=nrmti(is)
! convert the core density to complex spherical harmonics expansion
zrhomt(1:npmt(is))=0.d0
i=1
if (spincore) then
! spin-polarised core
  do ir=1,nri
    zrhomt(i)=rhocr(ir,ias,1)+rhocr(ir,ias,2)
    i=i+lmmaxi
  end do
  do ir=nri+1,nr
    zrhomt(i)=rhocr(ir,ias,1)+rhocr(ir,ias,2)
    i=i+lmmaxo
  end do
else
! spin-unpolarised core
  do ir=1,nri
    zrhomt(i)=rhocr(ir,ias,1)
    i=i+lmmaxi
  end do
  do ir=nri+1,nr
    zrhomt(i)=rhocr(ir,ias,1)
    i=i+lmmaxo
  end do
end if
! solve Poisson's equation in the muffin-tin
call zpotclmt(nr,nri,nrmtmax,rlmt(:,:,is),wprmt(:,:,is),zrhomt,zvclmt)
! add the nuclear Coulomb potential
i=1
do ir=1,nri
  zvclmt(i)=zvclmt(i)+vcln(ir,is)
  i=i+lmmaxi
end do
do ir=nri+1,nr
  zvclmt(i)=zvclmt(i)+vcln(ir,is)
  i=i+lmmaxo
end do
! compute the gradient of the potential
call gradzfmt(nr,nri,rlmt(:,-1,is),wcrmt(:,:,is),zvclmt,npmtmax,gzfmt)
end subroutine

