
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: mixpack
! !INTERFACE:
subroutine mixpack(tpack,n,v)
! !USES:
use modmain
use moddftu
! !INPUT/OUTPUT PARAMETERS:
!   tpack : .true. for packing, .false. for unpacking (in,logical)
!   n     : total number of real values stored (out,integer)
!   v     : packed potential (inout,real(*))
! !DESCRIPTION:
!   Packs/unpacks the muffin-tin and interstitial parts of the Kohn-Sham
!   potential and magnetic field into/from the single array {\tt v}. This array
!   can then be passed directly to the mixing routine. See routine {\tt rfpack}.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tpack
integer, intent(out) :: n
real(8), intent(inout) :: v(*)
! local variables
integer idm,k
n=0
! pack the Kohn-Sham potential and magnetic field
call rfpack(tpack,n,npmt,npmtmax,vsmt,vsir,v)
do idm=1,ndmag
  call rfpack(tpack,n,npcmt,npcmtmax,bsmt(:,:,idm),bsir(:,idm),v)
end do
! pack the DFT+U potential if required
if (tvmatmt) then
  k=2*lmmaxdm*nspinor*lmmaxdm*nspinor*natmtot
  if (tpack) then
    call dcopy(k,vmatmt,1,v(n+1),1)
  else
    call dcopy(k,v(n+1),1,vmatmt,1)
  end if
  n=n+k
end if
end subroutine
!EOC

