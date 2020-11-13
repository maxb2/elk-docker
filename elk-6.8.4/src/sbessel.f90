
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: sbessel
! !INTERFACE:
subroutine sbessel(lmax,x,jl)
! !INPUT/OUTPUT PARAMETERS:
!   lmax : maximum order of Bessel function (in,integer)
!   x    : real argument (in,real)
!   jl   : array of returned values (out,real(0:lmax))
! !DESCRIPTION:
!   Computes the spherical Bessel functions of the first kind, $j_l(x)$, for
!   argument $x$ and $l=0,1,\ldots,l_{\rm max}$. The recursion relation
!   $$ j_{l+1}(x)=\frac{2l+1}{x}j_l(x)-j_{l-1}(x) $$
!   is used either downwards for $x<l$ or upwards for $x\ge l$. For $x\ll 1$
!   the asymtotic form is used
!   $$ j_l(x)\approx\frac{x^l}{(2l+1)!!}. $$
!   This procedure is numerically stable and accurate to near machine precision
!   for $l\le 50$.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!   Modified to return an array of values, October 2004 (JKD)
!   Improved stability, August 2006 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: lmax
real(8), intent(in) :: x
real(8), intent(out) :: jl(0:lmax)
! local variables
integer l,lst
! rescale limit
real(8), parameter :: rsc=1.d150,rsci=1.d0/rsc
real(8) xi,j0,j1,jt,t1
if ((lmax.lt.0).or.(lmax.gt.50)) then
  write(*,*)
  write(*,'("Error(sbessel): lmax out of range : ",I8)') lmax
  write(*,*)
  stop
end if
if ((x.lt.0.d0).or.(x.gt.1.d5)) then
  write(*,*)
  write(*,'("Error(sbessel): x out of range : ",G18.10)') x
  write(*,*)
  stop
end if
! treat x << 1
if (x.lt.1.d-8) then
  jl(0)=1.d0
  t1=1.d0
  do l=1,lmax
    t1=t1*x/(2*l+1)
    jl(l)=t1
  end do
  return
end if
if (lmax.eq.0) then
  jl(0)=sin(x)/x
  return
end if
xi=1.d0/x
if (x.lt.lmax) then
! for x < lmax recurse down
  j1=1.d0
  j0=0.d0
! starting value for l above lmax
  lst=lmax+lmax/8+14
  do l=lst,lmax+1,-1
    jt=(2*l+1)*j1*xi-j0
    j0=j1
    j1=jt
! check for overflow
    if (abs(j1).gt.rsc) then
! rescale
      jt=jt*rsci
      j0=j0*rsci
      j1=j1*rsci
    end if
  end do
  do l=lmax,0,-1
    jt=(2*l+1)*j1*xi-j0
    j0=j1
    j1=jt
! check for overflow
    if (abs(j1).gt.rsc) then
! rescale
      jt=jt*rsci
      j0=j0*rsci
      j1=j1*rsci
      jl(l+1:lmax)=jl(l+1:lmax)*rsci
    end if
    jl(l)=j0
  end do
! rescaling constant
  t1=1.d0/((j0-x*jl(1))*cos(x)+x*j0*sin(x))
  jl(:)=t1*jl(:)
else
! for x >= lmax recurse up
  jl(0)=sin(x)*xi
  if (lmax.eq.0) return
  jl(1)=(jl(0)-cos(x))*xi
  if (lmax.eq.1) return
  j0=jl(0)
  j1=jl(1)
  do l=2,lmax
    jt=(2*l-1)*j1*xi-j0
    j0=j1
    j1=jt
    jl(l)=j1
  end do
end if
end subroutine
!EOC

