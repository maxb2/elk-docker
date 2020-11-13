
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module moddelf

contains

subroutine delfiles(evec,devec,eval,occ,pmat,epsi)
use modphonon
implicit none
! arguments
logical, optional, intent(in) :: evec,devec,eval,occ,pmat,epsi
! local variables
integer ios
character(256) fext
if (present(evec)) then
! delete the first-variational eigenvector file
  open(122,file=trim(scrpath)//'EVECFV'//trim(filext),iostat=ios)
  close(122,status='DELETE',iostat=ios)
! delete the second-variational eigenvector file
  open(126,file=trim(scrpath)//'EVECSV'//trim(filext),iostat=ios)
  close(126,status='DELETE',iostat=ios)
end if
if (present(devec)) then
! construct the dynamical matrix file extension
  call dynfext(iqph,isph,iaph,ipph,fext)
! delete the eigenvector derivative files
  open(222,file=trim(scrpath)//'DEVECFV'//trim(fext),iostat=ios)
  close(222,status='DELETE',iostat=ios)
  open(226,file=trim(scrpath)//'DEVECSV'//trim(fext),iostat=ios)
  close(226,status='DELETE',iostat=ios)
end if
if (present(eval)) then
! delete the first-variational eigenvalue file
  open(120,file='EVALFV'//trim(filext),iostat=ios)
  close(120,status='DELETE',iostat=ios)
! delete the second-variational eigenvalue file
  open(124,file='EVALSV'//trim(filext),iostat=ios)
  close(124,status='DELETE',iostat=ios)
end if
if (present(occ)) then
! delete the occupation number file
  open(130,file='OCCSV'//trim(filext),iostat=ios)
  close(130,status='DELETE',iostat=ios)
end if
if (present(pmat)) then
! delete the momentum matrix elements file
  open(150,file='PMAT.OUT',iostat=ios)
  close(150,status='DELETE',iostat=ios)
end if
if (present(epsi)) then
! delete the inverse epsilon file
  open(180,file='EPSINV.OUT',iostat=ios)
  close(180,status='DELETE',iostat=ios)
end if
end subroutine

end module

