module read_initFile

#include "definition.h"
  
  implicit none
  integer, parameter :: fileLen=100

contains

subroutine read_initFileReal(fileName,varName,varValue)

  implicit none
  character(len=*),intent(IN) :: fileName,varName
  real, intent(OUT) :: varValue

  integer :: i,openStatus,inputStatus
  real :: simInitVars
  character(len=MAX_STRING_LENGTH) :: simCharVars
  integer :: pos1,pos2

  open(unit = 10, file=fileName, status='old',IOSTAT=openStatus,FORM='FORMATTED',ACTION='READ')

  do i=1,fileLen
     read(10, FMT = 100, IOSTAT=inputStatus) simCharVars
     pos1 = index(simCharVars,varName)
     pos2 = pos1+len_trim(varName)
     if (pos2 > len_trim(varName)) then
        read(simCharVars(pos2+1:),*)simInitVars
        varValue = simInitVars
     endif
  end do

  close(10)

100 FORMAT(A, 1X, F3.1)

end subroutine read_initFileReal



subroutine read_initFileInt(fileName,varName,varValue)

  implicit none
  character(len=*),intent(IN) :: fileName,varName
  integer, intent(OUT) :: varValue

  integer :: i,openStatus,inputStatus
  integer :: simInitVars
  character(len=MAX_STRING_LENGTH) :: simCharVars
  integer :: pos1,pos2

  open(unit = 11, file=fileName, status='old',IOSTAT=openStatus,FORM='FORMATTED',ACTION='READ')
  
  do i=1,fileLen
     read(11, FMT = 101, IOSTAT=inputStatus) simCharVars
     pos1 = index(simCharVars,varName)
     pos2 = pos1+len_trim(varName)
     if (pos2 > len_trim(varName)) then
        read(simCharVars(pos2+1:),*)simInitVars
        varValue = simInitVars
     endif
  end do
  
  close(11)

101 FORMAT(A, 1X, I5)

end subroutine read_initFileInt


subroutine read_initFileBool(fileName,varName,varValue)

  implicit none
  character(len=*),intent(IN) :: fileName,varName
  logical, intent(OUT) :: varValue

  integer :: i,openStatus,inputStatus
  logical :: simInitVars
  character(len=MAX_STRING_LENGTH) :: simCharVars
  integer :: pos1,pos2

  open(unit = 12, file=fileName, status='old',IOSTAT=openStatus,FORM='FORMATTED',ACTION='READ')
  
  do i=1,fileLen
     read(12, FMT = 102, IOSTAT=inputStatus) simCharVars
     pos1 = index(simCharVars,varName)
     pos2 = pos1+len_trim(varName)

     if (pos2 > len_trim(varName)) then
        read(simCharVars(pos2+1:),*)simInitVars
        varValue = simInitVars
     endif
  end do

  close(12)

102 FORMAT(A, 1X, L)

end subroutine read_initFileBool




subroutine read_initFileChar(fileName,varName,varValue)

  implicit none
  character(len=*),intent(IN)  :: fileName,varName
  character(len=*),intent(OUT) :: varValue

  integer :: i,openStatus,inputStatus
  character(len=MAX_STRING_LENGTH) :: simInitVars
  character(len=MAX_STRING_LENGTH) :: simCharVars
  integer :: pos1,pos2

  open(unit = 13, file=fileName, status='old',IOSTAT=openStatus,FORM='FORMATTED',ACTION='READ')
  
  do i=1,fileLen
     read(13, FMT = 103, IOSTAT=inputStatus) simCharVars
     pos1 = index(simCharVars,varName)
     pos2 = pos1+len_trim(varName)

     if (pos2 > len_trim(varName)) then
        read(simCharVars(pos2+1:),*)simInitVars
        varValue = simInitVars
     endif
  end do

  close(13)

103 FORMAT(A, 1X, A)

end subroutine read_initFileChar

end module read_initFile
