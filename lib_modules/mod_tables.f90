!*
MODULE tables
!*
USE par
USE const
USE ISO_FORTRAN_ENV
IMPLICIT NONE

  TYPE T_FileTable
    INTEGER(IT) :: nflnm=0
    CHARACTER(LEN_FILENAME_SYB)  :: csysnm(50)
    CHARACTER(LEN_FILENAME) :: cflnm (50)
  END TYPE T_FileTable

  TYPE(T_FileTable), SAVE :: FT_GLOBAL

  CONTAINS
  
  ! Read FileTable
  SUBROUTINE read_filetable(cflnm,FT)
    
  USE par
  USE const
  IMPLICIT NONE

  ! Input parameters:
  !-----------------------------------------------------------------------------------
  CHARACTER(LEN=*) :: cflnm
  TYPE(T_FileTable) :: FT
    
    ! Local parameters:
    !-----------------------------------------------------------------------------------
    INTEGER(IT),PARAMETER :: iflun=1001
    CHARACTER(LEN_STRING) cword(10),cline_read
    INTEGER(IT) nword,ierr
    
    OPEN(UNIT=iflun,FILE=cflnm,STATUS='old',IOSTAT=ierr)
    IF (ierr .NE. 0) THEN
      WRITE(ERROR_UNIT,'(A)') '***ERROR(read_filetable): open file ',trim(cflnm)
      CALL EXIT(1)
    ENDIF
    
    DO_RDFL: DO
      READ(iflun,'(A)',END=100,iostat=ierr) cline_read
      IF(ierr.ne.0) EXIT
      
      IF(cline_read(1:1).ne.' ') THEN
        CYCLE
      ENDIF

      CALL split_str(.true.,cline_read,' ',' ',' ',nword,cword)
      
      IF (nword .LT. 2) THEN
        WRITE(OUTPUT_UNIT,'(A)')'%%%WARNING(read_filetable): line is invalid',TRIM(cline_read)
        CYCLE
      ENDIF
      
      FT%nflnm=FT%nflnm+1
      
      FT%csysnm(FT%nflnm)=cword(1)
      FT%cflnm(FT%nflnm)=cword(2)

    ENDDO DO_RDFL 
    
    100 CONTINUE

    FT_GLOBAL = FT
    CLOSE(iflun)
  
    RETURN
  
  END SUBROUTINE read_filetable

  
  FUNCTION get_flname(FT,csysnm,cflnm)
  USE par
  IMPLICIT NONE
  
  ! Input parameters:
  !-----------------------------------------------------------------------------------
  TYPE(T_FileTable) ::  FT
  CHARACTER(LEN=*) :: csysnm
  CHARACTER(LEN=*) :: cflnm
    
    ! Local parameters:
    !-----------------------------------------------------------------------------------
    LOGICAL(LG) :: get_flname
    INTEGER(IT) :: i
    
    get_flname=.FALSE.
    DO i=1,FT%nflnm
      IF (TRIM(csysnm) .EQ. TRIM(FT%csysnm(i))) THEN
        cflnm=TRIM(FT%cflnm(i))
        get_flname=.TRUE.
        EXIT
      ENDIF
    ENDDO
    
    RETURN
  
  END FUNCTION get_flname


  !*
  FUNCTION f_tableFileName(fileType)
  !*
  USE par
  IMPLICIT NONE
  
  !*
  ! Declare
  !!-------
  CHARACTER(LEN_FILENAME) :: f_tableFileName
  CHARACTER(LEN=*) :: fileType

    !*
    ! Local variable declare
    !!----------------------
    LOGICAL(LG) :: ifTypeExist

    !*
    ! Start of executable code
    !!------------------------
    
    ifTypeExist = get_flname(FT_GLOBAL,fileType,f_tableFileName)
    

    RETURN
  END FUNCTION f_tableFileName
 
  SUBROUTINE split_str(lnoempty,string,c_start,c_end,seperator,nword,word)
  USE par
  USE const
  USE ISO_FORTRAN_ENV
  IMPLICIT NONE

  !*
  ! The arguments
  !!--------------------------
  LOGICAL(LG) :: lnoempty
  CHARACTER(LEN=*) :: string, word(1:*)
  CHARACTER :: c_start,c_end,seperator
  INTEGER(IT) :: nword

    !*
    ! The local variables
    !!--------------------------
    INTEGER(IT) :: i0,i1,ilast,i
    CHARACTER(LEN_STRING) :: line

    !*
    ! Start the exectuable code
    !!--------------------------

    line=string
    IF (LEN(string) .GT. LEN_STRING) THEN
      WRITE(ERROR_UNIT,'(A,I5)') '***ERROR(split_line): input string length > ', LEN_STRING
      CALL exit(1)
    END IF

    i0=0
    i1=LEN_TRIM(line)+1
    IF (LEN_TRIM(c_start) .NE. 0) i0=INDEX(line,c_start)
    IF (LEN_TRIM(c_end  ) .NE. 0) i1=INDEX(line,c_end  )

    nword=0
    IF (i1.EQ.0 .OR. i0.GT.i1) RETURN

    ilast=i0+1
    line(i1:i1)=seperator
    DO i=i0+1,i1
      IF (line(i:i) .EQ. seperator) THEN
        IF (LEN_TRIM(line(ilast:i-1)) .NE. 0) THEN
          nword=nword+1
          word(nword)=line(ilast:i-1)
        ELSE
          IF (.NOT.lnoempty) THEN
            nword=nword+1
            word(nword)=' '
          END IF
        END IF
        ilast=i+1
      END IF
    END DO

    RETURN

  END SUBROUTINE


END MODULE tables
  
