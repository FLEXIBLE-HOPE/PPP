!$RDJASON
      SUBROUTINE ucl_grid_reader(DATA,FILNAM)
!********1*********2*********3*********4*********5*********6*********7**
! RDJASON   
!                                                                       
! FUNCTION         read in the radiation pressure model grid file 
!                  UCL MODEL 
!                                                                       
! I/O PARAMETERS:                                                       
!                                                                       
!   NAME    I/O  A/S   DESCRIPTION OF PARAMETERS                        
!   ------  ---  ---   ------------------------------------------------                 
!   DATA     O    A    JASON SOLAR RADIATION GRID VALUE ARRAY
!   IUNT     I    S    UNIT NUMBER OF JASON GRID FILE
!   FILNAM   I    S    FILENAME OF GRID FILE
!                          GRIDX2.txt
!                          GRIDY2.txt
!                          GRIDZ2.txt
! COMMENTS                                                              
!                                                                       
!********1*********2*********3*********4*********5*********6*********7**
      IMPLICIT REAL (A-H,O-Z), LOGICAL(L)                                                         
      CHARACTER*80 FILNAM
      CHARACTER*80 store
      INTEGER minlon, maxlon, minlat, maxlat
      !real, dimension(181,361) :: data
      double precision, dimension(181,361) :: data
      integer :: iunt
                                                                     
!********************************************************************** 
! START OF EXECUTABLE CODE                                              
!********************************************************************** 
			
      INQUIRE(EXIST=LEXIST,FILE=FILNAM)
      OPEN(unit=iunt,file=filnam,status='old',iostat=ierror)
      IF(ierror.ne.0)THEN
          WRITE(*,150) filnam,ierror
  150     FORMAT(' ','Error opening file: ',A,'IOSTAT = ',I6)
          error=-1
      ELSE
          READ(iunt,*)
          REWIND iunt

!      read in file header
           READ(iunt,*)store
           READ(iunt,*)maxcolumn,maxrow
           COLMAX=DFLOAT(maxcolumn)
           ROWMAX=DFLOAT(maxrow)
           READ(iunt,*)minlon,maxlon
           RMNLON=DFLOAT(minlon)
           RMXLON=DFLOAT(maxlon)
           READ(iunt,*)minlat,maxlat
           RMNLAT=DFLOAT(minlat)
           RMXLAT=DFLOAT(maxlat)
           READ(iunt,*)dummy, dummy
   ! write(6,*) ' dbg reading file ',maxcolumn,maxrow,minlat,minlon
   ! write(6,*) ' dbg dummy ',dummy

           DO 11 i=1,maxrow
              READ(iunt,*)(data(i,j),j=1,maxcolumn)              
              icount = icount + maxcolumn
   11     CONTINUE

      CLOSE(iunt)

      ENDIF


      RETURN
      END                                           
