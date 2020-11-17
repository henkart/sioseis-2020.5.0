      SUBROUTINE FINDV(T,VTP,N,V)
!     FINDV FINDS AND RETURNS A VELOCITY V CORRESPONDING TO A GIVEN T FROM A
!  GIVEN VELOCITY FUNCTION VTP (VELOCITY-TIME-PAIRS) N LONG.
!     THE FIRST VELOCITY OF VTP IS USED FOR ALL T'S BEFORE THE FIRST T OF VTP.
!  LIKEWISE, THE LAST VELOCITY OF VTP IS USED FOR ALL T'S AFTER THE LAAST
!  GIVEN IN THE VTP ARRAY.
!
!   COPYRIGHTED BY:
!  PAUL HENKART, SCRIPPS INSTITUTION OF OCEANOGRAPHY, APRIL 1980
!
c  mod 1 Apr 10 - remove arithmetic IF statements
c
      DIMENSION VTP(1111)
!
      DO 100 I=1,N,2
      VX=VTP(I)
      TX=VTP(I+1)
c      IF(T-TX)100,200,10
      IF( t - tx .LT. 0 ) GOTO 100
      IF( t - tx .EQ. 0 ) GOTO 200
   10 IF(I.EQ.1) GO TO 200
      IF(VTP(I-1).EQ.VTP(I+1)) GO TO 300                         ! WATCH OUT FOR DIVIDE BY 0 (CONSTANT V)
      VX=(T-VTP(I-1))/(VTP(I+1)-VTP(I-1))*(VTP(I)-VTP(I-2))+VTP(I-2)
      GO TO 200
  100 CONTINUE
  200 V=VX
      RETURN
  300 V=VTP(I-1)
      RETURN
      END