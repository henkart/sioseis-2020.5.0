      SUBROUTINE SEMST(ISTK,ISEM,ANS,IDIV,N,N2,NDOWS,char,IFMT)
!      TRACE SEMBL,TEMP,K
!    SEMST IS THE VERSION OF SEMSTK THAT RUNS IN HOST MEMORY RATHER THAN IN THE AP.
!  SEMST PERFORMS THE FOLLOWING TASKS: 1) COMPUTES THE ENERGY OF THE STACKED
! TRACE (THE ENERGY WINDOWS ARE N TRACE SAMPLES LONG AND ARE COMPUTED EVERY
! N2 TRACE SAMPLES).  2) DIVIDES THE STACKED TRACE ENERGY BY THE NUMBER OF
! WINDOWS THAT WENT INTO THE STACKED TRACE, FOR EACH WINDOW.  THIS IS THE MEAN
! ENERGY OF THE STACKED TRACE.  3) DIVIDES THE STACKED TRACE MEAN ENERGY BY THE
! MEAN ENERGY OF THE INPUT WINDOWS.  4) MULTIPLIES THIS SEMBLANCE BY char (9).
!
! ARGUMENTS:
! ISTK - THE APMEM INDEX OF THE STACKED TRACE.
! ISEM - THE APMEM INDEX OF THE MEAN ENERGY OF THE INPUT TRACES.
! ANS  - THE ARRAY TO RECEIVE THE REAL SEMBLANCE VALUES.
! IDIV - THE APMEM INDEX OF THE NUMBER OF WINDOWS FORMING THE STACK.
! N    - THE NUMBER OF SAMPLES IN EACH WINDOW.
! N2   - THE DISTANCE BETWEEN WINDOWS (USUALLY N/2).
! NDOWS - THE NUMBER OF WINDOWS.
! IFMT - THE OUTPUT FORMAT - BETTER BE 1 - INCLUDED FOR COMPATIBILITY WITH SEMSTK
!
! COPYRIGHTED BY PAUL HENKART, 23 JANUARY 1983
!
      COMMON /APMEM/ A(32766)
      DIMENSION ans(1)
!
      IF(IFMT.EQ.1) GO TO 2
      PRINT 1
1     FORMAT(' ***  ERROR  ***  SEMST CAN NOT DO A FORMAT ',
     *  'OTHER THAN 1.')
      STOP
2     CONTINUE
      IDIV1=IDIV
      ISEM1=ISEM
      K=ISTK-1
      DO 20 I=1,NDOWS
         ENERGY=0.
         IF(A(IDIV1).LE.0.) A(IDIV1)=1.
         DO 10 J=1,N                                                    ! DO A WINDOW LENGTH
            IF(A(K+J).EQ.0.) GO TO 10                                   ! NO SENSE IN PLAYING WITH ZEROES!
            TEMP=A(K+J)/A(IDIV1)                                        !  THE STACKED TRACE HAS NOT BEEN DIVIDED YET!!!
            ENERGY=ENERGY+TEMP*TEMP                                     !  SUM OF THE SQUARES!
10       CONTINUE
!         PRINT 111,I,IDIV1,ISEM1,A(IDIV1),A(ISEM1),ENERGY
!  111    FORMAT(1X,3I10,3F12.5)
         IF(A(ISEM1).EQ.0.) A(ISEM1)=1.
         SEMBL=(ENERGY/(A(ISEM1)/A(IDIV1))) * char
         ANS(I)=SEMBL
!         IDIV1=IDIV1+1
         ISEM1=ISEM1+1
         K=K+N                                                          ! FIND THE BEGINNING OF THE NEXT WINDOW
   20 CONTINUE
      RETURN
      END
