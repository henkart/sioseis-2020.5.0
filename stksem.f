      SUBROUTINE STKSEM(INWIN,STKWIN,IAPSCR,ENERGY,NSAMPS,ONE,DIVIS)
C    STKSEM IS A SPECIAL PURPOSE SUBROUTINE USED BY THE SEMBLENCE VELOCITY
C  ANALYSIS.  THIS IS AN AP VECTOR FUNCTION CHAINER ROUTINE THAT MUST BE CALLED
C  FOR EACH VELOCITY WINDOW.  IT STACKS THE INPUT WINDOW INTO A STACKED TRACE
C  BUFFER IN THE AP, AND COMPUTES THE ENERGY OF THE INPUT WINDOW AND SUMS IT TO
C  THE MEAN ENERGY ARRAY. STKSEM ALSO BUILDS AN ARRAY FOR THE DIVISOR.
C    THIS ROUTINE COULD BE COMBINED WITH CVNMO AND LOOPED BY THE VFC FOR EACH
C  EXCEPT THAT THE AP ADDRESS OF THE INPUT WINDOW IS NOT A LINEAR INCREMENT
C  FROM WINDOW TO WINDOW.  IF THERE WAS A WAY TO DO INDIRECT ADDRESSING THIS
C  WOULD SOLVE THE PROBLEM. (CVNMO COMPUTES THE INDEXES OF THE WINDOW CENTER
C  POINTS.  THESE INDEXES ARE NEEDED AS AP ADDRESSES BY THIS ROUTINE.)
C
C  ARGUMENTS:
C  INWIN  - THE AP ADDRESS OF THE INPUT WINDOW.
C  STKWIN - THE AP ADDRESS OF THE STACK WINDOW (WHERE THE INPUT WILL BE ADDED).
C  IAPSCR - THE AP ADDRESS OF A SCRATCH WORD THAT CAN BE USED.
C  ENERGY - THE AP ADDRESS OF THE SCALAR ENERGY OF THE INPUT WINDOW.
C  NSAMPS - THE NUMBER OF SAMPLES IN THE INPUT WINDOW.
C  ONE    - THE AP ADDRESS OF THE NUMBER ONE (1.).
C
C
C  WRITTEN AND COPYRIGHTED BY:
C  PAUL HENKART, SCRIPPS INSTITUTION OF OCEANOGRAPHY, 29 APRIL 1981
C  ALL RIGHTS ARE RESERVED BY THE AUTHOR.  PERMISSION TO COPY OR REPRODUCE THIS
C  SUBROUTINE, BY COMPUTER OR OTHER MEANS, MAY BE OBTAINED ONLY FROM THE AUTHOR.
C
      IMPLICIT INTEGER (A-Z)
      CALL VADD(INWIN,1,STKWIN,1,STKWIN,1,NSAMPS)
      CALL SVESQ(INWIN,1,IAPSCR,NSAMPS)                                  /* COMPUTE THE ENERGY
      CALL VSADD(ENERGY,1,IAPSCR,ENERGY,1,1)
      CALL VSADD(DIVIS,1,ONE,DIVIS,1,1)
      RETURN
      END
