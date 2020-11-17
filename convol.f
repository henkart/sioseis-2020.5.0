      SUBROUTINE CONVOL(ITYPE,BUFIN,NSAMPS,FILTER,NFPTS,OBUF,NOUT)
!     CONVOL PERFORMS TIME DOMAIN CONVOLUTION AND CORRELATION BETWEEN ARRAYS
! BUFIN AND FILTER.  THE OUTPUT IS RETURNED IN OBUF AND IS  VARIABLE IN LENGTH
! DEPENDING ON THE TYPE OF OPERATION REQUESTED.
!
! ARGUMENTS:
!        ITYPE  - TYPE OF CONVOLUTION OR CORRRELATION DESIGNATOR.
!               =-2,  PARTIAL CONVOLUTION.  THIS TYPE IS DESIGNED FOR ZERO PHASE
!                     FILTERING WHERE NEGATIVE TIME AFTER FILTERING IS NOT WANTED
!                     THIS ASSUMES THAT THE FIRST TIME OF THE DATA IS 0 AND THE
!                     FILTER IS SYMMETRICAL.
!               =-1,  FULL CONVOLUTION.  THE TIME REFERENCE OF THE FIRST OUTPUT
!                     POINT IS T1+T2, WHERE T1 IS THE TIME OF THE FIRST SAMPLE
!                     IN BUFIN AND T2 IS THE TIME OF THE FIRST SAMPLE IN ARRAY
!                     FILTER.  E.G.  IF A ZERO PHASE FILTER HAS TIME POINTS -19
!                     TO +19 AND THE FIRST SAMPLE IN BUFIN IS AT TIME 0, THEN
!                     THE TIME OF THE FIRST OUTPUT SAMPLE IS -19.
!               =+1,  FULL CORRELATION.  BOTH POSITIVE AND NEGATIVE LAGS ARE
!                     COMPUTED.  CORRELATIONS ARE EVEN FUNCTIONS, THEREFORE HALF
!                     CORRELATIONS GIVE THE SAME INFORMATION AND SHOULD BE USED
!                     MOST OF THE TIME.
!               =+2,  HALF CORRELATION.  ONLY THE RIGHT HALF OF THE CORRELATION
!                     IS PERFORMED.  THIS TYPE IS DESIGNED FOR VIBROSEIS
!                     DESWEEPING WHERE A FULL CORRELATION WOULD BE A TOTAL
!                     WASTE OF TIME!!!
!        BUFIN  - INPUT ARRAY TO BE FILTERED (CONVOLVED) OR CORRELATED.
!        NSAMPS - THE NUMBER OF SAMPLES IN BUFIN TO FILTER.
!        FILTER - AN ARRAY OF REAL FILTER WEIGHTS.  AN ARRAY OF VALUES TO BE
!                 CORRELATED.
!        NFPTS  - THE NUMBER OF FILTER POINTS IN THE FILTER ARRAY.
!        OBUF   - THE OUTPUT ARRAY.
!        NOUT   - THE NUMBER OF OUTPUT SAMPLES TO GENERATE.  FOR CONVOLUTION
!                 THIS SHOULD NORMALLY BE MSAMPS+MPTS.  FOR CORRELATION THIS IS
!                 THE NUMBER OF LAGS TO OUTPUT, IT MIGHT BE  MIN0(MSAMPS,MPTS).
!
!   NOTE:  LEADING AND TRAILING ZEROES DO NOT HAVE TO BE GIVEN.
!
!   COPYRIGHTED BY:
!       PAUL HENKART, SCRIPPS INSTITUTION OF OCEANOGRAPHY, FEBRUARY 1979
!       REWRITTEN IN NOVEMBER 1979 FOR THE FPS AP.
!
      DIMENSION BUFIN(1),FILTER(1),OBUF(1)
      COMMON /SIOAP/ IASGND,IRELSE,IN,IOUT,NEXTAD,LAPSIZ,IFREE
!           IASGND - AP ASSIGNMENT SWITCH
!                  =0,  AP IS NOT ASSIGNED - CALL APINIT
!                  =1,  AP IS ASSIGNED - DO NOT CALL APINIT
!          IRELSE - AP RELEASE SWITCH
!                 =0,  AP SHOULD BE RELEASED BEFORE RETURNING
!                 =1,  AP SHOULD BE LEFT ASSIGNED
!          IN     - A SWITCH INDICATING WHETHER THE DATA IS ALREDY IN THE AP
!                 =0,  THE DATA IS NOT IN THE AP.
!                 =1,  THE DATA IS IN THE AP AT ADDRESS IN
!          IOUT   - A SWITCH INDICATING WHETHER THE DATA SHOULD BE LEFT IN THE AP
!                 =0,  THE DATA SHOULD BE MOVED BACK TO THE HOST
!                 =1,  THE DATA CAN BE LEFT IN THE AP WITHOUT AN APGET
!****
!****
!****
      IPAD1=IN                                                         ! SET UP THE ADDRESS OF THE INPUT DATA IN THE AP
      IF(IN.EQ.0) IPAD1=1
      ITEMP=ITYPE+3
      GO TO (110,120,120,120,130), ITEMP
  110 NZEROS=NFPTS/2                                                   !  PARTIAL CONVOLUTION
  111 IDATA=NZEROS+IPAD1                                               !  THE ADDRESS OF THE DATA IN THE AP
      IPAD2=IDATA+NSAMPS                                               !  THE ADDRESS OF THE BACK END PAD
      IAPFLT=IPAD2+NZEROS                                              !  THE ADDRESS OF THE FILTER IN THE AP
      CALL VMOV(IN+NSAMPS-1,-1,IDATA+NSAMPS-1,-1,nsamps)               ! MOVE THE DATA FOR  PAD
      CALL VCLR(IPAD1,1,NZEROS)                                        ! CREATE THE PAD
      GO TO 180
  120 NZEROS=NFPTS-1                                                   !  FULL CONVOLUTION AND CORRELATION
      GO TO 111
  130 NZEROS=NFPTS-1                                                   !  HALF CORRELATION
      IDATA=IPAD1                                                      !  DON'T NEED A FRONT END PAD
      IPAD2=IDATA+NSAMPS                                               !  ADDRESS OF THE BACKEND PAD
      IAPFLT=IPAD2+NZEROS                                              !  ADDRESS OF THE FILTER IN THE AP
  180 CALL VCLR(IPAD2,1,NZEROS)                                        ! START THE BACK END PAD
      CALL APPUT(FILTER,IAPFLT,NFPTS,2)                                !  START THE FILTER INTO THE AP
      INC=1                                                            ! THE DATA INCREMENT FOR CORRELATION
      IADDR=IAPFLT                                                     ! THE DATA ADRESS FOR CORRELATION
      IF(ITYPE.GE.0) GO TO 200
      INC=-1                                                           ! THE INC FOR CONVOLUTION
      IADDR=IAPFLT+NFPTS-1
  200 CONTINUE
      CALL APWAIT                                                      ! WAIT FOR THE DATA TRANSFERS AND CLEARING TO FINISH
      CALL CONV(IPAD1,1,IADDR,INC,IPAD1,1,NOUT,NFPTS)
      IN=IPAD1
      RETURN
      END