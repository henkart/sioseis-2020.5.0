      SUBROUTINE RLSEAP(BUFOUT,LSAMPS)
!     RLSEAP WAITS FOR AN AP PROGRAM TO FINISH RUNNING THEN TRANSFERS THE DATA
!  FROM AP LOCATION IN TO PRIME MEMORY LOCATION BUFOUT.  OF COURSE RLSEAP WILL
!  DO NOTHING IF THE DATA IS NOT IN THE AP!.  THE AP STATUS IS CHECKED AND A
!  PRINTER MESSAGE IS DISPAYED IF THERE WAS AN AP OVERFLOW, OR UNDERFLOW, OR
!  DIVIDE BY ZERO FAULT.
!     IF THERE ISN'T AN AP ON THE SYSTEM, THE DATA IS IN COMMON /APMEM/.
!  THIS MEMORY JUST LOOKS LIKE A REAL AP DATA MEMORY SO THAT ALL THE INDEXES
!  AND DATA ADDRESSES CAN REMAIN THE SAME, REGARDLESS OF AP EXISTANCE.
!
!  ARGUMENTS:
!    BUF    - THE ARRAY IN THE PROGRAM TO RECEIVE THE ARRAY IN THE AP IF THE
!             DATA IS SUPPOSED TO COME OUT OF THE AP
!    NSAMPS - THE NUMBER OF SAMPLES IN THE ARRAY BUF.
!
!
!    PAUL HENKART, SCRIPPS OCEANOGRAPHY, NOVEMBER 1979
!
      COMMON /EDITS/ IERROR,IWARN,IRUN,NOW,ICOMPT,isite, maxsamps,
     & nbperw, ireal
      COMMON /APMEM/ A(32766)
      COMMON /SIOAP/ IASGND,IRELSE,IN,IOUT,NEXTAD,LAPSIZ,IFREE,IUSEAP,
     *  IDECIM
!          IASGND - AP ASSIGNMENT SWITCH
!                 =0,  AP IS NOT ASSIGNED - CALL APINIT
!                 =1,  AP IS ASSIGNED - DO NOT CALL APINIT
!          IRELSE - AP RELEASE SWITCH
!                 =0,  AP SHOULD BE RELEASED BEFORE RETURNING
!                 =1,  AP SHOULD BE LEFT ASSIGNED
!          IN     - A SWITCH INDICATING WHETHER THE DATA IS ALREDY IN THE AP
!                 =0,  THE DATA IS NOT IN THE AP.
!                 >0,  THE DATA IS IN THE AP AT ADDRESS IN
!          IOUT   - A SWITCH INDICATING WHETHER THE DATA SHOULD BE LEFT IN THE AP
!                 =0,  THE DATA SHOULD BE MOVED BACK TO THE HOST
!                 =1,  THE DATA CAN BE LEFT IN THE AP WITHOUT AN APGET
!          NEXTAD - THE NEXTA AVAILABLE FREE LOCATION IN THE AP
!          LAPSIZ - THE SIZE OF THE AP DATA MEMORY
!          IFREE  - A SIGNAL INDICATING THAT SOME SUBROUTINE HAS SOMETHING IN
!                   THE AP THAT MUST BE SAVED, THUS PREVENTING THE AP FROM BEING
!                   RELEASED.
!                 =0, NO SUBROUTINE HAS ANYTHING TO SAVE
!                 <>0, SOMETHINGS ARE IN THE AP THAT MUST BE SAVED - DO NOT
!                      RELEASE THE AP.
!          IUSEAP - A SIGNAL INDICATING WHETHER TO USE THE AP OR NOT
!                 =0, DON'T USE THE AP (THERE MIGHT NOT BE ONE!)
!                 =1, USE THE AP
!****
!****
      DIMENSION ISTAT(4),BUFOUT(1)
      PARAMETER (IFMT=2)                                                ! THE DATA FORMAT FOR APGET - 2=HOST FLOATING POINT
      INTEGER*4 IASGND,IRELSE,IN,IOUT,NEXTAD,LAPSIZ,IFREE,IUSEAP,IDECIM
      INTEGER*4 LSAMPS

!      PRINT *,' call rlseap, iout=',iout,' in=',in,' iasgnd=',iasgnd
      IF(IOUT.NE.0) RETURN                                              ! DON'T RELEASE THE AP IF DATA IS LEFT IN IT!!
      IF(IN.EQ.0) RETURN                                                !  DON'T DO ANYTHING IF THE DATA IS NOT IN THE AP!!!
      NSAMPS=LSAMPS                                                     ! CONVERT TO 16BIT INTEGER
      IF(IUSEAP.EQ.0) GO TO 100                                         ! IS THERE AN AP?
      CALL APWR 
      JIN=IN                                                            ! CONVERT TO 16 BIT INTEGER
      CALL APGET(BUFOUT,JIN,NSAMPS,IFMT)                                ! GET THE SCALED TRACE BACK
      IN=0                                                              ! SET IT TO DATA NOT IN AP! - ASSUME IT IS MODIFIED IF IT IS TAKEN OUT
      CALL APSTAT(IERR,ISTAT)                                           ! WAITS FOR AP AND RETURN AP STATUS
      IF(IERR.NE.0) PRINT 10,ISTAT
   10 FORMAT(' AP STATUS: OVERFLOW=',I2,' UNDERFLOW=',I2,' DIVIDE BY ',
     *    'ZERO=',I2,' FORMAT=',I2)
      IF(IRELSE.EQ.1) RETURN
      IF(IFREE.NE.0) RETURN
      CALL APRLSE                                                       !  RELEASE THE AP BACK TO THE SYSTEM
      NEXTAD=0                                                          ! INAP ONLY ASSIGNS NEXTAD IF IT IS 0!
   20 IASGND=0
      RETURN
!****
!****     NO AP
!****
  100 J=IN-1
      IF( nsamps .GT. maxsamps*3 ) THEN
          PRINT *,' ***  ERROR  ***  Call to rlseap for ',nsamps,
     &     ' will overwrite maxsam*3 (',maxsamps*3,') samples.'
          STOP
      ENDIF
      DO I=1,NSAMPS
  110    BUFOUT(I)=A(J+I)                                                  ! MOVE THE DATA OUT OF THE AP SIMULATOR MEMORY
      ENDDO
      IN=0
      IF(IRELSE.EQ.1)  RETURN
      IF(IFREE.NE.0) RETURN
      GO TO 20
      END
