      SUBROUTINE uadded(BUF,LBUF)
!                                PROCESS UADD
!                                ------- ----
!  
!  Document Date: 15 May 1996
!  
!  PROCESS UADD adds user given points to every specified seismic trace.
!  The user may specify the start time of the addition.
!  
!  A common use of process uadd is to subtract the average trace from all
!  the data traces in GPR (ground penetrating radar).
!  
!  Spatial interpolation is not available.
!
!  Each parameter list must be terminated with the word END.  The entire
!  set of UADD parameters must be terminated by the word END.
!  
!  THE PARAMETER DICTIONARY
!  --- --------- ----------
!  
!  PTS   - A list of points to be added to each seismic trace.  The 
!          values must be separated by a blank, tab or new line.
!          The points may be given on as many lines as necessary.
!          When the values are given in an exponential format,
!          positive exponents must have a + rather than a blank.
!          Required. e.g. pts 1 2 3 4 5 6 7 8 9 10 11
!  
!  FNO   - The first shot (or rp) to be added.  Shot (rp) numbers
!          must increase monotonically.
!          Preset = 0
!  
!  LNO   - The last shot (rp) number to be added.  LNO must be
!          larger than FNO in each list and must increase list to list.
!          Preset = 999999999
!  
!  STIME - The start time, in seconds, of the addition.  i.e. the time of
!          the first point given in PTS.
!          Preset = the delay of the trace . e.g. stime 1
!  
!  END    - Terminates each parameter list.
!  
!  
!  Written and copyrighted (c) by:
!  Paul Henkart, Scripps Institution of Oceanography, May 1996
!  ALL RIGHTS RESERVED.
!  
!  THE PARAMETER LIST ON DISC IS:
!  WORD 1)  FNO
!       2)  LNO
!       3)  NPTS
!       4)  STIME
!       5)  LPRINT
!       6)  NPTS POINTS    ****  THIS IS A VALIABLE FROM LIST TO LIST!!!
!
!
!  ARGUMENTS:
!  BUF    - A SCRATCH ARRAY AT LEAST 500 32 BIT WORDS LONG.
!  LBUF   - THE SAME ARRAY BUT THE 32 BIT INTEGER EQUIVALENT.  NEEDED
!           BECAUSE PRIME FORTRAN DOESN'T ALLOW EQUIVALENCING OF ARGUMENTS.
!
      PARAMETER (NPARS=6)                                               ! THE NUMBER OF USER PARAMETERS
      PARAMETER (MULTIV=6)                                              ! THE PARAMETER NUMBER OF THE FIRST MULTI-VALUED PARAMETER
      COMMON /UADD/ MUNIT,NLISTS
      DIMENSION BUF(111),LBUF(111)
      CHARACTER*6 NAMES(NPARS)
      CHARACTER*1 TYPES(NPARS)
      DIMENSION LENGTH(NPARS)
      CHARACTER*80 TOKEN
      DIMENSION VALS(NPARS),LVALS(NPARS)
      COMMON /EDITS/ IERROR,IWARN,IRUN,NOW,ICOMPT
      COMMON /sioln1/ cbuf
      COMMON /sioln2/ idum1, idum2, iprint
      CHARACTER*100 cbuf
      INTEGER FNO
!
!
      EQUIVALENCE (FNO,LVALS(1)),
     2            (LNO,LVALS(2)),
     3            (LPRINT,LVALS(3)),
     4            (NPTS,LVALS(4)),
     5            (STIME,VALS(5)),
     6            (PTS,VALS(6))
      DATA NAMES/'FNO   ','LNO   ','LPRINT','NPTS  ','STIME ',
     *           'PTS   '/
      DATA TYPES/'L','L','L','L','F','F'/
      DATA LENGTH/3,3,6,4,5,3/
!****
!****      SET THE PRESETS
!****
      FNO=0
      LNO=999999999
      NPTS=0
      STIME = -1.
      LPRINT=0
      IADDWB=0
      LLNO = -1
      NLISTS=0
      NS=0
      isavep = iprint
!****
!****   GET A PARAMETER FILE
!****
      CALL GETFIL(1,MUNIT,TOKEN,ISTAT)
!****
!****   THE CURRENT COMMAND LINE IN THE SYSTEM BUFFER MAY HAVE THE PARAMETERS.
!****   GET A PARAMETER LIST FROM THE USER.
!****
      NTOKES=1
  100 CONTINUE
      CALL GETOKE(TOKEN,NCHARS)                                         ! GET A TOKEN FROM THE USER PARAMETER LINE
      CALL UPCASE(TOKEN,NCHARS)                                         ! CONVERT THE TOKEN TO UPPERCASE
      IF(NCHARS.GT.0) GO TO 150
      IF(NOW.EQ.1) PRINT 140
  140 FORMAT(' <  ENTER PARAMETERS  >')
      CALL RDLINE                                                       ! GET ANOTHER USER PARAMETER LINE
      NTOKES=0
      GO TO 100
  150 CONTINUE
      NTOKES=NTOKES+1
      DO 190 I=1,NPARS                                                  ! SEE IF IT IS A PARAMETER NAME
      LEN=LENGTH(I)                                                     ! GET THE LEGAL PARAMETER NAME LENGTH
      IPARAM=I                                                          ! SAVE THE INDEX
      IF(TOKEN(1:NCHARS).EQ.NAMES(I)(1:LEN).AND.NCHARS.EQ.LEN) GO TO 200
  190 CONTINUE                                                          ! STILL LOOKING FOR THE NAME
      IF(TOKEN(1:NCHARS).EQ.'END'.AND.NCHARS.EQ.3) GO TO 1000           ! END OF PARAM LIST?
      IF(NS.NE.0) GO TO 230
      PRINT 191, TOKEN(1:NCHARS)
  191 FORMAT(' ***  ERROR  *** UADD DOES NOT HAVE A PARAMETER ',
     *  'NAMED ',A10)
      IERROR=IERROR+1
      GO TO 100
!****
!****    FOUND THE PARAMETER NAME, NOW FIND THE VALUE
!****
  200 NS=0
      iprint = isavep
      NPARAM=IPARAM
  210 CONTINUE                                                          !  NOW FIND THE VALUE
      CALL GETOKE(TOKEN,NCHARS)
      CALL UPCASE(TOKEN,NCHARS)
      NTOKES=NTOKES+1
      IF(NCHARS.GT.0) GO TO 230                                         ! END OF LINE?
      IF(NOW.EQ.1) PRINT 140                                            ! THIS ALLOWS A PARAMETER TO BE ON A DIFFERENT LINE FROM THE NAME
      CALL RDLINE                                                       ! GET ANOTHER LINE
      NTOKES=0
      GO TO 210
  230 CONTINUE
      IF(TYPES(NPARAM).NE.'A') GO TO 240
      IF(NAMES(NPARAM).EQ.'ADDWB'.AND.TOKEN(1:NCHARS).EQ.'YES')
     *    IADDWB=1
      GO TO 100
  240 CONTINUE
      CALL DCODE(TOKEN,NCHARS,AREAL,ISTAT)                              ! TRY AND DECODE IT
      IF(ISTAT.EQ.2) GO TO 420                                          ! =2 MEANS IT IS A NUMERIC
      IERROR=IERROR+1                                                   ! DCODE PRINTED AN ERROR
      GO TO 100
  420 IF(TYPES(NPARAM).EQ.'L') GO TO 500
      IF(NPARAM.LT.MULTIV) GO TO 490                                    !  IS IT A MULTIVALUED PARAMETER
      iprint = 0
      NS=NS+1                                                           !  THE TOKEN WAS A MULTI-VALUED PARAMETER
      npts = ns
      ITEMP=MULTIV
      BUF(NS+ITEMP-1)=AREAL
      GO TO 100
  490 VALS(NPARAM)=AREAL                                                !  FLOATING POINT VALUES
      GO TO 100
  500 CONTINUE                                                          !  32 BIT INTEGER VALUES
      LVALS(NPARAM)=AREAL
      GO TO 100
!****
!****   FINISHED A LIST, NOW DO THE ERROR CHECKING
!****
 1000 CONTINUE
      iprint = isavep
      IF(FNO.GT.LLNO) GO TO 1020
      PRINT 1010
 1010 FORMAT(' ***  ERROR  ***  SHOT AND RP NUMBERS MUST INCREASE.')
      IERROR=IERROR+1
 1020 IF(LNO.GE.FNO) GO TO 1030                                         ! DO THEY INCREASE IN THIS LIST
      PRINT 1010
      IERROR=IERROR+1
 1030 LLNO=LNO
!****
!****      WRITE THE PARAMETER LIST TO DISC
!****
      LBUF(1)=FNO
      LBUF(2)=LNO
      LBUF(3)=NPTS
      BUF(4) = stime
      LBUF(5)=LPRINT
      NWRDS=NPTS+5
      IF(IAND(LPRINT,1).EQ.1)  PRINT *,(LBUF(I),I=1,3),buf(4),
     *   lbuf(5),(BUF(J),J=6,NWRDS)
      CALL WRDISC(MUNIT,BUF,NWRDS)
      NLISTS=NLISTS+1
      LLNO=LNO
      NS=0
      NPTS=0
 2020 CALL GETOKE(TOKEN,NCHARS)                                         ! GET THE NEXT TOKEN
      CALL UPCASE(TOKEN,NCHARS)
      NTOKES=NTOKES+1
      IF(NCHARS.GT.0) GO TO 2030                                        !WAS IT THE END OF A LINE?
      IF(NOW.EQ.1) PRINT 140
      CALL RDLINE                                                       ! GET ANOTHER LINE
      NTOKES=0
      GO TO 2020
 2030 IF(TOKEN(1:NCHARS).NE.'END'.OR.NCHARS.NE.3) GO TO 150
      RETURN                                                            !  FINISHED ALL OF THE PARAMETERS!!!
      END
