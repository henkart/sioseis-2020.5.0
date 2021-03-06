      SUBROUTINE udeced(scr,lscr)
!                                                                           
!                        PROCESS UDECON  (USER DEFINED DECONVOLUTION)
!                        ------- ------
!
! DOCUMENT DATE: 26 February 1988
!
!     UDECON calulates and applies a Wiener deconvolution filter given a series
!  of input and desired output values or the crosscorrelation and autocorrelatio
!  functions.  The output appears in the "standard output" (print file) and may
!  be the input, output, crosscorrelation, autocorrelation, the Wiener filter
!  coefficients, the Wiener-Levinson error, the deconvolved input, or the
!  rms error of the deconvolved and input series.
!     No standard SIOSEIS input or output is available.  All input and output
!  is within this process.  i.e. the procs list should be PROCS UDECON END
!     EACH PARAMETER LIST MUST BE TERMINATED WITH THE WORD END.  THE ENTIRE SET
!  OF UDECON PARAMETERS MUST BE TERMINATED BY THE WORD END.
!
!  THE PARAMETER DICTIONARY
!  --- --------- ----------
!  INPUT  - A list of input values representing a seismic trace.  The input will
!           be correlated with itself and with the desired output.
!  DESIRE - A list of desired output values representing the desired output trac
!  AUTO   - A list of autocorrelation values that will be used instead of the
!           input and desired output.
!  CROSS  - A list of crosscorrelation values.
!  OPRINT - A bit switch indicating the printer output.  Bit x is 2**x.
!         bit 0,  A debug switch for the edit phase.
!         bit 1,  A debug switch for the execute phase.
!         bit 2,  The prediction filter is printed.
!         bit 3,  The Wiener-Levinsion error is printed.
!         bit 4,  The deconvolved input is printed.
!         bit 5,  The RMS error between the deconvolved input and the input
!                 is printed.
!             Preset = 0    e.g. oprint 32    sets bit 5 only
!                                oprint 96    sets bits 5 and 6
!                                oprint 127   sets all bits
!  NFPTS  - The number of filter points to use.  NFPTS must be larger than
!           the period to be removed.  NFPTS is also the number of lags in
!           the correlations.  If NFPTS is omitted when the autocorrelation
!           is given, NFPTS is set to the number of points in the autocor-
!           relation.
!  PREWHI - THE PERCENTAGE PREWHITENING TO ADD BEFORE FILTER DESIGN.
!           A HIGH LEVEL OF PREWHITENING REDUCES THE EFFECTIVENESS OF THE
!           FILTER.  SOME LEVEL OF PREWHITENING IS NEEDED IN ORDER FOR THE
!           FILTER TO BE STABLE.  PREWHITENING IS LIKE PERFORMING A BANDPASS
!           FILTER BEFORE DECON.
!           PRESET=25.
!  END    - TERMINATES EACH PARAMETER LIST.
!
!  WRITTEN AND COPYRIGHTED (C) BY:
!  PAUL HENKART, SCRIPPS INSTITUTION OF OCEANOGRAPHY, August 1986
!  ALL RIGHTS ARE RESERVED BY THE AUTHOR.  PERMISSION TO COPY OR REPRODUCE THIS
!  SUBROUTINE, BY COMPUTER OR OTHER MEANS, MAY BE OBTAINED ONLY FROM THE AUTHOR.
!
!  mod 22 Dec 97 - Change MIN0 to MIN
!
!  ARGUMENTS:
!  BUF    - A SCRATCH ARRAY AT LEAST 60 32 BIT WORDS LONG.
!  LBUF   - THE SAME ARRAY BUT THE 32 BIT INTEGER EQUIVALENT.  NEEDED
!           BECAUSE FORTRAN DOESN'T ALLOW EQUIVALENCING OF ARGUMENTS.
!
      PARAMETER (NPARS = 8)                                             ! THE NUMBER OF USER PARAMETERS
      CHARACTER*6 NAMES(NPARS)
      CHARACTER*1 TYPE(NPARS)
      DIMENSION LENGTH(NPARS)
      CHARACTER*80 TOKEN
      DIMENSION VALS(NPARS),LVALS(NPARS)
      EQUIVALENCE (vals(1),lvals(1))
      COMMON /EDITS/ IERROR,IWARN,IRUN,NOW,ICOMPT
      COMMON /UDECON/ MUNIT,NLISTS
      DIMENSION scr(111),lscr(111)
      INTEGER oprint
      REAL input
!
!
      EQUIVALENCE (LPRINT,LVALS(1)),
     2            (PREWHI,VALS(2)),
     3            (nfpts,lvals(3)),
     4            (oprint,lvals(4)),
     5            (input,vals(5)),
     6            (desire,vals(6)),
     7            (auto,vals(7)),
     8            (cross,vals(8))
      DATA NAMES/'LPRINT','PREWHI','NFPTS ','OPRINT','INPUT ',
     *           'DESIRE','AUTO  ','CROSS '/
      DATA LENGTH/6,6,5,6,5,6,4,5/
      DATA TYPE/'L','F','L','L',4*'F'/
!****
!****      SET THE PRESETS
!****
      LPRINT=0
      PREWHI=25.
      nfpts = 0
      ninpts = 0                                                        !  count the number of input points given
      ndesir = 0                                                        !  count the number of desired output points given
      nauto = 0                                                         !  count the number of autocorrelation lags given
      ncross = 0                                                        !  count the number of crosscorrelation lags given
      oprint = 0                                                        !  set the output print
      ntypes = 0
!****
!*****    GET A PARAMETER FILE
!****
      CALL GETFIL(1,MUNIT,token,ISTAT)
!****
!****   THE CURRENT COMMAND LINE IN THE SYSTEM BUFFER MAY HAVE THE PARAMETERS.
!****   GET A PARAMETER LIST FROM THE USER.
!****
      ipoint = 1
      NS=0
      NLISTS=0
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
      IF(TOKEN(1:NCHARS).EQ.'END'.AND.NCHARS.EQ.3) GO TO 1000           ! END OF PARAM LI
      IF(NS.NE.0) GO TO 230
      PRINT 191, TOKEN(1:NCHARS)
  191 FORMAT(' ***  ERROR  *** DECON DOES NOT HAVE A PARAMETER ',
     *  'NAMED ',A10)
      IERROR=IERROR+1
      GO TO 100
!****
!****    FOUND THE PARAMETER NAME, NOW FIND THE VALUE
!****
  200 CONTINUE
      NS = 0
      NPARAM=IPARAM
  210 CONTINUE                                                          !  NOW FIND THE VALUE
      CALL GETOKE(TOKEN,NCHARS)
      CALL UPCASE(TOKEN,NCHARS)
      NTOKES=NTOKES+1
      IF( NCHARS .LE. 0 ) THEN                                           ! END OF LINE?
          IF( NOW .EQ. 1 ) PRINT 140                                     ! THIS ALLOWS A PARAMETER TO BE ON A DIFFE
          CALL RDLINE                                                    ! GET ANOTHER LINE
          NTOKES=0
          GOTO 210
      ENDIF
 230  CALL DCODE(TOKEN,NCHARS,AREAL,ISTAT)                              ! TRY AND DECODE IT
      IF( ISTAT .NE. 2 ) THEN                                           ! =2 MEANS IT IS A NUMERIC
          IERROR = IERROR + 1                                           ! DCODE PRINTED AN ERROR
          GOTO 100
      ENDIF
      IF( TYPE(NPARAM) .EQ. 'L') THEN
          lvals(nparam) = areal
          GOTO 100
      ENDIF
      IF( type(nparam) .EQ. 'F' ) THEN
          ns = ns + 1
          IF( NAMES(NPARAM) .EQ. 'INPUT' ) THEN
              IF( ninpts .EQ. 0 ) THEN
                  ipoint = ipoint + 1
                  lscr(npars+ipoint) = 9991                             ! word 1 is the array type
                  ntypes = ntypes + 1
                  idumb = npars+ipoint
              ENDIF
              ipoint = ipoint + 1
              scr(npars+ipoint) = areal
              ninpts = ninpts + 1
              GOTO 100
          ENDIF
          IF( names(nparam) .EQ. 'DESIRE' ) THEN
              IF( ndesir .EQ. 0 ) THEN
                  ipoint = ipoint + 1
                  lscr(npars+ipoint) = 9992
                  ntypes = ntypes + 1
                  idumb = npars+ipoint
              ENDIF
              ipoint = ipoint + 1
              scr(npars+ipoint) = areal
              ndesir = ndesir + 1
              GOTO 100
          ENDIF
          IF( names(nparam) .EQ. 'AUTO' ) THEN
              IF( nauto .EQ. 0 ) THEN
                  ipoint = ipoint + 1
                  lscr(npars+ipoint) = 9993
                  ntypes = ntypes + 1
                  idumb = npars+ipoint
              ENDIF
              ipoint = ipoint + 1
              scr(npars+ipoint) = areal
              nauto = nauto + 1
              GOTO 100
          ENDIF
          IF( names(nparam) .EQ. 'CROSS' ) THEN
              IF( ncross .EQ. 0 ) THEN
                  ipoint = ipoint + 1
                  lscr(npars+ipoint) = 9994
                  ntypes = ntypes + 1
                  idumb = npars+ipoint
              ENDIF
              ipoint = ipoint + 1
              scr(npars+ipoint) = areal
              ncross = ncross + 1
              GOTO 100
          ENDIF
          vals(nparam) = areal
          GOTO 100
      ENDIF
!****
!****   FINISHED A LIST, NOW DO THE ERROR AND VALIDITY CHECKS
!****
 1000 CONTINUE                                                           ! MAKE SURE ALL SHOT & RP NUMBERS INCREASE
      IF( ninpts .NE. 0 ) THEN                                           ! make sure desire is given when input is
          IF( ndesir .EQ. 0 ) THEN
              PRINT *,' ***  ERROR  ***  DESIRE must be given with ',
     *                'INPUT.'
              ierror = ierror + 1
          ENDIF
          IF( nauto+ncross .NE. 0 ) THEN
              PRINT *,' ***  WARNING  ***  Incompatible inputs (INPUT,',
     *                ' AUTO, CROSS).'
              iwarn = iwarn + 1
          ENDIF
      ENDIF
      IF( ndesir .NE. 0 ) THEN                                           ! make sure input is given when desire is
          IF( ninpts .EQ. 0 ) THEN
              PRINT *,' ***  ERROR  ***  INPUT must be given with ',
     *                ' DESIRE'
              ierror = ierror + 1
          ENDIF
      ENDIF
      IF( nauto .NE. 0 ) THEN                                            ! make sure cross is given when auto is
          IF( ncross .EQ. 0 ) THEN
              PRINT *,' ***  ERROR  ***  CROSS must be given with AUTO.'
              ierror = ierror + 1
          ENDIF
      ENDIF
      IF( ncross .NE. 0 ) THEN                                           ! make sure auto is given when cross is
          IF( nauto .EQ. 0 ) THEN
              PRINT *,' ***  ERROR  ***  AUTO must be given with CROSS.'
              ierror = ierror + 1
          ENDIF
      ENDIF
      IF( ncross+nauto .NE. 0 ) THEN                                    ! is either given?
          IF( ncross .NE. nauto ) THEN
              PRINT *, ' ***  WARNING  ***  The autocorrelation and',
     *        ' crosscorrelation are not the same length.'
              iwarn = iwarn + 1
              IF( nfpts .EQ. 0 ) THEN
!                  nfpts = MIN0(nauto,ncross)
                  nfpts = MIN(nauto,ncross)
                  PRINT *,' SIOSEIS will use ',nfpts,' filter points.'
              ENDIF
          ENDIF
          IF( nfpts .EQ. 0 ) nfpts = nauto
          IF( nfpts .NE. nauto ) THEN
              PRINT *,' ***  WARNING  ***  The filter length is not ',
     *      'the same length as the autocorrelation.'
              iwarn = iwarn + 1
          ENDIF
      ENDIF
!****
!****      WRITE THE PARAMETER LIST TO DISC
!****
      lscr(1) = lprint
      scr(2) = prewhi
      lscr(3) = nfpts
      lscr(4) = oprint
      lscr(5) = ninpts
      lscr(6) = ndesir
      lscr(7) = nauto
      lscr(8) = ncross
      lscr(9) = ntypes
      IF( IAND(LPRINT,1) .EQ. 1 )  PRINT *,lscr(1),scr(2),
     *    (lscr(i),i=3,10)
      nwrds = npars + ipoint
      CALL WRDISC(MUNIT,scr,NWRDS)
      NLISTS=NLISTS+1
      NS=0
      ipoint = 0
      ninpts = 0
      ndesir= 0
      nauto = 0
      ncross = 0
 2020 CALL GETOKE(TOKEN,NCHARS)                                         ! GET THE NEXT TOKEN
      CALL UPCASE(TOKEN,NCHARS)
      NTOKES=NTOKES+1
      IF( NCHARS .LE. 0 ) THEN                                          ! WAS IT THE END OF A LINE?
          IF( NOW .EQ. 1 ) PRINT 140
          CALL RDLINE                                                   ! GET ANOTHER LINE
          NTOKES=0
          GOTO 2020
      ENDIF
      IF(TOKEN(1:NCHARS).NE.'END'.OR.NCHARS.NE.3) GO TO 150
      RETURN                                                            !  FINISHED ALL OF THE PARAMETERS!!!
      END
