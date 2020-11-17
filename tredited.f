      SUBROUTINE tredited( buf, lbuf, ibuf )
!*******
!*******      MAKE SURE TO USE THE RIGHT COMMON BLOCK NAME
!*******      despike and tredit 
!
!                       PROCESS DESPIKE
!                       ------- -------
!
!  Document date: 27 June 1995
!
!      Process DESPIKE removes spikes from the seismic trace.  Several
!  different despiking algorithms are available.  The use of parameters
!  THRES or FAC or QUART indicates which algorithim is used.
!
!      Only one fno/lno parameter list may be given.
!
!  PARAMETER DICTIONARY
!  --------- ----------
!  THRES  - The minimum and maximum amplitude threshold.  Amplitudes
!           that exceed the treshold are replaced through linear
!           interpolation of the adjacent "good" amplitudes.
!           Preset = 0 0    e.g.   thres -1.e6 1.e6
!  FAC    - The tolerance factor used in the Trehu/Sutton method
!           described in Marine Gephysical Researches 16: 91-103, 1994.
!           The algorithm is based on a five point moving window.
!           Preset = 0.   e.g. fac 5.
!  QUART  - Quartile, amplitudes above this quartile are replaced by the
!           signed quartile value. ( QUART must be between 1 and 100).
!           This algorithm came from CWP/SU.
!           Preset = 0   e.g. quart 99.
!  FNO    - The first shot/rp number the parameter list applies to.
!           Preset = the first shot/rp received.    e.g.   fno 101
!  LNO    - The last shot/rp number the parameter list applies to.
!           Preset = the last shot/rp received.    e.g.   lno 101
!  FTR    - The first trace number the parameter list applies to.
!           Preset = the first trace of each shot/rp.    e.g.   ftr 10
!  LTR    - The last trace number the parameter list applies to.
!           Preset = the last trace of each shot/rp.    e.g.   ltr 10
!  END    - Terminates the parameter list.
!
!
!  Copyright (C) 1992 The Regents of the University of California
!  ALL RIGHTS RESERVED.
!
!  Mod 14 June 1994 - Add parameter TYPE and the Trehu/Sutton method.
!  Modified 27 June 1995 - Michael Holzrichter, ODP - Added QUART
!  Mod 17 Aug 95 by pch to remove parameter TYPE by setting it
!     depending on whether THRES or FAC or QUART are given
!  Mod 1 Apr 96 - Threshold won't work when thres(1) + thres(2) = 0
!  Mod 22 July 97 - Add KILL, ALPHA, SET, ADDWB, VEL, SES, SEL
!  Mod 2 April 98 - Add MEDIAN.
!  Mod 27 May 98 - If thres(1) = -thres(2), it failed the edit!
!  Mod 1 Apr 99 - Add parameter MINVAL
!  Mod 11 Sept 00 - Add LIMITS, IHDR, LHDR, IHDR
!  Mod 15 Feb 02 - Add WINLEN, HCYCLE
!  Mod 24 Aug 05 - Add ENDMUTE
!  Mod 1 Sep 05 - Add KILL FLAGONLY
!  Mod 13 Oct 09 - Add PASS
!                - Add SOLRAT and make it the same as FAC
!                - Turn kill on when SES is given.
!  Mod 31 Jul 12 - Add  MAXVAL, MAXTYPE, MINTYPE, DEVPCT
!
      PARAMETER ( npars = 32 )                                          ! the number of user parameters
      CHARACTER*80 token
      CHARACTER*7 names(npars)
      DIMENSION scr(npars), lscr(npars), lbuf(111), buf(111)
      INTEGER*2 ibuf(111)
      EQUIVALENCE (scr(1),lscr(1))
      COMMON /edits/ ierror, iwarn, irun, now, icompt
!*******
!*******      MAKE SURE TO USE THE RIGHT COMMON BLOCK NAME
!*******      despike and tredit 
      COMMON /tredit/ thres(2), lprint, fno, lno, ftr, ltr, type, fac,
     &       quart, kill, alpha, set(4), addwb, vel, ses(2), sel(2),
     &       lunmed, valmin, limits(2), ihdr, lhdr, hdr, winlen, hcycle,
     &       endmute, pass(2), valmax, mintype, maxtype, devpct
      REAL limits
      INTEGER fno, lno, ftr, ltr, type, addwb
      DATA names / 'THRES ', 'LPRINT', 'FNO   ', 'LNO   ', 'FTR   ',
     &             'LTR   ', 'TYPE  ', 'FAC   ', 'QUART ', 'KILL  ',
     &             'ALPHA ', 'SET   ', 'ADDWB ', 'VEL   ', 'SES   ',
     &             'SEL   ', 'SETS  ', 'MEDIAN', 'MINVAL', 'LIMITS',
     &             'IHDR  ', 'LHDR  ', 'HDR   ', 'WINLEN', 'HCYCLE',
     &             'ENDMUTE','PASS  ', 'SOLRAT', 'MAXVAL', 'MINTYPE',
     &             'MAXTYPE', 'DEVPCT'/
!**** 
!****    Set the parameter presets and various variable presets
!****
      lprint = 0
      fno = 0
      lno = 9999999
      ftr = 0
      ltr = 999999
      type = 0
      fac = 0.
      quart = 0
      kill = 0
      alpha = 1.
      nthres = 0
      nsets = 0
      nses = 0
      nsel = 0
      nlimits = 0
      npass = 0
      DO i = 1, 2
         thres(i) = 0.
         set(i) = 0.
         set(i+2) = 0.
         ses(i) = 0.
         sel(i) = 0.
         limits(i) = -999999.
         pass(i) = 0.
      ENDDO
      addwb = 0
      vel = 0.
      lunmed = 0
      valmin = 0.
      ihdr = 0
      lhdr = 0
      hdr = 0.
      winlen = 0.
      hcycle = 0.
      endmute = 0.
      maxval = 0.
      mintype = 1
      maxtype = 1
      devpct = 0.
!****
!****     get the user's parameters
!****
      ntokes = 0                                                        ! count the tokens
  100 CONTINUE
      CALL getoke( token, nchars )                                      ! get a token and it's length
      IF( nchars .EQ. 0 ) THEN                                          ! anything there?
          CALL rdline                                                   ! nope, get another line
          ntokes = 0
          GOTO 100
      ENDIF
      CALL upcase( token, nchars )                                      ! convert parameter names to upper case
  110 ntokes = ntokes + 1
      DO 200 nparam = 1, npars
         IF( token(1:nchars) .EQ. names(nparam) ) THEN                  ! find the parameter name in our list
  120        CALL getoke( token, nchars )                               ! get the value
             IF( nchars .EQ. 0 ) THEN
                 CALL rdline
                 ntokes = 0
                 GOTO 120
             ENDIF
             ntokes = ntokes + 1
             CALL upcase( token, nchars )
!****        This next contortion is because THRES may have 1 or 2
!****        values.  If only 1 value, then the next thing should be
!****        be a valid parameter name.
             DO i = 1, npars
                IF( token(1:nchars) .EQ. names(i) ) GOTO 110
             ENDDO
             IF( token(1:nchars) .EQ. 'END' ) GOTO 200
             IF( names(nparam) .EQ. 'KILL' ) THEN
                 IF( token(1:1) .EQ. 'Y' ) kill = 1
                 IF( token(1:6) .EQ. 'INSIDE' ) kill = 2
                 IF( token(1:7) .EQ. 'OUTSIDE' ) kill = 3
                 IF( token(1:8) .EQ. 'FLAGONLY' ) kill = 4
                 IF( token(1:2) .NE. 'NO' .AND. kill .EQ. 0 ) THEN
                     PRINT *,' ***  ERROR  ***  Bad KILL value of: ',
     &                    token(1:nchars)
                     ierror = ierror + 1
                 ENDIF
                 GOTO 100
             ENDIF
             IF( names(nparam) .EQ. 'ADDWB' ) THEN
                 IF( token(1:1) .EQ. 'Y' ) addwb = 1
                 GOTO 100
             ENDIF
             IF( names(nparam) .EQ. 'MEDIAN' ) THEN
                 IF( token(1:1) .EQ. 'Y' ) THEN
                     CALL getfil( 1, lunmed, token, istat )
                 ELSE
                   PRINT *,' ***  ERROR  ***  MEDIAN must be YES or NO.'
                     ierror = ierror +1
                 ENDIF
                 IF( lbuf(7) .EQ. 0 ) THEN
                     PRINT *,' ***  WARNING  ***  MEDIAN "stack" ',
     &                     ' requires "gathered" data.'
                 ENDIF
                 GOTO 100
             ENDIF
             IF( names(nparam) .EQ. 'MINTYPE' ) THEN
                 IF( token(1:6) .EQ. 'ABSVAL' ) mintype = 1
                 IF( token(1:4) .EQ. 'SDEV' ) mintype = 2
                 IF( token(1:3) .EQ. 'VAR' ) mintype = 3
                 IF( token(1:4) .EQ. 'SKEW' ) mintype = 4
                 IF( token(1:4) .EQ. 'KURT' ) mintype =5
                 GOTO 100
             ENDIF
             IF( names(nparam) .EQ. 'MAXTYPE' ) THEN
                 IF( token(1:6) .EQ. 'ABSVAL' ) maxtype = 1
                 IF( token(1:4) .EQ. 'SDEV' ) maxtype = 2
                 IF( token(1:3) .EQ. 'VAR' ) maxtype = 3
                 IF( token(1:4) .EQ. 'SKEW' ) maxtype = 4
                 IF( token(1:4) .EQ. 'KURT' ) maxtype =5
                 GOTO 100
             ENDIF
             CALL dcode( token, nchars, areal, istat )                  ! convert the alpha number to an internal machine number
             IF( istat .NE. 2 ) ierror = ierror + 1                     ! was the an error decoding it?
             IF( names(nparam) .EQ. 'THRES' ) THEN
                 nthres = nthres + 1
                 thres(nthres) = areal
                 IF( nthres .EQ. 1 ) GOTO 120
             ENDIF
             IF( names(nparam) .EQ. 'LPRINT' ) lprint = NINT(areal)
             IF( names(nparam) .EQ. 'FNO' ) fno = NINT(areal)
             IF( names(nparam) .EQ. 'LNO' ) lno = NINT(areal)
             IF( names(nparam) .EQ. 'FTR' ) ftr = NINT(areal)
             IF( names(nparam) .EQ. 'LTR' ) ltr = NINT(areal)
             IF( names(nparam) .EQ. 'TYPE' ) type = NINT(areal)
             IF( names(nparam) .EQ. 'FAC' ) fac = areal
!*****   make solrat the same as fac
             IF( names(nparam) .EQ. 'SOLRAT' ) fac = areal
             IF( names(nparam) .EQ. 'QUART' ) quart = areal / 100.
             IF( names(nparam) .EQ. 'ALPHA' ) alpha = areal
             IF( names(nparam) .EQ. 'VEL' ) vel = areal
             IF( names(nparam) .EQ. 'MINVAL' ) valmin = areal
             IF( names(nparam) .EQ. 'IHDR' ) ihdr = NINT(areal)
             IF( names(nparam) .EQ. 'LHDR' ) lhdr = NINT(areal)
             IF( names(nparam) .EQ. 'HDR' ) hdr = areal
             IF( names(nparam) .EQ. 'WINLEN' ) winlen = areal
             IF( names(nparam) .EQ. 'HCYCLE' ) hcycle = areal
             IF( names(nparam) .EQ. 'ENDMUTE' ) endmute = areal
             IF( names(nparam) .EQ. 'MAXVAL' ) valmax = areal
             IF( names(nparam) .EQ. 'DEVTYPE' ) devtype = areal
             IF( names(nparam) .EQ. 'SET' .OR.
     &           names(nparam) .EQ. 'SETS' ) THEN
                 nsets = nsets + 1
                 set(nsets) = areal
                 GOTO 120
             ENDIF
             IF( names(nparam) .EQ. 'SES' ) THEN
                 nses = nses + 1
                 ses(nses) = areal
                 kill = 1
                 IF( nses .EQ. 1 ) GOTO 120
             ENDIF
             IF( names(nparam) .EQ. 'SEL' ) THEN
                 nsel = nsel + 1
                 sel(nsel) = areal
                 IF( nsel .EQ. 1 ) GOTO 120
             ENDIF
             IF( names(nparam) .EQ. 'LIMITS' ) THEN
                 nlimits = nlimits + 1
                 limits(nlimits) = areal
                 IF( nlimits .EQ. 1 ) GOTO 120
             ENDIF
             IF( names(nparam) .EQ. 'PASS' ) THEN
                 npass = npass + 1
                 pass(npass) = areal
                 IF( npass .EQ. 1 ) GOTO 120
             ENDIF
             GOTO 100
         ENDIF
  200 CONTINUE
      IF( token(1:nchars) .NE. 'END') THEN
          PRINT *,' ***  ERROR  ***  No such parameter as ',
     *      token(1:nchars)
          ierror = ierror + 1
          GOTO 100
      ENDIF
!****
!****    Do some ERROR checking
!****
!      every sample is checked for NaN whether anything else is done or not.
!      IF( thres(1) + fac + quart + ses(2)+ sel(2) + lunmed + valmin +
!     &    winlen .EQ. 0.) THEN
!          IF( limits(1) .EQ. -999999. ) THEN
!            PRINT *,' ***  ERROR  ***  No DESPIKE parameters specified.'
!              ierror = ierror + 1
!          ENDIF
!      ENDIF
      IF( type .EQ. 1 .AND. thres(1) .EQ. 0. .AND. thres(2) .EQ. 0) THEN
          PRINT *,' ***  ERROR  ***  THRES must be given.'
          ierror = ierror + 1
      ENDIF
      IF( quart .LT. 0 .OR. quart .GT. 1. ) THEN
          PRINT *,' ***  ERROR  ***  QUART must be between 1 and 100.'
          ierror = ierror + 1
      ENDIF
      IF( ses(2) .NE. 0 .AND. sel(2) .EQ. 0 ) THEN
          PRINT *,' *** ERROR ***  SEL must be given when SES is used.'
          ierror = ierror + 1
      ENDIF
      IF( sel(2) .NE. 0 .AND. ses(2) .EQ. 0 ) THEN
          PRINT *,' *** ERROR ***  SES must be given when SEL is used.'
          ierror = ierror + 1
      ENDIF
      IF( sel(2) + ses(2) .NE. 0 .AND. fac .EQ. 0. ) THEN
          PRINT *,' ***  ERROR  ***  FAC must be given with SES/SEL.'
          ierror = ierror + 1
      ENDIF
      IF( nsets .EQ. 1 .OR. nsets .EQ. 3 ) THEN
          PRINT *,' ***  ERROR  ***  There must be 2 or 4 SETS values.'
          ierror = ierror + 1
      ENDIF
      IF( nsets .EQ. 4 .AND. thres(2) .NE. 0 ) THEN
          PRINT *,
     &    ' ***  WARNING  ***  Only the first THRES value will be used.'
          iwarn = iwarn + 1
      ENDIF
      IF( pass(1) .LT. 0 ) THEN
          PRINT *, ' ***  ERROR  ***  Bad PASS.'
          ierror = ierror + 1
      ENDIF
      IF( pass(1) .GT. 0 .AND. pass(2) .LE. 0 ) THEN
          PRINT *, ' ***  ERROR  ***  Bad PASS.'
          ierror = ierror + 1
      ENDIF
      IF( limits(1) .NE. -999999. ) THEN
          IF( ihdr + lhdr + hdr .EQ. 0 ) THEN
              PRINT *,
     &' ***  ERROR  ***  IHDR or LHDR or HDR must be given with LIMITS.'
              ierror = ierror + 1
          ENDIF
      ENDIF
      IF( ihdr + lhdr + hdr .NE. 0 ) THEN
          IF( limits(1) .EQ. -999999. ) THEN
              PRINT *,
     &' ***  ERROR  ***  LIMITS must be given with IHDR or LHDR or HDR.'
              ierror = ierror + 1
          ENDIF
      ENDIF
      IF( winlen .LT. 0. ) THEN
          PRINT *,' ***  ERROR  ***  WINLEN must be positive.'
          ierror = ierror + 1
      ENDIF
      IF( hcycle .LT. 0. ) THEN
          PRINT *,' ***  ERROR  ***  HCYCLE must be positive.'
          ierror = ierror + 1
      ENDIF
      IF( winlen .NE. 0. .AND. hcycle .EQ. 0. ) hcycle = winlen / 2
      IF( IAND(lprint,1) .NE. 0 ) THEN
          PRINT *,' thres=',thres,' fno=',fno,' lno=',lno,
     &            ' ftr=',ftr,' ltr=',ltr,' fac=',fac,' quart=',quart
          PRINT *,' set=',set,' lunmed=',lunmed,' minval=',valmin
          PRINT *,' limits=',limits,' ihdr=',ihdr,' lhdr=',lhdr,
     &            ' hdr=',hdr,' winlen=',winlen,' hcycle=',hcycle
          PRINT *,' endmute=',endmute,' pass=',pass
      ENDIF
!****
!****    finish up the parameter reading
!****
 2000 CONTINUE
      CALL getoke( token, nchars )                                       ! get the next token
      IF( nchars .LE. 0 ) THEN
          IF( now .EQ. 1 ) PRINT *,' <  ENTER PARAMETERS  >'
          CALL rdline                                                          ! get a new line of parameters
          ntokes = 0
          GOTO 2000
      ENDIF
      CALL upcase( token, nchars )
      ntokes = ntokes + 1
      IF( token(1:nchars) .NE. 'END' .OR. nchars .NE. 3 ) GOTO 100
      RETURN
      END
