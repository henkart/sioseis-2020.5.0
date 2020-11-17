      SUBROUTINE gained( tgp, num_gain )
!
!                       PROCESS GAINS
!                       ------- -----
!
!  Document date: 14 July 1993
!
!      Process GAINS applies a gain function.  Chapter 4 of Claerbout's
!  "Imaging the Earth's Interior" mentions several of the gain
!  functions implemented.
!                                               
!      At least one parameter list must be given, even if no parameters 
!  are specified, in order that the parameter presets be set.  e.g.
!      gains
!              end
!      end
!
!  PARAMETER DICTIONARY
!  --------- ----------
!  TYPE   - The type of gain to apply.
!         = 1,  a(i) = a(i) * (t*1000.)**alpha    ( USGS gain )
!           where a(i) is the trace and t is the time of the trace
!           sample in seconds.
!           Restrictions:  All traces must have the same start time.
!           Preset = 1
!         = 2,  a(i) = a(i) * (ABS(range)/SIGN(rscale,range)) ** alpha
!           when ABS(range) .GE. rscale;
!           where range is the range in the SEGY header, rscale
!           and alpha are given by the user.  SIGN is the Fortran SIGN
!           function which means that ABS(rscale) is used when range is
!           positive and -ABS(rscale) is used when range is neagtive.
!         = 3,  a(i) = a(i) * t ** alpha
!         = 4,  a(i) = a(i) ** alpha
!         = 5,  a(i) = a(i) * e ** (alpha * t)
!         = 6,  a(i) = SGN(a(i)) * ABS(a(i)) ** alpha
!         = 7,  a(i) = SQRT(a(i*2-1)**2+a(i*2)**2)  (modulus of complex
!                      trace)
!         = 8,  a(i) = a(i) * ABS(range/rscale) ** alpha
!           when ABS(range) .GE. rscale;
!           where range is the range in the SEGY header, rscale
!           and alpha are given by the user. 
!
!  Additional Parameters:
!  ---------- -----------
!  ETIME  - The end time of the gain function.  Data after the end time 
!           will receive the gain of the end time.
!           Preset = the last time of the first trace.  e.g.  etime 4.
!  ALPHA  - The exponent used in TYPE 1 and 2 gain.
!           Preset = 1.
!  RSCALE - The range scalar used in TYPE 2 gain.
!           PRESET = 1.
!  SUBWB  - Subtract water bottom time switch.  Type 3 and 5 ONLY.
!         = YES, The water bottom time is subtracted from the data time
!                in the gain types that use time as a vaiable.  e.g.
!                a(i) = a(i) * t ** alpha     becomes
!                a(i) = a(i) * (t-wbt) ** alpha
!  FNO    - The first shot/rp number the parameter list applies to.
!           Data (shots/rps) before FNO WILL NOT HAVE GAINS APPLIED.
!           Preset = the first shot/rp received.    e.g.   fno 101
!  LNO    - The last shot/rp number the parameter list applies to.
!           Data (shots/rps) AFTER LNO WILL NOT HAVE GAINS APPLIED.
!           Preset = the last shot/rp received.    e.g.   lno 101
!
!
!  Copyright (C) 1990 Seismic Reflection Processors, Solana Beach, CA.
!  ALL RIGHTS RESERVED.
!
!  Mod 12 Feb 91 to add subwt and add types 3-6, and multiple lists.
!  mod 16 sep 92 - redfine type scaling
!  mod 14 Jul 93 - Add fno and lno to EXCLUDE data from being gained.
!  mod Aug 96 - Add gain type 8
!  mod Jul 00 - Add gain type 9 and parameters TGP and ADDWB.
!  mod Oct 00 - Add TADD and TMULT
!  mod 6 Nov 01 - Add error check for tgps not being in pairs.
!  mod 14 Jan 03 - Bad error check on ntgp
!  mod 9 Jun 03 - Add WINLEN (running average window length)
!  mod 1 Dec 03 - Make multiple lists work with TGP.
!  mod 28 Jan 04 - Bad error check from 1 Dec update.
!  mod 14 Apr. 04 - Allow SUBWB on TYPE 9 (TGP).
!  mod 5 Jan 07 - Add num_gains - allow 3 process gains
!  mod 18 Feb 10 - Add check for type 4 and non-integer alpha
!  mod 14 Jun 11 - Add type 10 (20log(t*v))
!                - Add parameter v
!
      PARAMETER ( npars = 14 )                                          ! the number of user parameters
      PARAMETER ( max_gains = 3 )
      CHARACTER*80 token
      CHARACTER*6 names(npars)
      DIMENSION scr(npars), lscr(npars), tgp(1)
      EQUIVALENCE (scr(1),lscr(1))
      COMMON /edits/ ierror, iwarn, irun, now, icompt
      COMMON /gains/ igunit(max_gains), nglists(max_gains),
     &    ngwrds(max_gains)
      DATA names / 'TYPE  ','ETIME ','LPRINT','ALPHA ','RSCALE',
     &             'SUBWB ','FNO   ','LNO   ','TGP   ','ADDWB ',
     &             'TADD  ','TMULT ','WINLEN','V     '/
      DATA lastlno/0/
!**** 
!****    Set the parameter presets and various variable presets
!****
      itype = 0
      etime = -1.
      lprint = 0
      alpha = 1.
      rscale = 1.
      isubwb = 0
      nglists(num_gain) = 0
      fno = 0
      lno = 9999999
      ntgp = 0
      iaddwb = 0
      tadd = 0.
      tmult = 1.
      winlen = 0.
      v = 1500.
      CALL getfil( 1, igunit(num_gain), token, istat )                            ! get a file for the GAINS parameters
!****
!****     get the user's parameters -  there must be something, at least an "end"
!****
      ntokes = 0                                                        ! count the tokens
  100 CONTINUE
      CALL getoke( token, nchars )                                      ! get a token and it's length
  101 CALL upcase( token, nchars )                                      ! convert parameter names to upper case
      IF( nchars .EQ. 0 ) THEN                                          ! anything there?
          CALL rdline                                                   ! nope, get another line
          ntokes = 0
          GOTO 100
      ENDIF
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
             IF( names(nparam) .EQ. 'SUBWB' ) THEN
                 IF( token(1:1) .EQ. 'Y' ) isubwb = 1
                 GOTO 100
             ENDIF
             IF( names(nparam) .EQ. 'ADDWB' ) THEN
                 IF( token(1:1) .EQ. 'Y' ) iaddwb = 1
                 GOTO 100
             ENDIF
             CALL dcode( token, nchars, areal, istat )                  ! convert the alpha number to an internal machine number
             IF( istat .NE. 2 ) ierror = ierror + 1                     ! was the an error decoding it?
             IF( names(nparam) .EQ. 'ETIME' ) etime = areal
             IF( names(nparam) .EQ. 'TYPE' ) itype = NINT(areal)
             IF( names(nparam) .EQ. 'LPRINT' ) lprint = NINT(areal)
             IF( names(nparam) .EQ. 'ALPHA' ) alpha = areal
             IF( names(nparam) .EQ. 'RSCALE' ) rscale = areal
             IF( names(nparam) .EQ. 'FNO' ) fno = NINT(areal)
             IF( names(nparam) .EQ. 'LNO' ) lno = NINT(areal)
             IF( names(nparam) .EQ. 'TADD' ) tadd = areal
             IF( names(nparam) .EQ. 'TMULT' ) tmult = areal
             IF( names(nparam) .EQ. 'WINLEN' ) winlen = areal
             IF( names(nparam) .EQ. 'TGP') THEN
                 ntgp = ntgp + 1
                 tgp(ntgp) = areal
                 itype = 9
             ENDIF
             IF( names(nparam) .EQ. 'V' ) v = areal
             GOTO 100
         ENDIF
  200 CONTINUE
      IF( token(1:nchars) .NE. 'END' .AND. ntgp .NE. 0 ) THEN
          CALL dcode( token, nchars, areal, istat )
          IF( istat .NE. 2 ) ierror = ierror + 1
          ntgp = ntgp + 1
          tgp(ntgp) = areal
          GOTO 100
      ENDIF
      IF( token(1:nchars) .NE. 'END' ) THEN
          PRINT *,' ***  ERROR  ***  No such parameter as ',
     *      token(1:nchars)
          ierror = ierror + 1
          GOTO 100
      ENDIF
!****
!****    Do some ERROR checking
!****
      IF( etime .LT. -1 .OR. etime .GT. 100. ) THEN
          PRINT *,' ***  ERROR  ***  ETIME must be between 0 and 100.'
          ierror = ierror + 1
      ENDIF
      IF( itype .EQ. 0 .AND. winlen .EQ. 0. ) THEN
          PRINT *,' ***  ERROR  ***  TYPE not given.'
          ierror = ierror + 1
      ENDIF
      IF( itype .EQ. 0 .AND. winlen .NE. 0 ) THEN
          PRINT *,' ***  WARNING  ***  Amplitude runnig average only.'
          iwarn = iwarn + 1
      ENDIF
      IF( itype .LT. 0 .OR. itype .GT. 10 ) THEN
          PRINT *,' ***  ERROR  ***  Illegal TYPE value of ',itype
          ierror = ierror + 1
      ENDIF
      IF( itype .EQ. 4 ) THEN
          itemp = NINT(alpha)
          IF( FLOAT(itemp) .NE. alpha ) THEN
              PRINT *,' ***  ERROR  *** Type 4 requires integer alpha.'
              PRINT *,' Did you mean to use TYPE 6?'
          ENDIF
      ENDIF
      IF( isubwb .NE. 0 .AND. itype .NE. 3 .AND. itype .NE. 5 .AND.
     &    itype .NE. 9 ) THEN
         PRINT *,' ***  ERROR  *** SUBWB valid with TYPE 3, 5 & 9 only.'
          ierror = ierror + 1
      ENDIF
      IF( iaddwb .NE. 0 .AND. itype .NE. 9 .AND. itype .NE. 10 ) THEN
          PRINT *,' ***  ERROR  *** ADDWB valid with TYPE 9 or 10 only.'
          ierror = ierror + 1
      ENDIF
!      IF( nglists(num_gain) .GT. 0 ) THEN
!          PRINT *,' ***  ERROR  ***  Only one GAIN list is permitted.'
!          ierror = ierror + 1
!      ENDIF
      IF( ntgp / 2 * 2 .NE. ntgp ) THEN
          PRINT *,' ***  ERROR  ***  TGPs must be in pairs.'
          ierror = ierror + 1
      ENDIF
      IF( nglists(num_gain) .GT. 0 ) THEN
          IF( fno .NE. lastlno + 1 ) THEN
              PRINT *,' ***  ERROR  ***  No spatial allowed in GAINS.'
              ierror = ierror + 1
          ENDIF
      ENDIF
!****
!****   Write the parameter list to disk
!****
      lscr(1) = itype
      scr(2) = etime
      lscr(3) = lprint
      scr(4) = alpha
      scr(5) = rscale
      lscr(6) = isubwb
      lscr(7) = fno
      lscr(8) = lno
      lscr(9) = ntgp
      lscr(10) = iaddwb
      scr(11) = tadd
      scr(12) = tmult
      scr(13) = winlen
      scr(14) = v
      ngwrds(num_gain) = npars
      nglists(num_gain) = nglists(num_gain) + 1
      CALL wrdisc( igunit(num_gain), scr, ngwrds(num_gain) )
      IF( IAND(lprint,1) .NE. 0 ) THEN
          PRINT *,' itype=',itype,' etime=',etime,' alpha=',alpha,
     &      ' isubwb=',isubwb, ' fno=',fno,' lno=',lno,' ntgp=',ntgp
          PRINT *,' tadd=',tadd,' tmult=',tmult,' winlen=',winlen
      ENDIF
      IF( ntgp .GT. 0 ) THEN
          CALL wrdisc( igunit(num_gain), tgp, ntgp )
          ntgp = 0
      ENDIF
      lastlno = lno
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
      IF( token(1:nchars) .NE. 'END' .OR. nchars .NE. 3 ) GOTO 101
      RETURN
      END
