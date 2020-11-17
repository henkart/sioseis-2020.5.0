      SUBROUTINE sadded
!
!                       PROCESS SADD
!                       ------- ----
!
!  Document date: 7 August 1992
!
!      Process SADD performs a scalar addition to the seismic traces, 
!  i.e. a constant number specified by the parameter SCALAR is added to
!  every amplitude.
!      Only those shot/rps and traces specifically given have the scalar
!  addition applied.
!      There is no spatial interpolation.
!
!  PARAMETER DICTIONARY
!  --------- ----------
!  SCALAR - The scalar to add to the trace.
!           Preset = 0.            e.g.   scalar 10.e3
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
      PARAMETER ( npars = 6 )                                           ! the number of user parameters
      CHARACTER*80 token
      CHARACTER*6 names(npars)
      DIMENSION scr(npars), lscr(npars)
      EQUIVALENCE (scr(1),lscr(1))
      COMMON /edits/ ierror, iwarn, irun, now, icompt
      COMMON /sadd/ lun, nlists, nwrds
      INTEGER fno, lno, ftr, ltr
      DATA names / 'SCALAR', 'LPRINT', 'FNO   ', 'LNO   ', 'FTR   ',
     &             'LTR   ' /
!**** 
!****    Set the parameter presets and various variable presets
!****
      scalar = 0.
      lprint = 0
      fno = 0
      lno = 9999999
      ftr = 0
      ltr = 999999
      nwrds = npars
      nlists = 0
      CALL getfil( 1, lun, token, istat )                               ! get a file for the SADD parameters
!****
!****     get the user's parameters
!****
      ntokes = 0                                                        ! count the tokens
  100 CONTINUE
      CALL getoke( token, nchars )                                      ! get a token and it's length
      CALL upcase( token, nchars )                                      ! convert parameter names to upper case
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
             CALL dcode( token, nchars, areal, istat )                  ! convert the alpha number to an internal machine number
             IF( istat .NE. 2 ) ierror = ierror + 1                     ! was the an error decoding it?
             IF( names(nparam) .EQ. 'SCALAR' ) scalar = areal
             IF( names(nparam) .EQ. 'LPRINT' ) lprint = NINT(areal)
             IF( names(nparam) .EQ. 'FNO' ) fno = NINT(areal)
             IF( names(nparam) .EQ. 'LNO' ) lno = NINT(areal)
             IF( names(nparam) .EQ. 'FTR' ) ftr = NINT(areal)
             IF( names(nparam) .EQ. 'LTR' ) ltr = NINT(areal)
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
!****
!****   Write the parameter list to disk
!****
      scr(1) = scalar
      lscr(2) = fno
      lscr(3) = lno
      lscr(4) = ftr
      lscr(5) = ltr
      lscr(6) = lprint
      nlists = nlists + 1
      CALL wrdisc( lun, scr, nwrds )
      IF( IAND(lprint,1) .NE. 0 ) THEN
          PRINT *,' scalar=',scalar,' fno=',fno,' lno=',lno,' ftr=',ftr,
     &            ' ltr=',ltr
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
