      SUBROUTINE catex( buf, lbuf, ibuf, scr, lscr, iscr, nready )      
!
!  ARGUMENTS:
!  buf   - The trace, with SEGY header as TYPE REAL
!  lbuf  - The trace, with SEGY header as TYPE INTEGER*4
!  ibuf  - The trace, with SEGY header as TYPE INTEGER*2
!
!  COPYRIGHT (C) The Regents of the University of California
!  ALL RIGHTS RESERVED.  
!  Written by Paul Henkart, SIO, 30 September 1992
!  mod july 95 - initialize interp to 0
! 8 Apr 09 - Use common numdat rather than segy header word ISAMPPTR
! 14 Apr 09 - Wasn't working at all
!
      DIMENSION buf(111), lbuf(111), ibuf(111),
     &          scr(100), lscr(111), iscr(111)
      INTEGER*2 ibuf, iscr
      COMMON /cat/ lun, nlists, nwrds
      COMMON /edits/ ierror, iwarn, irun, now, icompt, isite, maxsamps
      INTEGER fno, ftr, fno2, ftr2
      CHARACTER*80 token
      SAVE
      COMMON /readt/ itunit, numhdr, numdat, ihunit, ireeln, jntrcs,
     *               ifmt, nskip, secs, lrenum, isrcf, idtype,
     *               nfskip, jform, itxsi, itxdel, nfktrc, norigtr
      COMMON /segyptr/ llsegptr, lrseqptr, lshotptr, lshtrptr, lrpnptr,
     *                 lrptrptr, itridptr, ldisptr,  lwbdptr,  lsxcoptr,
     *                 lrxcoptr, idelmptr, istmptr,  iendmptr, isampptr,
     *                 isiptr,   iyrptr,   idayptr,  ihrptr,   iminptr,
     *                 isecptr,  igmtptr,  ldelsptr,  lsmusptr,lemusptr,
     *                 lsisptr,  lwbtsptr, lgatptr,  lssmsptr, lesmsptr,
     *                 lsbptr,   ifoldptr, icvleptr, lespnptr
      DATA mlists/0/, fno2/0/, ndone/0/, lastno/0/, ntrcs/0/, mdone/0/,
     &     interp/0/, lunscr/0/
!
!****
!****  nready 0 means the input trace is not modified and should be passed to
!****         the next process.
!****  nready < 0 means nothing is output and we need another input trace.
!****
      nready = 0
!****
!****    Get the first parameter list
!****
   10 IF( mlists .EQ. 0 ) THEN
          CALL podisc( lun, 1, 0 )                                      ! get the first parameter list from disk
          CALL rddisc( lun, scr, nwrds, istat )
          mlists = mlists + 1
          itype = scr(1)
          fno = lscr(2)
          lno = lscr(3)
          ftr = lscr(4)
          ltr = lscr(5)
          lprint = lscr(6)
          n = lscr(7)
      ENDIF
      IF( lbuf(lrptrptr) .NE. 0 ) THEN
          lnum = lbuf(lrpnptr)
          ltrno = lbuf(lrptrptr)
      ELSE
          lnum = lbuf(lshotptr)
          ltrno = lbuf(lshtrptr)
      ENDIF
      IF( IAND(lprint,2) .NE. 0 ) PRINT *,' lnum=',lnum,' ltrno=',ltrno,
     &    ' fno=',fno,' lno=',lno
      IF( lnum .LT. fno .AND. mlists .EQ. 1 ) RETURN
      IF( lnum .GT. lno .AND. mlists .GE. nlists ) RETURN
      IF( ltrno .LT. ftr .OR. ltrno .GT. ltr ) RETURN
!****
!****
!****
   20 IF( IAND(lprint,2) .NE. 0 ) THEN
          PRINT *,' itype=',itype,' fno=',fno,' lno=',lno,' ftr=',ftr,
     &            ' ltr=',ltr,' n=',n
      ENDIF
      IF( lnum .LT. fno ) THEN                                          ! assume shot numbers always increase
          IF( mlists .EQ. 1 ) THEN
              IF( interp .EQ. 0 ) RETURN
              GOTO 500
          ENDIF
      ENDIF
      IF( lnum .LE. lno ) GOTO 500                                      ! is the shot in this list
      IF( mlists .GT. nlists ) THEN
          IF( interp .EQ. 0 ) RETURN
          GOTO 500
      ENDIF
      IF( lnum .GE. fno2 ) THEN
          IF( mlists .GT. 1 ) THEN
              itype = itype2
              fno = fno2
              lno = lno2
              ftr = ftr2
              ltr = ltr2
              n = n2
          ENDIF
          IF( mlists .GE. nlists ) RETURN
          CALL rddisc( lun, lscr, nwrds, istat )                        ! read the next parameter list
          itype2 = scr(1)
          fno2 = lscr(2)
          lno2 = lscr(3)
          ftr2 = lscr(4)
          ltr2 = lscr(5)
          n2 = lscr(7)
          mlists = mlists + 1
          GOTO 20
      ENDIF
!****
!****   Do it 
!****
  500 CONTINUE
      nsamps = numdat
      nready = -1
      IF( ndone .EQ. 0 .AND. ntrcs .EQ. 0 ) THEN
          	IF( lunscr .NE. 0 ) CALL frefil( 3, lunscr, istat )
          CALL getfil( 1, lunscr, token, istat )
      ENDIF
!****
!****   Shot/rp concatenation here.  Save all the traces on disk first.
!****  This gets messy if we don't have the correct number of input
!****  traces per shot in common or the same number of traces in each rp
!****
      IF( itype .EQ. 1 ) THEN
          CALL wrdisc( lunscr, buf, nsamps+numhdr )                     ! write the header too.
          IF( lnum .NE. lastno ) THEN
              ndone = ndone + 1
              ntrcs = 0
          ENDIF
          ntrcs = ntrcs + 1
          lastno = lnum
          IF( ndone .EQ. n .AND. ntrcs .EQ. jntrcs ) nready = ntrcs
          mdone = 0
          RETURN
      ENDIF
!****
!****   Trace concatenation here
!****
      IF( itype .EQ. 2 ) THEN
          CALL wrdisc( lunscr, buf, nsamps+numhdr )                     ! write the header too.
          ndone = ndone + 1
          IF( ndone .EQ. n ) nready = 1
          RETURN
      ENDIF
!*****
!*****
!*****
!*****

      ENTRY getcat( buf, lbuf, ibuf, scr, lscr, iscr )
!****
!****   Return a shot/rp concatenated trace.  Search the scratch file
!****  for all the traces with the same trace number.  This assumes
!****  That the traces start with 1 and are strictly increasing by 1.
!****
      IF( itype .EQ. 1 ) THEN
          mdone = mdone + 1
          nsamps = 0
          ndone = 0
          ntrcs = 0
          CALL podisc( lunscr, 1, 0 )
  600     CONTINUE
          CALL rddisc( lunscr, scr, numhdr, istat )
          IF( istat .NE. numhdr ) THEN
              nsamps = numdat
              GOTO 800
          ENDIF
          IF( lscr(lrptrptr) .NE. 0 ) THEN
              ltrno = lscr(lrptrptr)
          ELSE
              ltrno = lscr(lshtrptr)
          ENDIF
!          itemp = iscr(isampptr)                                        ! number of samples in the trace
          CALL ushort2long( iscr(isampptr), ltemp )
          lbuf(lgatptr) = lscr(lgatptr)                                 ! use the end of gather flag
          IF( ltrno .NE. mdone ) THEN
              CALL rddisc( lunscr, scr, ltemp, istat )
              GOTO 600
          ENDIF
!****     nsamps is 0 on the first of the cats
          IF( nsamps .EQ. 0 ) THEN
              DO i = 1, numhdr
                 buf(i) = scr(i)
              ENDDO
          ENDIF
          CALL rddisc( lunscr, buf(numhdr+nsamps+1), ltemp, istat )
          nsamps = nsamps + ltemp
          numdat = nsamps
          GOTO 600
      ENDIF
!****
!****
      IF( itype .EQ. 2 ) THEN
          CALL podisc( lunscr, 1, 0 )
          CALL rddisc( lunscr, buf, numhdr, istat )
!          nsamps = ibuf(isampptr)
          CALL ushort2long( ibuf(isampptr), nsamps )
          CALL rddisc( lunscr, buf(numhdr+1), nsamps, istat )
          ndone = ndone - 1
          IF( ndone .GT. 0 ) THEN
              DO 700 i = 1, ndone
                 CALL rddisc( lunscr, scr, numhdr, istat )
!                 itemp = ibuf(isampptr)
                 CALL ushort2long( iscr(isampptr), ltemp )
                 CALL rddisc( lunscr, buf(numhdr+nsamps+1), ltemp,istat)
                 nsamps = nsamps + ltemp
  700         CONTINUE
          ENDIF
          CALL podisc( lunscr, 1, 0 )
          ndone = 0
      ENDIF
!****
!****  Return a trace in buf
!****
  800 CONTINUE
!      ibuf(isampptr) = nsamps
!      IF( nsamps .GT. 32767 ) THEN
!          ibuf(isampptr) = 32767
!          lbuf(isampptr) = nsamps
!      ENDIF
      numdat = nsamps
      CALL long2ushort( numdat, ibuf(isampptr) )
      IF( numdat .GT. maxsamps ) THEN
          PRINT *,' ***  WARNING  ***  The output number of samples ',
     &      numdat,' exceeds the SIOSEIS allocation of ', maxsamps
      ENDIF
      RETURN
      END
