      SUBROUTINE treditex( buf, lbuf, ibuf, scr, lscr, iscr, nready )      
!*******
!*******      MAKE SURE TO USE THE RIGHT COMMON BLOCK NAME
!*******      despike and tredit 
!
!  ARGUMENTS:
!  buf   - The trace, with SEGY header as TYPE REAL
!  lbuf  - The trace, with SEGY header as TYPE INTEGER*4
!  ibuf  - The trace, with SEGY header as TYPE INTEGER*2
!
!  COPYRIGHT (C) The Regents of the University of California
!  ALL RIGHTS RESERVED.  
!  Written by Paul Henkart, SIO, August 1992
!
!  Mod 14 June 1994 - Add type and type 2 (Trehu/Sutton)
!  Mod 18 Aug 95 by Mike Holzrichter, ODP and pch - Add quart, remove type
!  Mod 1 Apr 96 - Threshold won't work when thres(1) + thres(2) = 0
!  Mod 22 July 97 - Add KILL, ALPHA, SET, ADDWB, VEL, SES, SEL
!  Mod 2 Apr 98 - Add MEDIAN
!  Mod 28 May 98 - Several things for median stack
!  Mod 1 Apr 99 - Add MINVAL/VALMIN
!  Mod 9 Sept 00 - SES/SEL was going to Trehu
!  Mod 11 Sept 00 - Add LIMITS, IHDR, LHDR, IHDR
!  Mod ????? - increase iapsiz = 5000000
!  Mod 11 Sept 01 - Make sure SES and SEL don't get too small
!                   when ADDWB is used.
!  Mod 15 Feb 02 - Add WINLEN and HCYCLE
!  Mod 15 Jul 02 - SES/SEL had bad sel end time check
!  Mod 11 Oct 02 - DESPIKE didn't work at ALL due to bad common name.
!                  Was   COMMON /despik/ instead of /despike/
!  Mod 24 Aug 05 - Add TAILMUTE
!  Mod 9 Apr 09 - Use numdat as well as unsigned int for number of samples.
!  Mod 13 Oct 09 - Add offline bandpass filtering
!  Mod 2 Dec 09 - thres & kill didn't work
!  Mod 17 May 10 - Stop using real mute times in header words 47 & 48
!  Mod 18 Jul 10 - Replace NaN with 0.
!  Mod 1 Aug 12 - Add valmax, mintype, maxtype, devpct
!  Mod 17 Jun 14 - Do the offline filter for ALL pickers
!                - Use frequency filter with dbdrop 48 rathe than time domain filter
!
!
      PARAMETER (iapsiz = 5000000)
      REAL    APDATA(0:iapsiz)
      INTEGER*2 IAPDATA(0:iapsiz)
      INTEGER LAPDATA(0:iapsiz)
      COMMON /apmem/ apdata
      EQUIVALENCE (APDATA,IAPDATA)
      PARAMETER (NPTS = 55)
      DIMENSION buf(111), lbuf(111), ibuf(111), scr(111),
     &          lscr(111), iscr(111)
      INTEGER*2 ibuf, iscr
!*******
!*******      MAKE SURE TO USE THE RIGHT COMMON BLOCK NAME
!*******      despike and tredit 
      COMMON /tredit/ thres(2), lprint, fno, lno, ftr, ltr, type, fac,
     &       quart, kill, alpha, set(4), addwb, vel, ses(2), sel(2),
     &       lunmed, valmin, limits(2), ihdr, lhdr, hdr, winlen, hcycle,
     &       tailmute, pass(2), valmax, mintype, maxtype, devpct

      REAL limits
      INTEGER fno, lno, ftr, ltr, type, kill, addwb
      DIMENSION filpts(NPTS),ldummy(2)
      COMMON /readt/ itunit, numhdr, numdat, ihunit, ireeln, jntrcs,
     *               ifmt, nskip, secs, lrenum, isrcf, idtype,
     *               nfskip, jform, itxsi, itxdel, nfktrc, norigtr
      COMMON /sioap/ iasgnd, irelse, in, iout, nextad, lapsiz, ifree,
     *               iuseap
      COMMON /segyptr/ llsegptr, lrseqptr, lshotptr, lshtrptr, lrpnptr,
     *                 lrptrptr, itridptr, ldisptr,  lwbdptr,  lsxcoptr,
     *                 lrxcoptr, idelmptr, istmptr,  iendmptr, isampptr,
     *                 isiptr,   iyrptr,   idayptr,  ihrptr,   iminptr,
     *                 isecptr,  igmtptr,  ldelsptr,  lsmusptr,lemusptr,
     *                 lsisptr,  lwbtsptr, lgatptr,  lssmsptr, lesmsptr,
     *                 lsbptr,   ifoldptr, icvleptr, lespnptr
      SAVE
      DATA iaddr/0/, ntraces/0/, msamps/100000/, last/0/, mindel/100000/
      DATA iwarn1/0/, nfilpts/0/
!
!
      IF( IAND(lprint,2) .NE. 0 ) THEN
          PRINT *,' thres=',thres,' fno=',fno,' lno=',lno,' ftr=',ftr,
     &     ' ltr=',ltr,' type=',type,' fac=',fac,' quart=',quart
          PRINT *,' kill=',kill,' alpha=',alpha,' set=',set
          PRINT *,' ses=',ses,' sel=',sel,' vel=',vel,' addwb=',addwb
          PRINT *,' limits=',limits,' ihdr=',ihdr,' lhdr=',lhdr,
     &            ' hdr=',hdr,' winlen=',winlen,' hcycle=',hcycle
          PRINT *,' endmute=',endmute,' pass=',pass,' valmax=',valmax,
     &         ' mintype=',mintype,' maxtype=',maxtype,' devpct=',devpct
      ENDIF
      nready = 1

      IF( ibuf(itridptr) .EQ. 2 .AND. lbuf(51) .NE. -1 ) RETURN         ! forget dead traces
      IF( lbuf(lrptrptr) .EQ. 0 ) THEN
          lnum = lbuf(lrpnptr)
          ltrno = lbuf(lrptrptr)
      ELSE
          lnum = lbuf(lshotptr)
          ltrno = lbuf(lshtrptr)
      ENDIF
      IF( lnum .LT. fno .OR. lnum .GT. lno ) RETURN
      IF( ltrno .LT. ftr .OR. ltrno .GT. ltr ) RETURN
!      nsamps = ibuf(isampptr)
      nsamps = numdat
      delay = buf(ldelsptr)
      idelay = NINT(delay / buf(lsisptr)) + 1
      si = buf(49)
!****
!****  Put all the data into apdata(nextad) then filter it if requested
!****
!**** force use of the "ap" so it's memory gets allocated and copied to apdata(nextad)
      CALL inap( buf(numhdr+1), numdat )
!**** get the data out of the ap into scr array apdata(nextad)
      iout = 0
      CALL rlseap( apdata(nextad), numdat )
!****  Well, rlseap set in to 0, so CALL rlseap( buf(numhdr+1), numdat )
!**** won't make a copy into buf, which median assumes.  Copy it explicitly.
      DO i = 1, numdat
         buf(numhdr+i) = apdata(nextad+i-1)
      ENDDO
!****
!****  Check for NaN, which also checks as zero.
!****
      DO i = 1, numdat
         IF( buf(numhdr+i) .EQ. 0 )  buf(numhdr+i) = 0.
      ENDDO
!****
!****  Do MEDIAN first because it's so unlike the rest because it needs
!****  a bunch of traces before it can do any work.
!**** 
      IF( lunmed .NE. 0 ) THEN
!**** Store all of the gather on disc until the end-of-gather flag.
          nwrds = nsamps + numhdr
          CALL wrdisc( lunmed, buf(1), nwrds )
          iaddr = iaddr + nwrds
          ntraces = ntraces + 1
          msamps = MIN( msamps, nsamps )
          mindel = MIN( mindel, idelay )
          last = MAX( last, idelay+nsamps )
          IF( lbuf(51) .NE. -1 ) THEN
              nready = 0
              RETURN
          ENDIF
          iaddr = 0
          CALL podisc( lunmed, 0, 0 )
          index = 1
          nosamps = last - mindel
          IF( nosamps*ntraces .GT. iapsiz ) THEN
              PRINT *,' ***  ERROR  ***  Not enough memory for MEDIAN.'
              PRINT *,' Decrease the record length.'
              PRINT *,' Make sure there is an end-of-gather flag.'
              STOP
          ENDIF
!****     Now put all of the traces in memory so they all start at the
!****     same delay and have the same length.
          DO i = 1, ntraces
             CALL rddisc( lunmed, scr, numhdr, istat)
!             nsamps = iscr(isampptr)
             CALL ushort2long( iscr(isampptr), nsamps )
             delay = scr(ldelsptr)
             idelay = NINT(delay / scr(lsisptr)) + 1
             IF( idelay .GT. mindel ) THEN
                 DO j = mindel, idelay
                    lapdata(index) = -1
                    index = index + 1
                 ENDDO
             ENDIF
             CALL rddisc( lunmed, apdata(index), nsamps, istat)
!****        fill/pad the ends with -1 if too short or muted and we'll
!****        exclude them later.
!             IF( scr(lemusptr) .NE. 0 ) THEN                            ! end mute time in seconds
!                 n = (scr(lemusptr)-delay) / scr(lsisptr) + 1
             IF( iscr(iendmptr) .NE. 0 ) THEN
                 n=(FLOAT(iscr(iendmptr))/1000.-delay)/scr(lsisptr)+1
                 DO j = 0, n-1
                    lapdata(index+j) = -1
                 ENDDO
             ENDIF
             index = index + nsamps
             n = nosamps - (idelay - mindel + nsamps)
             IF( n .GT. 0 ) THEN
                 DO j = 0, n-1
                    lapdata(index+j) = -1
                 ENDDO
                 index = index + n
             ENDIF
          ENDDO
          DO i = 1, nosamps
             n = 0
             jndex = i
             DO j = 1, ntraces
                IF( lapdata(jndex) .NE. -1 ) THEN
                    n = n + 1
                    scr(n) = apdata(jndex)
                ENDIF
                jndex = jndex + nosamps
             ENDDO
             IF( n .LE. 1 ) THEN
                 IF( n .EQ. 0 ) buf(numhdr+i) = 0.
                 IF( n .EQ. 1 ) buf(numhdr+i) = scr(1)
             ELSE
                 CALL sort( n, scr )
                 n2 = n / 2
                 IF( 2*n2 .EQ. n ) THEN
                     buf(numhdr+i) = .5 * (scr(n2)+scr(n2+1))
                 ELSE
                     buf(numhdr+i) = scr(n2+1)
                 ENDIF
             ENDIF
          ENDDO
          IF( lbuf(lrptrptr) .NE. 0 ) THEN
              lbuf(lrptrptr) = 1
          ELSE
              lbuf(lshtrptr) = 1
          ENDIF
          ibuf(isampptr) = nosamps
          CALL long2ushort( nosamps, ibuf(isampptr) )
          numdat = nosamps
          index = 1
          msamps = 100000
          ntraces = 0
          last = 0
          mindel = 100000
          CALL podisc( lunmed, 0, 0 )
          nready = 1
          RETURN
      ENDIF
!****
!****   Do all the work from scr rather than buf so we can do the
!**** detection on the trace that has been scaled.  But the interpolation
!**** or kill is done in buf.  i.e. Don't clobber the trace!
!****
      DO i = 1, nsamps
         scr(i) = buf(numhdr+i)
      ENDDO
      IF( pass(2) .NE. 0. ) THEN
          srate = 1. / si
          idbdrop = 48
          CALL filters( numdat, scr, srate, pass(1), pass(2),
     &             idbdrop, 0, 0, 0)
      ENDIF
      IF( alpha .NE. 0.) THEN
          DO i = 1, nsamps
             scr(i) = scr(i) ** alpha
          ENDDO
      ENDIF
!****
!****
      IF( set(2) .EQ. 0. ) THEN
          is = 1
          ie = nsamps
      ELSE
          time1 = set(1)
          time2 = set(2)
          time3 = set(3)
          time4 = set(4)
          IF( addwb .NE. 0 ) THEN
              time1 = time1 + buf(lwbtsptr)
              time2 = time2 + buf(lwbtsptr)
              time3 = time3 + buf(lwbtsptr)
              time4 = time4 + buf(lwbtsptr)
          ENDIF
          IF( vel .NE. 0 ) THEN
              range = FLOAT(lbuf(ldisptr))
              time1 = SQRT(time1*time1 + (range*range)/(vel*vel))
              time2 = SQRT(time2*time2 + (range*range)/(vel*vel))
              time3 = SQRT(time3*time3 + (range*range)/(vel*vel))
              time4 = SQRT(time4*time4 + (range*range)/(vel*vel))
          ENDIF
          is = NINT((time1-delay) / buf(lsisptr)) + 1
          ie = NINT((time2-delay) / buf(lsisptr)) + 1
          is2 = NINT((time3-delay) / buf(lsisptr)) + 1
          ie2 = NINT((time4-delay) / buf(lsisptr)) + 1
          IF( is .LT. 1 ) is = 1
          IF( ie .GT. nsamps ) ie = nsamps
          IF( is2 .LT. 1 ) is2 = 1
          IF( ie2 .GT. nsamps ) ie2 = nsamps
      ENDIF
!****
!****   type 1 - remove a spike if it exceeds a threshold
!****
      IF( thres(1) .NE. 0. .AND. thres(2) .NE. 0. ) THEN
          IF( IAND(lprint,2) .NE. 0 ) PRINT *,' Type 1 - Threshold'
   10     CONTINUE
          DO 100 istart = is, ie
             IF( scr(istart) .LT. thres(1) .OR. 
     &           scr(istart) .GT. thres(2) ) THEN
                 IF( tailmute .NE. 0. ) THEN
!****            zero from the beginning of the spike to the end
                     iadd = tailmute / buf(49)
                     jstart = istart + iadd
                     IF( IAND(lprint,8) .NE. 0 ) PRINT *,
     &                   ' zeroing from sample ',jstart,' to the end.'
                     DO i = jstart, nsamps
                        buf(numhdr+i) = 0.
                     ENDDO
                     RETURN
                 ENDIF
!****            If outside the 2 thresholds, find the next one between
!****            the two thresholds.
                 DO 20 j = istart + 1, ie
                    IF( scr(j) .GE. thres(1) .AND.
     &                  scr(j) .LE. thres(2) ) THEN
                        iend = j - 1
                        GOTO 110
                    ENDIF
   20            CONTINUE
!****            Nothing in the rest of the window lies between thresholds
!****            here if bad to the end of the trace.  Zero it.
                 DO j = istart, nsamps
   50               buf(numhdr+j) = 0.
                 ENDDO
                 IF( IAND(lprint,8) .NE. 0 ) PRINT *,' zeroing from ',
     &               istart,' to ',nsamps
                 RETURN
             ENDIF
  100     CONTINUE
!***      Got here if finished the whole trace and all is ok
          RETURN
  110     CONTINUE
!****     Got here if we found the begining and end of the "spikes"
          IF( kill .NE. 0 ) GOTO 9000
          IF( istart .EQ. 1 ) THEN
              DO i = 1, iend
  120            buf(numhdr + i) = 0.
              ENDDO
              IF( IAND(lprint,8).NE.0) PRINT *,' zeroing from 1 to ',iend
          ELSE
              n = iend - istart + 1
              x = ( buf(numhdr+iend+1)-buf(numhdr+istart-1) ) / FLOAT(n)
              IF( IAND(lprint,8) .NE. 0 ) PRINT *,' interpolating ',n,
     &               ' between ',istart-1,' to ',iend+1,' x=',x
              DO 150 i = istart, iend
                 buf(numhdr+i) = buf(numhdr+i-1) + x
  150         CONTINUE
          ENDIF
          is = iend + 1
          GOTO 10
      ENDIF
!****
!****    type 2 - Trehu/Sutton
!****    SES/SEL also uses fac
!****
      IF( fac .NE. 0. .AND. ses(1)+ses(2)+sel(1)+sel(2) .EQ. 0. ) THEN
          IF( IAND(lprint,2) .NE. 0 ) PRINT *,' Type 2 - Trehu/Sutton'
          DO 300 i = is+2, ie-2
             adiff1 = (ABS(scr(i-2)-scr(i-1))+ABS(scr(i+1)-scr(i+2)))/2.
             adiff2 = (ABS(scr(i-1)-scr(i)) + ABS(scr(i)-scr(i+1))) / 2.
             IF( adiff2 .GT. fac * adiff1 ) THEN
                 IF( kill .NE. 0 ) GOTO 9000
                 buf(numhdr+i) = (buf(numhdr+i+1) + buf(numhdr+i-1)) /2.
             ENDIF
  300     CONTINUE
      ENDIF
!****
!****    CWP/SU  Quartile
!****
      IF( quart .NE. 0. ) THEN
          IF( IAND(lprint,2) .NE. 0 ) PRINT *,' Type 3 - Quartile'
          k = NINT( quart * FLOAT(ie) )
          IF (K .LT. is ) K = is
          IF (K .GT. ie) K = ie
!         COMPUTE ABSOLUTE VALUES
          DO 310 I = 1, nsamps
             scr (I) = ABS (buf(numhdr+i))
  310     CONTINUE
!         FIND THE VALUE OF QUANTH SAMPLE BY PARTIALLY SORTING
          ILOW = is
          IHI  = ie
          DO WHILE (ILOW .LT. IHI)
             ak = scr (K)
             I  = ILOW
             J  = IHI
             DO WHILE (I .LE. J)
                DO WHILE (scr (I) .LT. ak)
                   I = I+1
                ENDDO
                DO WHILE (scr (J) .GT. ak)
                   J = J-1
                ENDDO
                IF (I .LE. J) THEN
                   TEMP = scr (I)
                   scr (I) = scr (J)
                   scr (J) = TEMP
                   I = I + 1
                   J = J - 1
                ENDIF
             ENDDO
             IF (J .LT. K) ILOW = I
             IF (I .GT. K) IHI  = J
          ENDDO
!         PERFORM THE CLIP
          ak = scr (K)
          DO 400 I = 1, nsamps
             IF( scr(i) .GT. ak) THEN
                 IF( kill .NE. 0 ) GOTO 9000
                 IF (buf(numhdr+i) .GT. 0.0) THEN
                    buf(numhdr+i) = ak
                 ELSE
                    buf(numhdr+i) = -ak
                 ENDIF
             ENDIF
  400    ENDDO
         RETURN
      ENDIF
!****
!****   SES/SEL   (short over long average)
!****
      IF( SES(2) .NE. 0. ) THEN
          IF( IAND(lprint,2) .NE. 0 ) PRINT *,' Type 4 - SES/SEL'
          time1 = ses(1)
          time2 = ses(2)
          time3 = sel(1)
          time4 = sel(2)
          IF( addwb .NE. 0 ) THEN
              time1 = time1 + buf(lwbtsptr)
              time2 = time2 + buf(lwbtsptr)
              time3 = time3 + buf(lwbtsptr)
              time4 = time4 + buf(lwbtsptr)
          ENDIF
          IF( vel .NE. 0 ) THEN
              range = FLOAT(lbuf(ldisptr))
              time1 = SQRT(time1*time1 + (range*range)/(vel*vel))
              time2 = SQRT(time2*time2 + (range*range)/(vel*vel))
              time3 = SQRT(time3*time3 + (range*range)/(vel*vel))
              time4 = SQRT(time4*time4 + (range*range)/(vel*vel))
          ENDIF
!****     This is stupid.  NINT((time1-delay) / buf(lsisptr))
!****     gives:    Note: IEEE floating-point exception flags raised:
!****           Inexact;  Invalid Operation;
!****     when the start/end time and the delay are less than a sample
!****     interval apart.  argh.
          IF( time1-delay .LT. buf(lsisptr) ) THEN
              is = 1
          ELSE
              is = NINT((time1-delay) / buf(lsisptr)) + 1
          ENDIF
          IF( time2-delay .LT. buf(lsisptr) ) THEN
              ie = 1
          ELSE
              ie = NINT((time2-delay) / buf(lsisptr)) + 1
          ENDIF
          IF( is .LT. 1 ) is = 1
          IF( ie .GT. nsamps ) ie = nsamps
!          print *,time1,delay,buf(lsisptr),is
!          print *,time2,delay,buf(lsisptr),ie
          IF( is .GE. ie ) THEN
              IF( iwarn1 .LT. 10 ) THEN
                  iwarn1 = iwarn1 + 1
                  PRINT *,' SES window is too short, no editing.'
              ENDIF
              RETURN
          ENDIF 
!          IF( pass(2) .NE. 0. .AND. kill .GT. 0 ) THEN
!              IF( nfilpts .EQ. 0 ) THEN
!                  CALL bpass( pass(1), pass(2), filpts, NPTS, si, scr )
!                  nfilpts = NPTS
!              ENDIF
!              ldummy(1) = 1
!              ldummy(2) = numdat
!              CALL tvfilt( apdata(nextad), scr, apdata(nextad+numdat),
!     &             filpts, nfilpts, 1., ldummy, 1, numdat )
!          ENDIF
          ave1 = 0.
          DO i = is, ie
             ave1 = ave1 + ABS(scr(i))
          ENDDO
          ave1 = ave1 / FLOAT(ie-is+1)
          IF( time3-delay .LT. buf(lsisptr) ) THEN
              is = 1
          ELSE
              is = NINT((time2-delay) / buf(lsisptr)) + 1
          ENDIF
          IF( time4-delay .LT. buf(lsisptr) ) THEN
              ie = 1
          ELSE
              ie = NINT((time4-delay) / buf(lsisptr)) + 1
          ENDIF
          IF( is .LT. 1 ) is = 1
          IF( ie .GT. nsamps ) ie = nsamps
!          print *,time3,delay,buf(lsisptr),is
!          print *,time4,delay,buf(lsisptr),ie
          IF( is .GE. ie ) THEN
              IF( iwarn1 .LT. 10 ) THEN
                  iwarn1 = iwarn1 + 1
                  PRINT *,' SEL window is too short, no editing.'
              ENDIF
              RETURN 
          ENDIF
          ave2 = 0.
          DO i = is, ie
             ave2 = ave2 + ABS(scr(i))
          ENDDO
          ave2 = ave2 / FLOAT(ie-is+1)
          IF( IAND(lprint,4) .NE. 0 ) THEN
              PRINT *,time1,time2,time3,time4
              PRINT *,' SES abs ave = ',ave1,' SEL abs ave = ',ave2
          ENDIF
          IF( ave1 / ave2 .GT. fac ) GOTO 9000
          RETURN
      ENDIF
!****
!****   Kill the trace if the average exceeds a threshold.  Only 1
!****   thres may be given.
!****
      IF( thres(1) .NE. 0. .AND. thres(2) .EQ. 0.) THEN
          IF( IAND(lprint,2) .NE. 0 ) PRINT *,' Type 5 - thres'
          ave1 = 0.
          DO i = is, ie
             ave1 = ave1 + ABS(scr(i))
          ENDDO
          ave1 = ave1 / FLOAT(ie-is+1)
          IF( IAND(lprint,4) .NE. 0 )
     &        PRINT *,' index ',is,' to ',ie,' ave= ',ave1
          IF( ave1 .GT. thres(1) ) GOTO 9000
          IF( set(3)+set(4) .NE. 0 ) THEN
              ave1 = 0.
              DO i = is2, ie2
                 ave1 = ave1 + ABS(scr(i))
              ENDDO
              ave1 = ave1 / FLOAT(ie2-is2+1)
              IF( IAND(lprint,4) .NE. 0 )
     &            PRINT *,' index',is2,' to ',ie2,' ave= ',ave1
              IF( ave1 .GT. thres(1) ) GOTO 9000
              RETURN
          ENDIF
      ENDIF
!****
!****   Detect dead traces if no value exceed VALMIN
!                 IF( token(1:6) .EQ. 'ABSVAL' ) maxtype = 1
!                 IF( token(1:4) .EQ. 'SDEV' ) maxtype = 2
!                 IF( token(1:3) .EQ. 'VAR' ) maxtype = 3
!                 IF( token(1:4) .EQ. 'SKEW' ) maxtype = 4
!                 IF( token(1:4) .EQ. 'KURT' ) maxtype =5
!****
      IF( valmin .NE. 0. .AND. mintype .EQ. 1 ) THEN
          IF( IAND(lprint,2) .NE. 0 ) PRINT *,' Type MINVAL',valmin
          DO i = is, ie
             IF( ABS(scr(i)) .GE. valmin ) RETURN
          ENDDO
          GOTO 9000
      ENDIF
      IF( valmax .NE. 0 .AND. mintype .EQ. 1 ) THEN
          IF( IAND(lprint,2) .NE. 0 ) PRINT *,' Type MAXVAL'
          big = 0.
          DO i = is, ie
             IF( ABS(scr(i)) .GE. big ) big = ABS(scr(i))
          ENDDO
          IF( big .LE. valmax ) RETURN
          GOTO 9000
      ENDIF
      IF( valmin .NE. 0. .AND. mintype .EQ. 2 ) THEN
          CALL moment( scr(is), ie-is+1, ave, sdev, var, skew, curt )
          IF( IAND(lprint,16) .NE. 0 ) PRINT *,' std deviation=',sdev
          IF( sdev .GE. valmin ) RETURN
          GOTO 9000
      ENDIF
      IF( valmax .NE. 0. .AND. maxtype .EQ. 2 ) THEN
          CALL moment( scr(is), ie-is+1, ave, sdev, var, skew, curt )
          IF( IAND(lprint,16) .NE. 0 ) PRINT *,' std deviation=',sdev
          IF( sdev .LE. valmax ) RETURN
          GOTO 9000
      ENDIF
      IF( valmin .NE. 0. .AND. mintype .EQ. 3 ) THEN
          CALL moment( scr(is), ie-is+1, ave, sdev, var, skew, curt )
          IF( IAND(lprint,16) .NE. 0 ) PRINT *,' variance=',var
          IF( var .GE. valmin ) RETURN
          GOTO 9000
      ENDIF
      IF( valmax .NE. 0. .AND. maxtype .EQ. 3 ) THEN
          CALL moment( scr(is), ie-is+1, ave, sdev, var, skew, curt )
          IF( IAND(lprint,16) .NE. 0 ) PRINT *,' variance=',var
          IF( var .LE. valmax ) RETURN
          GOTO 9000
      ENDIF
      IF( valmin .NE. 0. .AND. mintype .EQ. 4 ) THEN
          CALL moment( scr(is), ie-is+1, ave, sdev, var, skew, curt )
          IF( IAND(lprint,16) .NE. 0 ) PRINT *,' skewness=',skew
          IF( skew .GE. valmin ) RETURN
          GOTO 9000
      ENDIF
      IF( valmax .NE. 0. .AND. maxtype .EQ. 4 ) THEN
          CALL moment( scr(is), ie-is+1, ave, sdev, var, skew, curt )
          IF( IAND(lprint,16) .NE. 0 ) PRINT *,' skewness=',skew
          IF( skew .LE. valmax ) RETURN
          GOTO 9000
      ENDIF
      IF( valmin .NE. 0. .AND. mintype .EQ. 5 ) THEN
          CALL moment( scr(is), ie-is+1, ave, adev, sdev, var,skew,curt)
          IF( IAND(lprint,16) .NE. 0 ) PRINT *,' kurtosis=',curt
          IF( curt .GE. valmin ) RETURN
          GOTO 9000
      ENDIF
      IF( valmax .NE. 0. .AND. maxtype .EQ. 5 ) THEN
          CALL moment( scr(is), ie-is+1, ave, sdev, var, skew, curt )
          IF( IAND(lprint,16) .NE. 0 ) PRINT *,' kurtosis=',curt
          IF( curt .LE. valmax ) RETURN
          GOTO 9000
      ENDIF
!****
!****   Kill/keep by header value
!****         KILL 2 (INSIDE) means to kill the trace if the header
!****               value is inside the limits.
!****         KILL 3 (OUTSIDE)  means to kill the trace if the header
!****               value is outside the limits.
!****
      IF( limits(1) .NE. -999999. ) THEN
          IF( IAND(lprint,2) .NE. 0 ) PRINT *,' Type LIMITS'
          IF( ihdr .NE. 0 ) THEN
              index = ihdr
              value = ibuf(index)
          ENDIF
          IF( lhdr .NE. 0 ) value = lbuf(lhdr)
          IF( hdr .NE. 0. ) THEN
              index = NINT(hdr)
              value = buf(index)
          ENDIF
          IF( kill .EQ. 2 ) THEN
              IF( value .GE. limits(1).AND.value.LE.limits(2)) GOTO 9000
          ENDIF
          IF( kill .EQ. 3 ) THEN
              IF( value .LT. limits(1) .OR.value.GT.limits(2)) GOTO 9000
          ENDIF
      ENDIF
!****
!****  LTZ - Local Trace Zeroing
!****
      IF( winlen .NE. 0. ) THEN
          nwin = NINT(winlen/buf(lsisptr))
          nwin2 = NINT(winlen/buf(lsisptr)/2.)
          nhcycle = hcycle/buf(lsisptr)
!	 print *,nwin, nwin2, ' nhcycle=', nhcycle,' nsamps=',nsamps
          ave = 0.
          DO i = 1, nwin
             ave = ave + buf(numhdr+i)
          ENDDO
!****     find the "average" window amplitude (sum of amplitudes) for all windows
          apdata(nextad) = ave
!	 print *,' index=',(nextad),' ave=',ave
          DO i = 1, nsamps - nwin - 1
             ave = ave - buf(numhdr+i) + buf(numhdr+nwin+i)
             apdata(nextad+i) = ave
!	 print *,' index=',(nextad+i), buf(numhdr+i),
!     &  buf(numhdr+nwin+i),' ave=',apdata(nextad+i)
          ENDDO
!****     find the zero crossings of the window averages
          nz = 0
          DO i = nextad, nextad+nsamps-nwin
!	 print *,' i=',i,apdata(i) * apdata(i+1),apdata(i),apdata(i+1)
             IF( apdata(i) * apdata(i+1) .LE. 0 ) THEN
                 nz = nz + 1
                 lscr(nz) = i
!	 print *,' i=',i,' nz=',nz,' crossing at ',lscr(nz)
             ENDIF
          ENDDO
!****     zero the window if zero crossing exceeds hcycle (half cycle)
          IF( nz .GT. 0 ) THEN
              IF( lscr(1) .GT. nhcycle ) THEN
!	 print *,' zeroing from 1 to ',lscr(1)
                  DO j = 1, lscr(1)
                     buf(numhdr+j) = 0.
                  ENDDO
              ENDIF
              DO i = 1, nz-1
                 IF( lscr(i+1) - lscr(i) .GT. nhcycle ) THEN
!	 print *,' zeroing from',lscr(i),' to ',lscr(i+1)
                     DO j = lscr(i), lscr(i+1)
                        buf(numhdr+j) = 0.
                     ENDDO
                 ENDIF
              ENDDO
          ELSE
              PRINT *,' ***  Warning  ***  Trace had no zero crossings.'
          ENDIF
      ENDIF

      RETURN
!****
!****   KILL the trace
!****
 9000 CONTINUE
      IF( IAND(lprint,8) .NE. 0 ) PRINT *,' Killing shot ',lbuf(3),
     &   ' trace ', lbuf(4), ' rp ', lbuf(6), ' trace ', lbuf(7)
      IF( kill .NE. 4 ) THEN
          DO i = 1, nsamps
             buf(numhdr+i) = 0.
          ENDDO
      ENDIF
      ibuf(itridptr) = 2
      RETURN

      END
