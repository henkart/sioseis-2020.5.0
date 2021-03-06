      SUBROUTINE diex( buf, lbuf, ibuf, scr, lscr, iscr, istop )
!  diex is the execution phase of SIOSEIS process DISKIN.  died read the
!  users parameters, wrote them to disk, and set up common.
!     This reads "pure" SEGY files, so all the machine idiosyncrasies
!  must be handled, such as DEC byte swaps and Cray not having INTEGER*2 or
!  INTEGER*4.
!     Set istop to 0 if this is not the last trace.
!     Set istop to 1 if this is the last trace.
!     Set istop to -1 if past the last trace. (ie there is no trace in buf))
!     The Cray requires some real skulduggery because it does not have 2 or 4
!  byte integers, to make matters worse, the diskio is done in 8 byte words.
!  Variable idarn is used to indicate when the data part of the trace is not
!  an even multiple of 2.  When this occurs, we must drop the last data sample
!  of the trace so the next disk read will contain the first word of the trace
!  header and we can then shift the header after unpacking it!
!
! jform = 1 means SEGY
!       = 2       SSC
!       = 3       IRIS
!       = 4       LDGO DSS-240 (pre 1992 Ewing)
!       = 5       SWAPPED
!       = 6       SU, NOHEAD
!       = 7       ODEC
!       = 8       KNUDSEN
!       = 9       EDGETECH
!       = 10      UTIG-OBS
!       = 11      EDGETECH5
!       = 12      WAV
!       = 13      ASCII
!       = 14      BINARY
!       = 15      NIUST

!
!  Copyright (C) The Regents of the University of California
!  Written by Paul Henkart, Scripps Oceanography, La Jolla, Ca 92093 Feb 1988
!  ALL RIGHTS RESERVED.
!
!  MODIFICATIONS: 
!  13 Jun 20 - Use unsigned integer from binary header for position before read on random access
!  3 Jun 20 - trinc preset must have changed and caused error with sort
!  18 Jul 18 - Convert ODEC nav to arcseconds from floating point
!  5 Nov 15 - Added FORMAT NIUST
!  11 Apr 15 - Add iaddr64 for 64 bit disc address for realtime and random
!  8 Jul 14  - Clean up shot/rp searches when fno/lno are given, especially when
!                the shot numbers are not increasing.
!            - Add a "suggestion" when nothing is found.
!            - Rev 1 random access didn't work on MCS (ntraces > 1)
!  June 14   - Use subroutine swp_tr for trace header byte swaps so it's consistent with disko
!  31 Mar 14 - Add error if trinc is used when process SORT is used.
!            - add si to common/readt/ - will eventually remove all ref to buf(49)
!            - add segyrev to common/readt/
!            - Allow Rev 1 scalar for delay
!            - set ibuf(108) to 0 if rev 0 so that disko doesn't pass it with rev 1
!            - make sure the fixed file flag is set to 0 on rev 0
!  6 Feb 14 - Do the same byte swaps as DISKOX
!           - When converting from segy rev 0 to rev 1, be careful of the delay scalar.
!  19 Apr 12 - multi list didn't work when diex over-rode random preset of 1
!  14 Feb 12 - Set itrno to 1 when setting (4) = 1 when shot tr and rp tr = 0 (see 8 Sept 11)
!  8 Dec 11 - When realtime, change the order of getting the disk address relative 
!             to the read on the binary part.
!  28 Nov 11 - When realtime, change the order of getting the disk address relative
!              to the read on the trace header.
!  8 Sept 11 - set the shot trace number to 1 when it's 0 along with the rp num and rp tr no.
!  1 Sept 11 - Redo message when Rev 1 and fixed and shot not in file.
!  1 Sept 11 - Set the user parameter SORT to SHOT when binary header says it is and buf(7) > 0
!  12 Jul 11 - User can override the sort (shot vs cdp) when doing random disk I/O
!  6 June 11 - Change FORMAT KNUDSEN so that their segy "unspecified" get swapped right
!  27 Apr 11 - FORMAT ASCII didn't set the REAL si (buf(49))
!  21 Mar 11 - Add LNTODO, read the last N traces of the file.
!  15 Feb 11 - zero end of trace when secs > numdat
!  6 Dec 10 - 2GB file on PPC failed the random check
!  27 May 10 - Add FORMAT ASCII and BINARY.
!  24 May 10 - comment out setting water depth/time in words 50 & 54
!  17 May 10 - Comment out saving real mute times in header words 47 & 48
!  19 Nov 09 - Save the data domain (time, depth, FK, f, etc) on the first trace
!              and reset it every other time.
!  25 Aug 09 - Use new podisc_un with unsigned disk address for random
!  25 Aug 09 - Turn off warning message when random and fixed length flag = 0
!  28 Jul 09 - Print message if random is set and first trace wanted is > filesize
!  23 Jul 09 - Turn random off if file is > 2GB since random uses disk addresses
!  17 Jul 09 - Byte swap trace header words 46-60 rather than 46-53
!  19 Jun 09 - Add RANDOM
!  17 Mar 09 - Require NTRCS when NOINC is used on shot data.'
!  9 Feb 09 - Set ODEC source xy coordinate type (word 45) to 3 (decimal degs)
!  5 Feb 09 - Allow/honor IPATH of 200 characters
!  11 Jan 09 - Die gracefully with message if too many samples for sioseis
!  11 Dec 08 - Allow WAV files.
!  26 Sep 08 - Use trace header word 58 (number of samples) with c unsigned short
!            - Eliminate Cray code.
!  18 Jul 08 - SET set numdat incorrectly
!  3 June 08 - Try and detect and correct an old bug in segdin and segddin
!              in converting degrees to arcseconds.
!  29 Apr 08 - Use podisc (lun, 1, 0 ) to rewind since (lun, 0, 0 ) is wrong!
!  20 Sept 07 - Parameter NSAMPS wasn't working.
!  10 Sept 07 - fday not preset BEFORE byte swap on byte swapped machines.
!  18 July 07 - Still confusion about when to swap the binary header!
!  18 June 07 - Check for zero'd trace header and skip it if it is.
!  14 Feb 07 - format swapped swapped the headers on little endian machines
!  7 Feb 07 - Change print statement: ' Bad delay or set. delay=',delay,' set=',set
!  5 Dec 06 - Intel and format swapped didn't swap enough header bytes
!  4 Oct 06 - ODEC non-Linux needs to byte swap year/day/hr/min/sec
!           - ODEC writes bad buffers sometimes
!  mod 27 Sep 06 - Linux ODEC should NOT byte swap the binary header
!  1 Aug 06 - iday > lday crashed rather than stopping
!  27 Apr 06 - gfortran chokes on internal reads that were necessary for VMS
!  22 Sep 05 - Put path in common so others know the file name
!  1 Sep 05 - Ask the operator for C or Q when realtime times out.
!  8 Aug 05 - Wait rather than stop if REALTIME was set.
!  6 Aug. 05 - Set FDAY to the header if fgmt was given and fday wasn't.
!  14 Sep 04 - Save trace 0 on luntr0 if it happens.
!  11 Aug 04 - Save Edgetech course and speed in segy short words 48 & 49
!  8 Jul 04 - Refine fgmt/lgmt/fday/lday on day boundaries.
!  14 Jun 04 - Don't use renum & retrac from list 2 on (use from list 1 only)
!  29 Jan 04 - Knudsen (jform=8) was doing a byte swap
! 16 July 03 - Add jform 10, UTIG-OBS, for the 32 bit integer delay
!  2 June 03 - Do not save the SEG-Y headers in diex because died did
!              and headed.f may have modified them.
!  mod 28 May 03 - Tighten the check for Rev 1 so that crap in the
!                  binary header doesn't trigger the extended headers.
!  1 May 03 - Get by the SEG-Y Rev 1 Textual Header Extension records.
!  1 Apr 03 - If the sample interval is 0, get it from the binary header
!  21 Mar 03 - Create the LDEO header if there's a trace 0
!  7 Jan 03 - Edgetech on PC needs byte swap on NMEA string before moving!
!  22 Nov 02 - Allow delays > 20 secs (was trigger that delay was in samples)
!            - Add the UTIG 32 bit integer delay scheme.
!  19 Nov 02 - Set intrcs if EDGETECH was given.
!  15 Nov 02 - The 18 Oct change didn't do byte swap on binary header!
!  18 Oct 02 - Honor the new segy headers when opath changes
!  09 Aug 02 - New trace 28 choked, so create jform 4
!  03 May 02 - More HP compiler work arounds.
!  01 May 02 - The Edgetech offset (range - lbuf(10)) isn't legal.
!  11 Jul 01 - Add si 41 being 1./ 24000. and 83 being 1. / 12000. for Edgetech
!  21 Jun 01 - Need long integer in argument to rddisc when reading trace 0
!            - ldgo_tr0 was changed in sept 00 for REAL*8 args.
!  30 May 01 - do NTRGAT, then RENUM, then RETRAC
!  18 Feb 01 - Add/change some of the header bytes to swap on PC/DEC
!  23 Jan 01 - Clean up ODEC (jform = 7).
!  Mod ????? - increase iapsiz = 5000000
!  22 Dec 00 - Make the sort "stack" if the rp number, rp trace number AND
!              the "cdp fold" (short integer word 17) is non-zero
! 26 Oct 00 - DEC/INTEL need to swap SEG-Y words 20 and 22 also
! 2 Oct 00 - fno/lno and sort didn't work right
! 20 Jul 00 - Skip ldgo_tr0 if section 11 is bad
! 20 June 00 - Take care of sample intervals of 1/4096 and 1/8292
!  8 June 00 - The SEG-Y headers were being reread from the data file
!              rather than from the header file, which caused byte
!              swap problems.
!  5 June 00 - Add parameter NSHOTS
! 15 Feb 00 -  Old trace 0 trace length bug workaround rejected new correct length.
! 29 Oct 99 - The 14 July change for negative delay was bad because
!             Geometrics uses negative delays legally.
! 29 Sep 99 - Add ALT (alternate between IPATH and IPATH2), alt 2 means alternate
!             by trace.  alt 1 means alternate by record.
!  22 Sep 99 - Add EdgeTech just to make things easier on the user.
! 14 July 99 - Knudson uses an unsigned 16 bit integer in the delay and
!              does it in nanoseconds.
! 7 June 99 - lprint 8 to print the disk address before each trace header
! 12 Jan 99 - set si to 1/15000 if number of microsecs is 67
! 12 Jan 99 - set si to 1/16000 if number of microsecs is 63
! 10 Dec 98 - Sample interval of 31 really 3.125E-05
! 31 Aug 98 - Add parameter nsamps for overriding the SEG-Y trace header
! 26 Aug 98 - Create format ODEC, so format SWAPPED is true
! 19 June 98 - The 11 Mar, seg-y header removal wasn't done quite
!              right, so the "new" binary header was staring at 0!
! 11 Mar 98 - Remove writing SEG-Y headers to disk - diskin edit did it!
! 11 Mar 97 - Add kludge for 3.5, add 160 samples when format is SWAPPED
! 9 Mar 97  - Add data format "SWAPPED" (jform=5)
!           - Convert the LDEO lat/long to seconds of arc as per SEG-Y
! 25 Feb 97 - Save the LDEO trace 0 in luntr0 so geom can get it.
! 17 Oct 96 - Get around HP compiler bug on 16 bit integer SEG-Y
! 31 Aug 96 - fgmt/lgmt didn't work with giving ltr!
! 22 May 96 - Change trace id of 100 to trace id 1.
! 19 Apr 96 - Allow trace id of 100 - It means GPR.
! 18 Mar 96 - Add alltr, allno, notype, noindex, itrtype, itrindex
! 27 Sep 94 - Unset variable vmsvax changed to 4.
! 15 Jun 94 - forgat and fno/lno together didn't work.
!  9 Mar 94 - Do renum BEFORE ntrgat, rather than AFTER.
!  9 Jan 94 - Allow noinc 99999 and trinc 99999
! 27 Aug 93 -  set si to 1/2048 if number of microsecs is 488
! 25 Aug 93  - Use Iris convention of long integer word 58 if short
!              word 58 (number of samples per trace) is 32767
!  9 Aug 93  - Fix some fuzzy logic when ftr, ltr, trinc are given
!  21 May 93 - Argh. Some 16 bit data had an odd number of samples which
!              screwed up calling rddisc with 4 byte word count. Changed
!              to use a byte count - VMS diskio doesn't have byte I/O
!  19 May 93 - Argh.  Force the LDGO DSS-240 trace 0 trace length to be
!              (2048 - 240) / 4 bytes long, regardless of what the
!              header. (on the field tapes it is 0)
!  16 Feb 93 - Add the LDGO DSS-240 logic.  Convert trace 0 and then
!              throw it away
!            - Throw DSS-240 shot 0 away (The DSS-240 is screwed up).
!  30 Nov 92 - Don't change jntrcs as 23 OCT change unless jntrcs = 0
!              As per Alistair's problem.
!  10 Nov 92 - trinc didn't work except when fno and lno were used.
!  23 Oct 92 - Renum didn't work when ftr/ltr given
!            - set jntrcs to (ltr-ftr)/trinc +1 when ftr & ltr given.
!  7 July 1992 - Add stack as a sort type and allow renum on stack
!  12 June 92 - Pad the data when set is bigger than the data.
!  24 Jan 92 - Allow Iris 16 bit, 32 bit, and host real formats.
!              If Iris data has < 32767 samples, make it into SEGY.
!  17 Jan 92 - REWIND parameter didn't always work
!  6 Dec 91 - VMS mods for sort (vms diskio can't write characters)
!  26 Nov 91 - set the dead trace flag when setting istop = -1
!  27 Oct 91 - add the "end-of-sort" flag stuff.
!  29 July 91 - Add spath - read using order of the sort file.
!  21 Jan 91 - The DSS 240 trace data has a trace of 0, a trace id of 28, and
!              nsamps of 0 (which should be 964)
!  9 Jan 91 - ipath2 didn't work
!  5 Nov 90 - The delay parameter didn't work
!  3 May 90 - Add the IRIS PASSCAL pseudo SEGY format (jform=3)
!  23 Mar 90 - Add the UTIG 16 bit floating point format in place of the
!             SEGY 32 bit gain data (type 4)
! 20 Mar 90 - Allow a trace id (SEGY TRACE HEADER) of 9 to mean seismic data
!           - Also goto read another trace rather than another parameter
!             list if it's not seismic data (GOTO 200 rather than GOTO 100)
! 1 Mar 90 - Make SSC files (jform=2) IEEE rather that IBM (ifmt 5 not ifmt 1)
! 1 Feb 90 - VMS files might be zero filled to be a multiple of 512
!          - VMS also has trouble opening a new input file using the
!            same unit number.
! 10 Jan 90 - Honor secinc (fgmt and fsec must be given!!)
! 7 Jan 90 - Allow the data part of the trace to be 1 sample short (The Cray has
!            a problem when there are an odd number of samples and there is only
!            one trace in the file).
! 16 Nov 89 - the sort parameter didn't work
! 14 Nov 89 - added espn to SEGY header word 5 when SSC 
! 6 Nov 89 - noinc on rp sorted data used jntrcs
! 5 Nov 89 - Use jntrcs rather than intrcs when checking for ltr
! 23 Oct 89 - a) toss out traces when ibuf(15) isn't 0, 1, 2, 3, 6
!             b) use 2000 samples when ssc header is zero
!                and ntrtot = 99999 if SSC header says 0
!             c) correct SEGY binary header sample interval when SSC input
!                (see 29 Aug)
! 12 Oct 89 - a) correct bad SSC forgat calculation.(see 10 Jul update)
!             b) set numdat and ifmt in common when SSC
! 19 Sep 89 - Honor fsec and lsec
! 31 Aug 89 - Honor lgmt even if lday isn't given!
! 30 Aug 89 - Extend the check for zero SEGY trace header meaning EOF to
!             machines besides VMS! So non VMS SIOSEIS can read VMS files.
! 29 Aug 89 - Add the sample interval (in mics) and nsamps to the
!    pch      binary header when SSC input (Vista requires it!)
! 28 Aug 89 - 1) added forgat .GT. 1 logic
! 29 Jul 89 - 1) changed rewind header length from 850 to 900 words - pch
! 10 Jul 89 - 1) Added logic to honor fno/lno/noinc when forgat 1 (foreign gathers)
!    pch      2)  Added logic for SSC forgat
! 20 Jun 89 - Modified for "new" SSC IDF files
!    pch
! 30 Apr 89 - Fixed fgmt/lgmt/igmt (igmt=ihour*100+imin) rather than ihour*60)
!    pch
! 13 Apr 89 - Changed the arguments of mrgfk from iptype to idtype.
!    pch      iptype controls the way the user is reading data
!             (e.g. gmt vs shot vs rp vs all)
!             idtype is the domain of the data (time, frequency, fk,...)
! 12 Apr 89 - Changed the logics of the new sort parameter because Graham
!    pch      had a file with a bizarre binary header.
! 11 Apr 89 - Had the SSC delay and distance swapped.
!
      INTEGER prime, vaxunix, apollo, vaxvms, cray, convex, ieee
      PARAMETER ( prime = 1, vaxunix = 2, apollo = 3, vaxvms = 4,
     &            cray = 5, convex = 6, ieee = 7 )
      PARAMETER ( NIUST = 15 )
      DIMENSION buf(1000), scr(1000)
      INTEGER lbuf(1000), lscr(1000)
      INTEGER*2 ibuf(1000), iscr(1000), rshift
      INTEGER*2 i15, i32767
      INTEGER ldgodss(6)
      INTEGER*8 filesize, iaddr64, ltemp64, l2temp64, l3temp64
      REAL*8 ship_lat, ship_long, tail_lat, tail_long
      INTEGER fno, ftr, fday, fgmt, gmtinc, fsec, secinc, renum, decimf,
     *        trinc, rewind, allno, alltr, notype, noindex, itrtype,
     &        itrindex, random, addr_tr1
      COMMON /sioap/ iasgnd,irelse,in,iout,nextad,lapsiz,ifree,iuseap
      COMMON /segdin/ idummy(3), luntr0, luntr0a, lunldeo
      CHARACTER*210 sect11
      CHARACTER*210 token
      EQUIVALENCE (token,sec11)
      INTEGER*2 isect11(77)
      EQUIVALENCE (isect11,sect11)
      COMMON /edits/ ierror, iwarn, irun, now, icompt, isite, maxsamps,
     & nbperw, ireal
      COMMON /diskin/ ipunit, nlists, npwrds, lun2
      COMMON /diskin1/ spath, path
      CHARACTER*200 path, oldpat
      CHARACTER*80 spath
      INTEGER*2 ishort(2)
      EQUIVALENCE (ishort(1),llong)
      COMMON /readt/ itunit, numhdr, numdat, ihunit, ireeln, jntrcs,
     *               ifmt, nskip, secs, lrenum, isrcf, idtype,
     *               nfskip, jform, itxsi, itxdel, nfktrc, norigtr,
     *               nrskip, nfiles, irewind, delay, segyrev, si
      COMMON /binhdr/ ibinhdr(200)
      COMMON /sort1/ sdummy1(3), lkey1, sdummy2(3), lkey2, sdummy3(6)
      INTEGER*2 ibinhdr
      PARAMETER (iapsiz = 5000000)
      REAL    APDATA(0:iapsiz)
      COMMON /apmem/ apdata
      DIMENSION set(2), newsort(2)
      EQUIVALENCE (laddress,newsort(1)), (iflag,newsort(2))
!**** the following 2 are illegal ANSII standard because you can't
!**** equivalence an integer to a character.  VMS needs to do it in
!**** order to write a character string to disk using my diskio.
      CHARACTER*80 cheadr(40)
!****                   REMOVE THE EQUIVALENCE ON THE CRAY
      EQUIVALENCE( cheadr(1), ivms)
!      include Alistair's fk include stuff  include ',rggbl.inc'
! global include file for mrgFK to be added to diex.ftn
      COMMON /SIOLN1/ CBUF
      CHARACTER*200 CBUF
      COMMON /sioln2/ ICHAR, NCHARS, iprint, lunpo
      INTEGER retrac
      CHARACTER*5 sort
      integer TRNSPARM,MRGINIT,MRGTRC
      parameter (TRNSPARM = 1)
      parameter (MRGINIT  = 2)
      parameter (MRGTRC   = 3)
      integer IPTFKU                                                    ! Special data type for FK User type
      parameter (IPTFKU = 4)
!.. Bug fix for old disk file
      integer CHCKBIN
      parameter(CHCKBIN = 5)
      SAVE 
      DATA mlists/0/, oldpat/' '/, mtrgat/0/, idarn/0/, nhread/0/
      DATA ncdp/0/, ndead/0/, sort/'shot'/, numcdp/1/, iforcnt/0/
      DATA ntraces/0/, mtraces/0/, lunsort/0/, ibad/0/, lastshot/0/
      DATA nsampss/0/, ndone/0/, iprint1/0/, i15/15/, i32767/32767/
      DATA random/-1/
!
      in = 0                                                            ! assume the data will NOT go into the ap
      icount = 0
!****  (if it does, set in = 1 when it goes into the ap)
      IF( mlists .NE. 0 ) GOTO 200
      idomain = ibinhdr(31)
      CALL podisc( ipunit, 1, 0 )                                       ! rewind the parameter file
      CALL getfil( 2, idunit, path, istat )                             ! get a disk file unit number
  100 CALL rddiscb( ipunit, path, 200, istat )
      CALL rddisc( ipunit, scr(26), npwrds-25, istat )
      IF( istat .NE. npwrds-25 ) THEN
          PRINT *,' ***  ERROR  ***  Program error in subroutine diex.',
     *     npwrds, istat
          STOP
      ENDIF
      mlists = mlists + 1
!      WRITE( path, '(25A4)' ) (lscr(i),i=1,25)
      lprint = lscr(26)
      fno = lscr(27)
      nextno = fno
      lno = lscr(28)
      noinc = lscr(29)
      ftr = lscr(30)
      nexttr = ftr
      ltr = lscr(31)
      fday = lscr(32)
      lday = lscr(33)
      fgmt = lscr(34)
      lgmt = lscr(35)
      gmtinc = lscr(36)
      fsec = lscr(37)
      lsec = lscr(38)
      secinc = lscr(39)
      IF( mlists .EQ. 1 ) renum = lscr(40)
      secs = scr(41)
      jfmt = lscr(42)
      sinew = scr(43)
      decimf = lscr(44)
      delayj = scr(45)
      intrcs = lscr(46)
      ntrgat = lscr(47)
      trinc = lscr(48)
      ialt = lscr(49)
      set(1) = scr(50)
      set(2) = scr(51)
      rewind = lscr(52)
      jform = lscr(53)
      forgat = scr(54)
      mintrs = lscr(55)
      IF( mlists .EQ. 1 ) THEN
          retrac = lscr(56)
          lretrac = retrac
      ENDIF
      jsort = lscr(57)
      ascii = scr(58)
      allno = lscr(60)
      alltr = lscr(61)
      notype = lscr(62)
      noindex = lscr(63)
      itrtype = lscr(64)
      itrindex = lscr(65)
      nsamps = lscr(66)
      ntodo = lscr(67)
      IF( lscr(68) .EQ. -1 ) THEN
!****     the user did not give random, did diex set it before?
          IF( random .EQ. -1 ) random = lscr(68)
      ELSE
!***      let the user override and program decision
          random = lscr(68)
      ENDIF
      IF( random .NE. 0 .AND. ftr .LE. 0 .AND. fno+lno .NE. 0 ) ftr = 1
      lntodo = lscr(69)
      IF( spath .NE. ' ' ) THEN
          IF( icompt .EQ. vaxvms ) CALL frefil( -2, lunsort, istat )
          CALL getfil(4, lunsort, spath, istat )
          IF( istat .NE. 0 ) THEN
               PRINT *,' ***  ERROR  *** Could not open file', spath
               STOP
          ENDIF
          IF( trinc .NE. 99999 ) THEN
              IF( lkey1 .EQ. 4 .OR. lkey2 .EQ. 4 ) THEN
                  PRINT *,' *** ERROR  ***  DISKIN parameter TRINC does 
     &not work when SORT is used.'
                  STOP
             ENDIF
          ENDIF
      ENDIF
      IF( lun2 .NE. 0 ) THEN                                            ! lun2 is used for the kludging together of 2 files
          lun1 = idunit
          CALL podisc(lun2, 2, 900 )                                    ! get by the file headers
      ENDIF
!**** fk stuff needs some special attention!
      CALL mrgfk( TRNSPARM, buf, lbuf, ibuf, scr, lscr, iscr, istop,
     *            ihunit, iptype )
!**** figure whether to read the file as is, or by shots/rps, or by gmt
!**** iptype = 1 means read all data, the user is not controlling
!****        = 2 means read user specified shots/rps
!****        = 3 means read user specified gmt data
      IF( fno+lno+fgmt+fday .EQ. 0 .AND. lday .EQ. 366 .AND.
     *    lgmt .EQ. 2500 ) iptype = 1
      IF( fno+lno .NE. 0 ) iptype = 2
      IF( fday+fgmt.NE.0.OR.lday .NE. 366.OR.lgmt.NE.2500) iptype = 3
      IF( IAND(lprint,2) .NE. 0 ) THEN
          PRINT *, oldpat
          PRINT *, path
          PRINT *, fno, lno, noinc, ftr, ltr, fday, lday, fgmt, lgmt
          PRINT *, gmtinc, fsec, lsec, secinc, renum, secs, jfmt
          PRINT *, sinew, decimf, delayj, intrcs, ntrgat, iptype, set
          PRINT *, rewind, jform, forgat, mintrs, retrac, jsort, ascii
          PRINT *, allno, alltr, notype, noindex, itrtype, itrindex
          PRINT *, nsamps, jntrcs, random
      ENDIF
      IF( path .NE. oldpat ) THEN
          IF( mlists .NE. 1 ) CALL frefil( -2, idunit, istat )
          idarn = 0
!****  another VMS kludge - I don't understand why, but VMS won't let me
!***  use the same unit number on a new file, even though it was closed
          mode = -4
          IF( icompt .EQ. vaxvms .AND. oldpat .NE. ' ') mode = 4
          CALL getfil( mode, idunit, path, istat )
          IF( istat .NE. 0 ) THEN
              PRINT *,' ***  ERROR  ***  Can not open file ',path
              STOP
          ENDIF
          IF( random .NE. 0 ) CALL filsiz( path, filesize )
          oldpat = path
          CALL podisc( ihunit, 1, 0 )                                   ! rewind the header file
!****     If it's an SSC disc file, do the SSC header here
          IF( jform .EQ. 2 ) THEN
              CALL rddisc( idunit, buf, 128, istat )
              IF( icompt .EQ. vaxunix .OR. icompt .EQ. vaxvms ) 
     &            CALL swap16(ibuf,128)
              nsampss = lbuf(1)
              IF( nsampss .EQ. 0 ) THEN
                  nsampss = 2000
                  PRINT *,' *** WARNING  ***  DISKIN is using 2000',
     *            ' samples.'
              ENDIF
              isi = buf(2) * 1000.                                      ! millisecond sample interval
              si = buf(2) / 1000.                                            ! mils to secs
              jntrcs = lbuf(3)
              ntrtot = lbuf(5)                                          ! number of traces total in the SSC file
              IF( ntrtot .EQ. 0 ) ntrtot = 99999
              ifmt = 5                                                  ! set the input data format common to IEEE
              DO i = 1, 3                                           ! get past 3 more 512 byte blocks
                 CALL rddisc( idunit, buf, 128, istat )
              ENDDO
          ENDIF
          IF( path(1:9) .EQ. '/dev/null' ) THEN
              ndeadtr = jntrcs
              IF( sort .NE. 'shot' ) THEN
                  lbuf(7) = 0
                  ndeadtr = 1
              ENDIF
              IF( ftr .NE. -12345 .AND. ltr .NE. 0 ) ndeadtr=ltr-ftr+1
              ndeadrec = 1
              itemp = 1
              IF( noinc .NE. 99999 ) itemp = noinc
              IF( fno .NE. 0 .AND. lno .NE. 0 ) 
     &            ndeadrec = (lno - fno) / itemp + 1
              ndead_done = 0
              ndead = ndeadtr * ndeadrec
              GOTO 200
          ENDIF
!****     The IRIS PASSCAL and SU formats do not have the EBCDIC or binary
!****     headers, so create them so that the output is true SEGY
          IF( jform .EQ. 3 ) THEN
              CALL rddisc( idunit, buf, 60, istat )                     ! get the first trace header
              CALL podisc( idunit, 1, 0 )                               ! set the IRIS file to the beginning
              IF( icompt .EQ. vaxunix .OR. icompt .EQ. vaxvms )
     *            CALL swap16( ibuf(58), 2 )
              nsampss = lbuf(58)
!              IF( nsampss .LT. 32767 ) ibuf(58) = nsampss
              IF( nsampss .LT. 32767 ) THEN
                  CALL long2ushort(nsampps,ibuf(58))
              ELSE
                  PRINT *,' Too many samples for SEG-Y.'
              ENDIF
              numdat = nsampss
              isi = ibuf(59)
              si = FLOAT(isi) / 1000000.
              jntrcs = 1
              IF( jfmt .EQ. 0 ) THEN
                  IF( icompt .EQ. vaxunix .OR. icompt .EQ. vaxvms )
     *                CALL swap16( ibuf(103), 1 )
                  ifmt = ibuf(103)
                  IF( ibuf(103) .EQ. 0 ) ifmt = 3                       ! 16 bit integer
                  IF( ibuf(103) .EQ. 1 ) ifmt = 2                       ! 32 bit integer
              ENDIF
          ENDIF
!****     Create the EBCDIC and Binary SEGY headers for SSC and IRIS and WAV
          IF( jform .EQ. 2 .OR. jform .EQ. 3 .OR. jform .EQ. 6 .OR.
     &       jform .EQ. 12 .OR. jform .EQ. 13 .OR. jform .EQ. 14 ) THEN
!****        died did some of the WAV stuff
              DO 130 i = 1, 40                                          ! build a SEGY EBCDIC header
                 WRITE( cheadr(i), 120 ) i
  120            FORMAT('c ',I2,' ')
                 CALL ascebc( cheadr(i), 80, cheadr(i) )
  130         CONTINUE
              CALL wrdisc( ihunit, ivms, 800 )                          ! write 3200 character from cheadr to disk
              DO i = 1, 200                                             ! build the binary header
                 ibuf(i) = 0
              ENDDO
!****         THE EDIT PHASE MUST SET THESE VARIABLES!!!!!!!!
              ibuf(7) = jntrcs
              ibuf(9) = isi
              ibuf(10) = isi
              ibuf(11) = nsampss
              ibuf(12) = nsampss
              ibuf(13) = ifmt
!              IF( icompt .EQ. 2 .OR. icompt .EQ. 4 ) THEN
!                  IF( jform.NE.5 .AND. jform.NE.7) CALL swap16(ibuf,200)
!              ELSE
!                  IF(jform.EQ.5 .OR. jform.EQ.7) CALL swap16(ibuf,200)
!              ENDIF
              CALL wrdiscb( ihunit, buf, 400 )
!****         skip reading the headers since there aren't any!
              GOTO 200
          ENDIF
!**** diskin edit already put the header of the first file on disk, so it exists!
!**** DO NOT SAVE THESE HEADERS (process header may have modified them).
          CALL podisc( ihunit, 1, 0 )
!****     imvs is equivalenced to cheader
          CALL rddiscb( idunit, ivms, 3200, istat )                     ! get the EBCDIC header
!          CALL wrdiscb( ihunit, ivms, 3200 )
          CALL rddiscb( idunit, scr, 400, istat )                       ! get the binary header
!          IF( icompt .EQ. 2 .OR. icompt .EQ. 4 .OR. jform .EQ. 5) THEN
!              IF( (icompt .EQ. 2 .OR. icompt .EQ. 4) .AND. jform .EQ. 5)
!     &            GOTO 141
!              IF( jform .NE. 7 ) CALL swap16( iscr, 200 )
!  141         CONTINUE
!          ENDIF
!****     ODEC & SWAPPED are in little endian (little byte first)
          IF( icompt .NE. 2 .AND. icompt .NE. 4 ) THEN
              IF( jform .EQ. 5 .OR. jform .EQ. 7 ) CALL swap16(iscr,200)
          ELSE
              IF( jform .NE. 5 .AND. jform .NE. 7) CALL swap16(iscr,200)
          ENDIF
          jntrcs = iscr(7)
          IF( intrcs .NE. 0 ) THEN                                      ! user overriding intrcs?
              jntrcs = intrcs
              iscr(7) = intrcs
          ENDIF
          ifmt = iscr(13)                                               ! user overriding ifmt?
          IF( jfmt .NE. 0 ) THEN
              ifmt = jfmt
              iscr(13) = jfmt
          ENDIF
          CALL ushort2long( iscr(11), nbytes1 )
          nbytes1 = nbytes1 * 4
          IF( ifmt .EQ. 3 ) nbytes1 = nbytes1 / 2
          isort = iscr(15)                                              ! sort: 0 or 1 = shot, 2 = cdp
          idtype = iscr(31)                                             ! data type
          nfktrc = iscr(32)                                             ! number of wavenumbers in the job
          itxsi = iscr(33)                                              ! tx domain sample interval in micros
          itxdel = iscr(34)                                             ! tx domain time delay in mils
          norigtr = iscr(36)                                            ! the number of original traces in tx (before tx2fk)
          segyrev = REAL(iscr(151)) / 256.
          IF( segyrev .LT. 0 .OR. segyrev .GT. 5 ) THEN
              PRINT *,' ***  WARNING  ***  Bad SEG-Y Revision of ',
     &                segyrev
              PRINT *,' Setting segyrev to 0.'
              segyrev = 0
          ENDIF
!****     watch out for rev 0 files that disko will convert to rev 1
          IF( segyrev .GE. 1. ) THEN
              ifixed = iscr(152)
          ELSE
              ifixed = 0
          ENDIF
          ntextual = iscr(153)    ! Number of Extended textual blocks
          IF( ntextual .NE. 0 .AND. segyrev .LT. 1 ) THEN
              PRINT *,' ***  WARNING  ***  Not SEG-Y Rev. 1'
              ntextual = 0
          ENDIF
          addr_tr1 = 3200 + 400
          IF( ntextual .GT. 0 ) addr_tr1 = ntextual * 3200 + addr_tr1
          CALL mrgfk( MRGINIT, buf, lbuf, ibuf, scr, lscr, iscr, istop,
     *                idunit, iptype )                                  ! check for fk polar 
          CALL podiscb( idunit, 1, addr_tr1 )
          CALL podiscb( ihunit, 1, 3600 )
!****     DO NOT REWRITE THE BINARY HEADER!!
!          CALL wrdisc( ihunit, scr, 100 )
          IF( ntextual .EQ. -1 ) THEN
  190         CALL rddiscb( idunit, cheadr, 3200, istat )
              ntextual = ntextual + 1
              DO i = 1, 40
                 IF( cheadr(i)(1:16) .EQ. '((SEG: EndText))' ) GOTO 191
              ENDDO
              GOTO 190
  191         CONTINUE
          ENDIF
          IF( ntextual .GT. 0 ) addr_tr1 = ntextual * 3200 + addr_tr1
      ELSE
          IF( rewind .EQ. 0 ) GOTO 300
!****  that little contortion was due to the trace header already being
!**** in memory if we read a new fno/lno list but the file name was the same.
      ENDIF
      IF( random .NE. 0 .AND. fno .NE. 0 ) THEN
          IF( segyrev .LT. 1 .OR. ifixed .EQ. 0 ) THEN
!              PRINT *,' ***  WARNING  ***  Can not do random I/O - rev',
!     &            rev,' fixed length flag = ',ifixed
              random = 0
          ELSE
              CALL rddisc( idunit, buf, numhdr, istat )
              CALL podisc( idunit, 2, -numhdr )
              IF( icompt .EQ. 2 .OR. icompt .EQ. 4 )
     &            CALL swp_trhdr( ibuf, lbuf )
              no_first = lbuf(3)
              no_wanted = fno
              IF( no_wanted .LE. 0 ) no_wanted = lbuf(3)
              itr_wanted = ftr
!             User may override the sort (shot vs rp) with the SORT parameter.  1=shot
!             isort comes from the binary header, jsort comes from the user.
!      print *,' isrt=',isrt,' jsort=',jsort,' 6=',lbuf(6),' 7=',lbuf(7)
              IF(isort.EQ. 1 .AND. lbuf(7) .NE. 0.AND.jsort .EQ. 0) THEN
                  jsort = 1
                  sort = 'cdp'
              ENDIF
              IF( lbuf(7) .NE. 0 .AND. jsort .NE. 1 ) THEN
                  no_first = lbuf(6)
                  IF( no_wanted .LE. 0 ) no_wanted = lbuf(6)
              ENDIF
          ENDIF
      ENDIF
      IF( IAND(lprint,2) .NE. 0 ) THEN
          PRINT *,' binary hdr sort=',isort,
     &    ' format= ',ifmt,' data starts at',addr_tr1
          PRINT *,' rev=',rev,' fixed=',ifixed,' ntext=',ntextual
      ENDIF
!
  200 CONTINUE
      ibinhdr(31) = idomain
      IF( ndead .NE. 0 ) THEN
          IF( sort .EQ. 'shot' ) THEN
              IF( ndead_done .EQ. 0 ) THEN
                  IF( fno .NE. 0 ) THEN
                      lbuf(3) = fno
                      fno = fno + 1
                  ELSE
                      lbuf(3) = lbuf(3) + 1
                  ENDIF
              ENDIF
              ndead_done = ndead_done + 1
              IF( ftr .NE. -12345 .AND. ltr .NE. 0 ) THEN
                  lbuf(4) = ftr + ndead_done - 1
                  IF( ndead_done .EQ. ltr-ftr+1 ) ndead_done = 0
              ELSE
                  lbuf(4) = ndead_done
                  IF( ndead_done .EQ. jntrcs ) ndead_done = 0
              ENDIF
          ELSE
!              ncdp = ncdp + 1
!              lbuf(7) = ncdp
              lbuf(7) = lbuf(7) + 1
              IF( ndeadtr .EQ. 1 ) lbuf(51) = -1                        ! set the end of gather flag on the last one
          ENDIF
          ndead = ndead - 1
          ibuf(15) = 2                                                  ! set the dead trace flag to dead trace
          DO i = 1, nsampss
             buf(numhdr+i) = 0.
          ENDDO
          RETURN
      ENDIF
      IF( rewind .NE. 0 ) CALL podiscb( idunit, 1, addr_tr1 )
!      IF( rewind .NE. 0 ) THEN
!***      before c, disk addresses started with 1 rather than 0, so when
!***      sioseis did both, we read the headers rather than jjust junping there.
!          CALL podisc( idunit, 1, 0 )                                   ! position after the SEGY headers
!          n =900 
!          CALL rddisc( idunit, scr, n, istat )
!          rewind = 0
!      ENDIF
      n = numhdr
      IF( jform .EQ. 2 ) n = 256                                        ! SSC trace header is 1024 bytes
      IF( iptype .EQ. IPTFKU ) THEN                                     ! if in FK domain then
          CALL mrgfk( MRGTRC, buf, lbuf, ibuf, scr, lscr, iscr, istop,
     *                idunit, iptype )                                  ! fk domain diskin is done here!
          RETURN
      ENDIF
!****
!****  If it's a "sort" file, position the disk to the right trace first
!****
      IF( lunsort .NE. 0 ) THEN
          itemp = lunsort
          IF( ntraces .EQ. 0 ) THEN                                     ! Is this the first traces?
!****         Make a copy of the sort file because something else might want it too
              CALL getfil( 1, lunsort, cheadr(1), istat )               ! open a scratch file
              CALL podisc( itemp, 1, 0 )                                ! rewind
              CALL rddisc( itemp, ivms, 20, istat )
              CALL wrdisc( lunsort, ivms, 20 )
              CALL rddisc( itemp, ntraces, 1, istat )
              CALL wrdisc( lunsort, ntraces, 1, istat )
              DO i = 1, 4
                 CALL rddisc( itemp, apdata, ntraces, istat )
                 CALL wrdisc( lunsort, apdata, ntraces )
              ENDDO
              CALL frefil( 2, itemp, istat )
              CALL podisc( lunsort, 1, 0 )                              ! rewind
              CALL rddisc( lunsort, ivms, 20, istat )
              CALL rddisc( lunsort, ntraces, 1, istat )
          ENDIF
!         newsort is equivalenced to laddress and iflag
          CALL rddisc( lunsort, newsort, 2, istat )
          mtraces = mtraces + 1
          IF( mtraces .EQ. ntraces ) istop = 1
          CALL podisc( idunit, 1, laddress )
      ENDIF
      IF( IAND(lprint,8) .NE. 0 ) THEN
          CALL adrdisc64( idunit, iaddr64 )
          PRINT *,' address = ',iaddr64,' iflag=',iflag,' n= ',n
      ENDIF
!**** Random access requires the user to give FNO
      IF( random .NE. 0 .AND. fno .NE. 0 ) THEN
!****   watch out for 32 bit arithmetic, so force 64
          ltemp64 = ( (no_wanted - no_first ) * jntrcs +
     &            itr_wanted - 1 )
          l2temp64 = (nbytes1+240)
          l3temp64 = 3600
          iaddr64 = ltemp64 * l2temp64 + l3temp64
          IF( iaddr64 .GT. filesize .OR. iaddr64 .LT. 3600 ) THEN
              PRINT *,' ***  ERROR  ***  SEGY Rev 1 indicates shot',
     &         no_wanted,' is not in the file with the first shot of',
     &         no_first
              PRINT *,' file size =',filesize,' wanted byte ',iaddr64
              PRINT *,' Try again with RANDOM 0.'
              PRINT *,' See DISKIN parameter RANDOM.'
              CALL EXIT
          ENDIF
          CALL podiscb64( idunit, 1, iaddr64 )
          itr_wanted = itr_wanted + 1
          IF( itr_wanted .GT. jntrcs ) THEN
              no_wanted = no_wanted + 1
              itr_wanted = ftr
          ENDIF
      ENDIF
      IF( jform .EQ. 12 ) THEN
!****     WAV file input - create a trace header
          DO i = 1, 60
             lbuf(i) = 0
          ENDDO
          CALL podiscb( idunit, 0, 0 )
          CALL rddiscb( idunit, scr, 44, istat )
!****     WAV is little endian (Intel)
          IF( icompt .NE. 2 .AND. icompt .NE. 4 ) THEN
              CALL swap16( iscr(12), 1 )
              CALL swap16( iscr(17), 2 )
              CALL swap32( lscr(7), 1 )
              CALL swap32( lscr(11), 1 )
          ENDIF
          intrcs = iscr(12)
          ltemp = lscr(7)
          si = 1. /  FLOAT(lscr(7))
          isi = NINT(si * 1000000. )
          lbuf(3) = 1
          lbuf(4) = ndone + 1
          ibuf(15) = 1
          ibuf(59) = isi
          nbperw = 4
          nbytes = lscr(11)
          IF( iscr(17) .EQ. 2 ) THEN
              ifmt = 3
              nbperw = 2
!*****        yucko.  nbytes gets divided by 2 if ifmt 3 later, so fake it out
              nbytes = nbytes * 2
          ENDIF
          numdat = lscr(11) / nbperw
          ibuf(58) = numdat
          IF( IAND(lprint,4) .NE. 0 ) THEN
              PRINT *,' WAV: nchans=',iscr(12),' bits per sample=',
     &         iscr(18),' sr=',lscr(7),' nbytes per trace=',lscr(11)
              PRINT *,' numdat=',numdat,' nbperw=',nbperw
          ENDIF
!          CALL rddiscb( idunit, buf(numhdr+1), nbytes, istat )
          istop = 1
          GOTO 330
      ENDIF
      IF( jform .EQ. 13 .OR. jform .EQ. 14 ) THEN
          IF( ndone .GT. 0 ) THEN
              CALL rddiscb( ipunit, path, 200, istat )
              CALL rddisc( ipunit, scr(26), npwrds-25, istat )
              mlists = mlists + 1
          ENDIF
          si = sinew
          DO i = 1, 60
             lbuf(i) = 0
          ENDDO
          isi = NINT(si * 1000000. )
          lbuf(3) = 1
          lbuf(4) = ndone + 1
          ndone = ndone + 1
          ibuf(15) = 1
          ibuf(59) = isi
          buf(49) = si
          nbperw = 4
          numdat = 0
          IF( mlists .EQ. nlists ) istop = 1
      ENDIF
      IF( jform .EQ. 13 ) THEN   !  ASCII file
          OPEN( UNIT=idunit, FILE=path, STATUS='OLD', 
     &          IOSTAT=istat, FORM='FORMATTED' )
          IF( istat .NE. 0 ) THEN
              PRINT *,' ***  ERROR  ***  Could not open file:  ',path
              istop = -1
              RETURN
          ENDIF
          DO i = 1, 10000
             CAll rline( idunit )
             IF( nchars .LE. 0 ) THEN
                 CALL long2ushort( numdat, ibuf(58) )
                 CLOSE( UNIT=idunit )
                 GOTO 350
             ENDIF
  208        CALL getoke( token, nchars1 )
             IF( nchars1 .GT. 0 ) THEN
                 numdat = numdat + 1
                 CALL dcode( token,nchars1,buf(numhdr+numdat),istat)
                 IF( istat .NE. 2 ) THEN
                     PRINT *,' ***  ERROR  ***  File ',path,
     &               'sample ',numdat,token(1:nchars),
     &               ' is not numeric.'
                     GOTO 350
                 ENDIF
                 GOTO 208
             ENDIF
          ENDDO
      ENDIF
      IF( jform .EQ. 14 ) THEN    ! binary file input
          CALL frefil( -2, idunit, istat )
          CALL getfil( 4, idunit, path, istat )
          IF( istat .NE. 0 ) THEN
              PRINT *,' ***  ERROR  ***  Could not open file:  ',path
              istop = -1
              RETURN
          ENDIF
          CALL rddisc( idunit, buf(numhdr+1), 200000, numdat )
          IF( numdat .LE. 0 ) THEN
              PRINT *,' ***  ERROR  ***  Binary file ',path,
     &           ' failed with status = ',numdat
              istop = -1
              RETURN
          ENDIF
          CALL long2ushort( numdat, ibuf(58) )
          GOTO 350
      ENDIF
  209 CONTINUE
      IF( ireal .NE. 0 ) CALL adrdisc64( idunit, iaddr64 )
      CALL rddisc( idunit, buf, n, istat )                              ! read the trace header
      IF( lunsort .NE. 0 .AND. iflag .NE. -99999 ) lbuf(51) = iflag
!****
!****  This kludge is for VMS.  VMS file lengths must be in multiples
!****  512 bytes.  DISKOX fills the last buffer with zeroes, so on
!****  VMS machines, an SEGY header that is all zeroes is end of file!
!****
      IF( icompt .EQ. vaxvms ) THEN
          DO 210 i = 1, 60
             IF( lbuf(i) .NE. 0 ) GOTO 211
  210     CONTINUE
          istat = -1                                                    ! it's a VMS end of file!
  211     CONTINUE
      ENDIF
      IF( istat .NE. n) THEN                                            ! end of file?
!****     REAL-TIME?
          IF( ireal .NE. 0 ) THEN
!             adrdisc return address in BYTES
!ccc              CALL adrdisc( idunit, iaddr )
              CALL frefil( -2, idunit, istat )
              IF( icount .EQ. 1 .OR. icount .EQ. 6 ) 
     &            PRINT *,' Waiting for a shot.'
              CALL SLEEP(5)
              icount = icount + 1
              CALL getfil( -4, idunit, path, istat )
              CALL filsiz( path, filesize )
              CALL podiscb64( idunit, 1, iaddr64 )
              IF( icount .LT. 12 ) GOTO 209
              PRINT *,' Timeout after 60 seconds.'
!              PRINT *,' Control-c to stop.'
!              icount = -380
!              GOTO 209
              istop = -1
              RETURN
          ENDIF
          IF( mlists .EQ. nlists ) THEN                                 ! any more user lists?
              IF( ndone .EQ. 0 ) THEN
                  IF( random .EQ. 1 ) PRINT *,
     &' ***  SUGGESTION:  Nothing found. Try diskin parameter random 0.'
                  IF( allno .EQ. 0 ) PRINT *,
     &' ***  SUGGESTION:  Nothing found. Try diskin parameter allno yes'
              ENDIF
              istop = -1                                                ! return a signal saying no data in buf and STOP
              ibuf(15) = 2                                              ! set the dead trace flag
!              CALL frefil( 2, idunit, istat )
              RETURN
          ENDIF
          GOTO 100
      ENDIF
!**** If it is an SSC file rather than SEGY do it here!
      IF( jform .EQ. 2 ) THEN                                           ! SSC format?
          IF( nhread .GE. ntrtot ) THEN                                 ! there are ntrtot traces in this SSC file
              istop = -1
              ibuf(15) = 2                                              ! set the dead trace flag
              RETURN
          ENDIF
          nhread = nhread + 1                                           ! number of headers read
          IF( icompt .EQ. vaxunix .OR. icompt .EQ. vaxvms ) THEN
              CALL swap32( lbuf, 1 )
              CALL swap16( ibuf(3), 36 )
              CALL swap32( lbuf(22), 1 )
              CALL swap16( ibuf(118), 1 )
          ENDIF
          lcdp = lbuf(1)
          lespn = lbuf(3)
          icdptr = lbuf(4)
          ishot = ibuf(10)
          istrno = ibuf(12)
          ldist = buf(9)
          idelay = buf(30)
          iday = ibuf(100)
          ihour = ibuf(102)
          min = ibuf(104)
          isec = ibuf(106)
          DO i = 1, 60
             lbuf(i) = 0
          ENDDO
          IF( istrno + icdptr .EQ. 0 ) THEN                             ! are both the shot and rp trace number zero?
              IF( lcdp .NE. 0 ) THEN
                  icdptr = 1
              ELSE
                  istrno = 1
              ENDIF
          ENDIF
          lbuf(3) = ishot
          lbuf(4) = istrno
          lbuf(5) = lespn
          lbuf(6) = lcdp
          lbuf(7) = icdptr
          ibuf(15) = 1
          lbuf(10) = ldist
          ibuf(55) = idelay
          ibuf(58) = nsampss
          numdat = nsampss                                               ! set common (some process use numdat rather than ibuf(58))
          ibuf(59) = isi
          ibuf(80) = iday
          ibuf(81) = ihour
          ibuf(82) = min
          ibuf(83) = isec
          GOTO 295
      ENDIF
      IF( icompt .EQ. 2 .OR. icompt .EQ. 4 .OR. jform .EQ. 5 ) THEN
          IF( (icompt .EQ. 2 .OR. icompt .EQ. 4) .AND. jform .EQ. 5 ) 
     &        GOTO 225
          IF( jform .NE. 7 ) CALL swp_trhdr( ibuf, lbuf)
  225     CONTINUE
      ENDIF
      IF( icompt .NE. 2 .AND. icompt .NE. 4 .AND. jform .EQ. 7 )    ! 7 = ODEC
     &    CALL swp_trhdr( ibuf, lbuf)
      IF( jform .EQ. 7 ) THEN
!  ODEC has the shot in lbuf(2) and the trace in lbuf(1), but increments the
!  trace number rather than the shot number, so swap the meaning!
          lbuf(3) = lbuf(2)
          lbuf(4) = lbuf(1)
!****        ibuf(45) = 3        !  the nav are IEEE floating decimal degrees - I want code 5!!!!
          ibuf(36) = -100
!****     arcsec = 60sec/min * 60min/deg = 3600. sec/deg
          lbuf(19) = NINT(buf(19)*3600.*100.)
          lbuf(20) = NINT(buf(20)*3600.*100.)
          lbuf(21) = NINT(buf(21)*3600.*100.)
          lbuf(22) = NINT(buf(22)*3600.*100.)
          ibuf(45) = 2
!  ODEC writes corrupt buffers occassionally, the length is fine, the contents are not.
          IF( lbuf(1) .NE. 1 .AND. ibuf(15) .NE. 1 ) THEN
              nbytes = lastnum * 2 + 320
              CALL podiscb( idunit, 2, nbytes )
              GOTO 209
          ENDIF
          lastnum = ibuf(58)
      ENDIF
      IF( jform .EQ. 9 .OR. jform .EQ. 11 ) THEN                       ! EdgeTech
!****    move the date because diskin writes over long words 45-47
          lbuf(3) = lbuf(1)
          lbuf(4) = lbuf(4) + 1
          IF( jform .EQ. 11 ) lbuf(4) = ibuf(13)
          lbuf(7) = 0
!****     The distance of fish behind antenna (cm)
          lbuf(10) = ibuf(19)
!****     Xstar ver 5 has 5120 bytes per trace (5120-240 = 4880, or 2440 INT*2 words)
          IF( jform .EQ. 11 ) nsamps = 2440
          IF( ibuf(100) .NE. 0 ) THEN 
              IF( icompt .EQ. vaxunix .OR. icompt .EQ. vaxvms )
     &            CALL swap16( ibuf(94), 7 )
!****         Use the NMEA time only if it's there.  o.w. leave it as the cpu time
              IF( ibuf(94) + ibuf(95) + ibuf(96) .NE. 0 ) THEN
                  ibuf(79) = ibuf(100)                                  ! year
                  ibuf(80) = ibuf(99)                                   ! day of year
                  ibuf(81) = ibuf(94)                                   ! hour
                  ibuf(82) = ibuf(95)                                   ! minute
                  ibuf(83) = ibuf(96)                                   ! second
              ELSEIF( iprint1 .LT. 5 ) THEN
                  PRINT *,
     &            ' NMEA time missing - using topside cpu time.'
              ENDIF
              IF( ibuf(98) .NE. 0 ) THEN
                  ibuf(63) = ibuf(97)                                   ! course
                  ibuf(64) = ibuf(98)                                   ! speed in tenths of knots
              ENDIF
          ENDIF
      ENDIF
      IF( ibuf(15) .EQ. 28 ) THEN
          IF( luntr0 .EQ. 0 ) CALL getfil(1,luntr0,cheadr,istat)
          CALL podiscb( luntr0, 1, 0 )
      ENDIF
      IF( ibuf(15) .EQ. 28 .AND. jform .EQ. 4 ) THEN
!**** The Lamont DSS240 system stores stuff such as GMT in trace 0
!**** and the shot number comes from trace 0 (what's in (3) is the file
!**** number.  Trace 0 also counts in the header word for the number of
!**** traces per shot (if there are traces 0-184, then there are 185 tr)
!         There was a screw up early on where trace 0 length was wrong.
!         Trace 0 was always 1024 bytes, unless there was a mistake.
!         Then, we had to change the trace 0 trace length to be the same 
!         as the data traces because the commercial packages couldn't
!         deal with different trace lengths.  SEG-D adopted the old
!         DMS 2000 trace 0 convention and puts a 1 in word 13 and 1
!         214 in word 14 (I*2 word past the segy header).  It should
!         have put an 11 in word 13.  sigh.
!
          IF( ibuf(58) .NE. 1024 .AND. ibuf(79) .LT. 99 )
     &        ibuf(58) = (2048 - 240) / 4
          ltemp = ibuf(58)                                              ! argument must be long integer!
          CALL rddisc( idunit, buf(61), ltemp, istat )                  ! get the rest of trace 0
          IF( istat .NE. ltemp ) THEN
              IF( mlists .NE. nlists ) GOTO 100
              ibuf(15) = 2
              istop = -1
              RETURN
          ENDIF
          CALL wrdiscb( luntr0, buf(61), 4096 )
          IF( jform .EQ. 0 ) intrcs = intrcs - 1
!****     Get the shot number from "section 11", so find section 11
          index = 0
          jform = 4
          IF( ibuf(58) .EQ. 1024 .OR. ibuf(134) .EQ. 214 ) THEN         ! the SIOSEIS created the trace 0
              j = 120 + 13
              DO i = 1, 15
                 IF( ibuf(j+1) .EQ. 214) index=j
!                 IF( ibuf(j) .EQ. 11 .AND. ibuf(j+1) .EQ. 214) index=j
                 j = j + ibuf(j+1)/2
              ENDDO
              DO i = 1, 77
                 isect11(i) = ibuf(index+i+1)
              ENDDO
!****  Hex(3E) gets written if the Lamont system isn't hooked up
              IF( sect11(1:1) .EQ. '>' ) GOTO 245
              CALL ldgo_tr0( sect11(1:1),ldgoshot,ldgodss(1),ldgodss(2),
     &           ldgodss(3), ldgodss(4), ldgodss(5),ldgodss(6),
     &           ship_lat, ship_long, wdepth,
     &           tail_lat, tail_long, tail_dist, tail_bear )
          ELSE
              ldgoshot = 0
              iunits = 1
              itemp = lbuf(62)
              DO 230 i = 1, 8
                 ldgoshot = ldgoshot + IAND(itemp,15) * iunits
                 itemp = LRSHIFT(itemp,4)
                 iunits = iunits * 10
  230         CONTINUE
              ldgodss(1) = 0                                            ! DSS doesn't have year!
              DO i = 1, 4
                 ldgodss(i+1) = ibuf(134+i)
              ENDDO
!**** The following won't compile on the HP, so comment out on HP
              ldgodss(6) = IAND(RSHIFT(ibuf(658),8),i15) * 100
     &               + IAND(ibuf(658),i15) * 10
     &               + IAND(RSHIFT(ibuf(659),8),i15)
          ENDIF
!****     if the shot number is zero, the data is not of that time and
!****     is at least 2 shots old.  Kill it via the dead trace flag and
!****     renumber the shot to what it should have been.
          kill = 1
          IF( ldgoshot .EQ. 0 ) THEN
              lastshot = lastshot + 1
              PRINT *,' shot number 0 becoming shot ',lastshot,
     &        ' at ', (ldgodss(i),i=1,4)
              ldgoshot = lastshot
              kill = 2
          ENDIF
          lastshot = ldgoshot
  245     CALL rddiscb( idunit, buf, 240, istat )                       ! get the next trace header because we already read the trace
          IF( istat .NE. 240 ) THEN
              ibuf(15) = 2
              istop = -1
              RETURN
          ENDIF
      ENDIF
      IF( jform .EQ. 4 .AND. sect11(1:1) .NE. '>' ) THEN
          lbuf(3) = ldgoshot
          lbuf(16) = wdepth
          ibuf(15) = kill
          ibuf(79) = ldgodss(1)
          ibuf(80) = ldgodss(2)
          ibuf(81) = ldgodss(3)
          ibuf(82) = ldgodss(4)
          ibuf(83) = ldgodss(5)
          ibuf(84) = ldgodss(6)
          ibuf(36) = -10
          lbuf(19) = NINT(ship_lat*60.*60.*10.)
          lbuf(20) = NINT(ship_long*60.*60.*10.)
          ibuf(45) = 2
!          buf(50) = wdepth / 750.
!          buf(54) = wdepth
      ENDIF
!**** If no samples and trace id is 0, we have a bad/zero trace header
      IF( ibuf(58) .EQ. 0 .AND. ibuf(15) .EQ. 0 ) GOTO 209
      IF( ltr .GE. 0 .AND. ftr .GE. 0 .AND. jntrcs .EQ. 0 .AND.
     &    trinc .NE. 99999 ) jntrcs = (ltr - ftr) / trinc + 1
!      numdat = ibuf(58)
!      IF( ibuf(58) .EQ. 32767 ) numdat = lbuf(58)                       ! Iris needs long integer
!      IF( IAND(ibuf(58),32768) .NE. 0 ) 
!      IF( ibuf(58) .LT. 0 )
!     &    numdat = IAND(ibuf(58),i32767) + 32768
!      IF( ibuf(58) .LT. 0 ) numdat = lbuf(29)
      CALL ushort2long( ibuf(58), numdat )
      iday = ibuf(80)
      ihour = ibuf(81)
      imin = ibuf(82)
      isec = ibuf(83)
  295 IF( IAND(lprint,2) .NE. 0 ) 
     *    PRINT *,' shot',lbuf(3),' trace ',lbuf(4),' rp ',lbuf(6),
     &    ' trace ',lbuf(7),' numdat ',numdat,' tr id ',ibuf(15)
      IF( nsamps .GT. 0 ) THEN                                         ! nsamps is a user param
          numdat = nsamps
!          ibuf(58) = nsamps
          CALL long2ushort( nsamps, ibuf(58) )
      ENDIF
      nwrds = numdat
      nbytes = nwrds * 4                                                ! watch out for an odd number of 16 bit words!
!**** The DSS 240 trace data has a trace 0 with an id of 28 and nsamps of 0,
!**** but it really is a 2048 byte record, so we already read 240 of the 2048.
!**** The 1996 mar/sigma data recorded on the DMS2000 saved trace 0 during
!**** segdin as a 4096 bytes trace (1024 samples) excluding the trace header.
      IF( ibuf(15) .EQ. 28 .AND. numdat .NE. 1024 .AND. jform .EQ. 4 )
     &     nwrds = 452
!****  if the data isn't seismic, skip it
!****  the trace id (ibuf(15)) is 0 if people are lazy and didn't set it
!**** 1 = seismic data, 2=dead, 3 = dummy, 6 = sweeep, 9 = UTIG code
!****  28 = LDEO trace 0 stuff, 100 means GPR (ground penetrating radar)
      IF( ibuf(15) .EQ. 100 ) ibuf(15) = 1
      IF( ibuf(15) .NE. 0 .AND. ibuf(15) .NE. 1 .AND.
     *    ibuf(15) .NE. 2 .AND. ibuf(15) .NE. 3 .AND.
     *    ibuf(15) .NE. 6 .AND. ibuf(15) .NE. 9 ) GOTO 310
!****   VMS ftp fills the end of a 512 byte block with zeroes, so the
!****  rddisc might not know when end of file is!!  To boot, when this
!****  file is transfered to another machine, the problem is there too!
!****  So, if numdat, the number of samples in the trace, is zero, make
!****  it the end of the job!  All hell is going to break loose with
!****  zero samples!   - why?
      IF( numdat .LE. 0 ) THEN
          IF( mlists .NE. nlists ) GOTO 100
!          istop = -1
          ibuf(15) = 2                                                  ! set the dead trace flag
          RETURN
      ENDIF
!****
!****  If LNTODO, back off the end of file N traces
!****  Caution:  When doing realtime it is possible that the inode is updated and
!****  the trace read before the trace is fully written.  This is especially possible
!****  if LNTODO = 1 or the data is being written across a slow network.  It happened
!****  with segd - the later part of the trace or shot was garbage - and of course
!****  not repeatable.
!****
      IF( lntodo .NE. 0 ) THEN
          ltemp64 = 3600
          iaddr64 = filesize - ltemp64
          ltemp64 = (nbytes+240)
          l2temp64 = lntodo
          iaddr64 = iaddr64 / ltemp64 * ltemp64
          iaddr64 = iaddr65 - l2temp * ltemp64
          ltemp64 = 3600
          iaddr64 = iaddr64 + ltemp64
          IF( iaddr64 .LT. 3600 ) iaddr64 = 3600
          lntodo = 0
          CALL podiscb64( idunit, 1, iaddr64 )
          GOTO 209
      ENDIF
!****
!****   The right trace has been read, now set the sort type (shot vs rp vs stack)
!****
      IF( lbuf(7) .EQ. 0 ) THEN                                         ! is it shot or rp sorted?
          sort = 'shot'
!      ELSEIF( ibuf(17) .EQ. 0 ) THEN                                    ! cdp fold
!          sort = 'cdp'
      ELSEIF( lbuf(3)+lbuf(4)+lbuf(10) .EQ. 0 ) THEN
          sort = 'stack'
      ELSE
          sort = 'cdp'
      ENDIF
!**** well, that works most of the time, but Palmer sets ibuf(17) to 1 on
!**** shot data.  diskox doesn't set the binary header because it doesn't
!**** know whether it's before or after gather or stack.  gather and stack
!**** can't set it because there may be a diskox before it .  Sigh.  So,
!**** stack normally (parameter HEADER, which preset=NORMAL) sets the shot
!**** number, shot trace number, and range to zero.
      IF( JSORT .ne. 0 ) THEN                                           ! is the user overiding?
          IF( jsort .EQ. 1 ) sort = 'shot'
          IF( jsort .EQ. 2 ) sort = 'cdp'
          IF( jsort .EQ. 3 ) sort = 'stack'
      ENDIF
      IF( sort .EQ. 'shot' ) THEN
          no = lbuf(3)
          itrno = lbuf(4)
      ELSE
          no = lbuf(6)
          itrno = lbuf(7)
      ENDIF
      IF( notype .GT. 0 ) THEN
          IF( notype .EQ. 1 ) no = ibuf(noindex)
          IF( notype .EQ. 2 ) no = lbuf(noindex)
          IF( notype .EQ. 3 ) no = buf(noindex)
      ENDIF
      IF( itrtype .GT. 0 ) THEN
          IF( itrtype .EQ. 1 ) itrno = ibuf(itrindex)
          IF( itrtype .EQ. 2 ) itrno = lbuf(itrindex)
          IF( itrtype .EQ. 3 ) itrno = buf(itrindex)
      ENDIF
!**** set the shot trace number to 1 when it's 0 along with the rp num and rp tr no.
      IF( lbuf(3) .NE. 0 .AND. lbuf(4)+lbuf(6)+lbuf(7).EQ.0) THEN
          lbuf(4) = 1
          itrno = 1    ! noinc needs this
      ENDIF
      IF( forgat .NE. 0. .OR. ntrgat .GE. 0 ) lbuf(51) = 0
      IF( IAND(lprint,2) .NE. 0 ) THEN
          PRINT *,' no=',no,' itrno=',itrno,' fno=',fno,' lno=',lno,
     *            ' ftr=',ftr,' ltr=',ltr,' nextno=',nextno
          PRINT *,' nexttr=',nexttr,' iptype=',iptype,' idtype=',
     *            idtype,' jsort=',jsort
      ENDIF
!****
!****  find the right trace ( a little patchwork since I had GOTOs into
!****       the middle of an IF-THEN-ENDIF block (310), which the CRAY
!****       objects to)
!****
  300 CONTINUE
      IF( iptype .EQ. 1 ) THEN                                          ! reading any and all data?
          IF( ftr .EQ. -12345 .AND. ltr .EQ. 0 ) THEN
              IF( trinc .EQ. 99999 )  GOTO 330                          ! ftr, ltr and trinc as presets
          ENDIF
          IF( itrno .LT. ftr ) GOTO 310     ! skip
          IF( itrno .LT. nexttr .AND. trinc .NE. 99999 ) GOTO 310     ! skip
          IF( itrno .GT. ltr ) THEN
              nexttr = ftr
              GOTO 310     ! skip
          ENDIF
          IF( itrno .NE. nexttr .AND. trinc .NE. 99999 ) THEN
              nexttr = nexttr + trinc
              IF( itrno .GE. ltr ) nexttr = ftr
              GOTO 310
          ENDIF
          IF( trinc .NE. 99999 ) nexttr = nexttr + trinc
          IF( itrno .GE. ltr ) nexttr = ftr
          GOTO 330
      ENDIF
      IF( iptype .EQ. 2 ) THEN                                          ! reading by shots/rps?
          IF( fno+lno .EQ. 0 .AND. noinc .EQ. 99999 ) GOTO 330     ! keep this trace - read everything
          IF( no .GT. lno .AND. lno .NE. 0 ) THEN
              IF( random .EQ. 0 ) GOTO 310   !  keep looking
              IF( nlists .EQ. mlists ) THEN
                  IF( allno .NE. 0 ) GOTO 330    ! keep this trace because we haven't reached lno
                  istop = -1
                  ibuf(15) = 2                                          ! set the dead trace flag
                  CALL frefil( 2, idunit, istat )                       ! close the file, release the unit
                  RETURN
              ENDIF
              GOTO 100    ! get another parameter list
          ENDIF
!****     nextno and nexttr don't work if it's a sort input
          IF( no .LT. nextno .AND. lunsort .EQ. 0 ) GOTO 310   ! skip the trace
          IF( itrno .LT. ftr ) GOTO 310
          IF( ltr .NE. 0 ) THEN
              IF( itrno .GT. ltr ) GOTO 310
              IF( itrno .EQ. ltr .AND. forgat .EQ. 0 ) THEN
                  IF( noinc .NE. 99999 ) nextno = nextno + noinc
                  nexttr = ftr
                  GOTO 305
              ENDIF
          ELSE
              IF( sort .EQ. 'shot' ) THEN
                  IF( itrno .EQ. jntrcs ) THEN
                      IF( noinc .NE. 99999 ) nextno = nextno + noinc
                      nexttr = ftr
                      GOTO 305
                  ELSE
                      IF( noinc .NE. 99999 .AND. jntrcs .EQ. 0 ) THEN
          PRINT *,' ***  ERROR  ***  SEG-Y header value for ntrcs is 0.'
                 PRINT *,' NOINC requires SEG-Y header ntrcs to be set.'
                  PRINT *,' Use DISKIN parameter NTRCS to set it.'
                          CALL EXIT
                       ENDIF
                  ENDIF
              ELSE
                  IF( lbuf(51) .EQ. -1 .OR. sort .EQ. 'stack' ) THEN
                      IF( noinc .NE. 99999 ) nextno = nextno + noinc
                      nexttr = ftr
                      GOTO 305
                  ENDIF
              ENDIF
          ENDIF
  301     IF( trinc .EQ. 99999 ) GOTO 330
          IF( trinc .NE. 1 ) THEN
              IF( itrno .NE. nexttr .AND. lunsort .EQ. 0 ) GOTO 310
              nexttr = nexttr + trinc
              IF( ltr .GT. 0 ) THEN
                  IF( nexttr .GT. ltr ) THEN
                      IF( noinc .NE. 99999 ) nextno = nextno + noinc
                      nexttr = ftr
                  ENDIF
              ELSE
                  IF( nexttr .GT. jntrcs .OR. lbuf(51) .EQ. -1 ) THEN
                      IF( noinc .NE. 99999 ) nextno = nextno + noinc
                      nexttr = ftr
                   ENDIF
              ENDIF
          ENDIF
  305 CONTINUE
      ENDIF
      IF( iptype .EQ. 3 ) THEN                                          ! reading by gmt?
          IF( fday .EQ. 0 ) fday = ibuf(80)
          igmt = ihour*100 + imin
          IF( iday .LT. fday ) GOTO 310
          IF( igmt .LT. fgmt ) THEN
              IF( fday .NE. 0 .AND. iday .LE. fday ) GOTO 310
          ENDIF
          IF( iday .EQ. fday .AND. igmt .EQ. fgmt .AND. isec .LT. fsec)
     *        GOTO 310
          fgmt = 0                                                      ! successive days should start at the beginning
          IF( secinc .NE. 0 ) THEN
              IF( isec .NE. fsec ) GOTO 310
              fsec = fsec + secinc
              IF( fsec .GE. 60 ) fsec = fsec - 60
          ENDIF
          IF( iday .GT. lday ) THEN
              IF( mlists .EQ. nlists ) THEN
                  istop = -1
                  RETURN
              ENDIF
              GOTO 100
          ENDIF
          IF( iday .GE. lday .OR. lday .EQ. 366) THEN 
              IF( igmt .GT. lgmt .OR. iday .GT. lday .OR.
     *               (igmt .EQ. lgmt .AND. isec .GT. lsec) ) THEN
                  IF( mlists .EQ. nlists ) THEN
                      istop = -1
                      ibuf(15) = 2                                          ! set the dead trace flag
                      CALL frefil( 2, idunit, istat )                       ! close the file, release the unit
                      RETURN
                  ENDIF
                  GOTO 100
              ENDIF
          ENDIF
          IF( ftr .EQ. -12345 .AND. ltr .EQ. 0 .AND. trinc .EQ. 1 ) 
     &        GOTO 330                                                  ! ftr, ltr and trinc as presets
          IF( itrno .LT. nexttr ) GOTO 310
          IF( itrno .GT. ltr .AND. ltr .NE. 0 ) THEN
              nexttr = ftr
              GOTO 310
          ENDIF
          IF( trinc .EQ. 99999 ) GOTO 330
          IF( itrno .NE. nexttr ) THEN
              nexttr = nexttr + trinc
              IF( itrno .GE. ltr ) nexttr = ftr
              GOTO 310
          ENDIF
          nexttr = nexttr + trinc
          GOTO 330
      ENDIF
      GOTO 330                                                          ! it must be a trace to keep!
!****
!****    skip the data portion of the trace
!****
  310 CONTINUE
      IF( ifmt .EQ. 3 .OR. ifmt .EQ. 4 ) THEN                           ! 16 bit integer?
          nwrds = nwrds / 2
          nbytes = nbytes / 2
      ENDIF
      IF( jform .EQ. 7 ) nbytes = nbytes + 320
      IF( jform .EQ. 2 ) nwrds = ((nsampss-1) / 128 + 1) * 128          ! SSC data is in 512 bytes blocks
  320 IF( icompt .NE. 4 ) THEN                                          ! VMS VAX doesn't have byte oriented I/O
          CALL rddiscb( idunit, scr, nbytes, istat )
      ELSE
          CALL rddisc( idunit, scr, nwrds, istat )
          istat = istat * 2                                             ! kludge!!
      ENDIF
      IF( istat .NE. nbytes ) THEN
          IF( ireal .NE. 0 ) THEN
              CALL adrdisc64( idunit, iaddr64 )
              CALL frefil( -2, idunit, istat )
              CALL SLEEP(1)
              icount = icount + 1  
              IF( icount .EQ. 1 ) PRINT *,' Waiting for a shot.'
              CALL getfil( -4, idunit, path, istat )
              CALL podiscb64( idunit, 1, iaddr64 )
              IF( icount .LT. 100 ) GOTO 320
          ENDIF
          istop = -1
          ibuf(15) = 2                                                  ! set the dead trace flag
          CALL frefil( 2, idunit, istat )                               ! close the file, release the unit
          RETURN
      ENDIF
      GOTO 200
!****
!****   We got a trace to keep
!****
  330 CONTINUE
      ndone = ndone + 1
      IF( ntodo .EQ. ndone ) istop = 1
      IF( ibuf(59) .EQ. 0 ) ibuf(59) = ibinhdr(9)
      si = REAL(ibuf(59))/1000000.
      IF( ibuf(59) .EQ. 7812 ) si = 7.8125E-03                          ! 1/128 is really .0078125, not .007812
      IF( ibuf(59) .EQ. 488 ) si = 1. / 2048.
      IF( ibuf(59) .EQ. 244 ) si = 1. / 4096.
      IF( ibuf(59) .EQ. 122 ) si = 1. / 8192.
      IF( ibuf(59) .EQ. 83 ) si = 1. / 12000.
      IF( ibuf(59) .EQ. 67 ) si = 1. / 14836.79525                      ! is really 67.4
      IF( ibuf(59) .EQ. 63 ) si = 1. / 16000.
      IF( ibuf(59) .EQ. 54 ) si = 1. / 18518.51852
      IF( ibuf(59) .EQ. 45 ) si = 1. / 22222.22222
      IF( ibuf(59) .EQ. 41 ) si = 1. / 24000.
      IF( ibuf(59) .EQ. 39 ) si = 1. / 25510.2048
      IF( ibuf(59) .EQ. 31 ) si = 3.125E-05                             ! 1./32000.
!****   Knudsen sample intervals are:
!   3ms chirp:  fs = 25510.2048:    interval = 39
!   6ms chirp:  fs = 22222.22222:    interval = 45
!   12ms chirp:  fs = 18518.51852:    interval = 54
!   24ms chirp:  fs = 14836.79525:    interval = 67
!*****
!*****  Try and detect and correct a bug corrected in 2005 in segdex/segddex
!*****
      IF( ibuf(36) .EQ. -10 ) THEN
          IF( ABS(lbuf(19))/36000. .GT. 180. .OR.
     &        ABS(lbuf(20))/36000. .GT. 90. ) THEN
              lbuf(19) = NINT(DFLOAT(lbuf(19)) / .6 )
              lbuf(20) = NINT(DFLOAT(lbuf(20)) / .6 )
              lbuf(21) = NINT(FLOAT(lbuf(21)) / .6 )
              lbuf(22) = NINT(FLOAT(lbuf(22)) / .6 )
              ibuf(36) = -100
          ENDIF
      ENDIF
!
      buf(49) = si
      IF( sinew .GT. 0 ) THEN
          si = sinew                                                    ! is the user overriding the header?
          buf(49) = si
          temp = si * 1000000. + .5                                     ! convert the new to mics
          ibuf(59) = temp
      ENDIF
!****
!****  Do NTRGAT, then RENUM, then RETRAC
!****  If NTRGAT is given then the rp and trace numbers will be different
!****  IF RENUM is given, then
      IF( ntrgat .GT. 0 ) THEN
          IF( lbuf(6) .EQ. 0 ) THEN
!****      when ntrgat, use the rp number if already there.
!****      if not there use the shot number if it exists.
!****      if no shot number either, start with 1
              IF( lbuf(4) .EQ. 0 ) THEN
                  lbuf(6) = numcdp
              ELSE
                  lbuf(6) = lbuf(3)
              ENDIF
          ENDIF
          mtrgat = mtrgat + 1
          lbuf(7) = mtrgat
          sort = 'cdp'
          IF( mtrgat .EQ. ntrgat ) THEN
              lbuf(51) = -1
              mtrgat = 0
              numcdp = numcdp + 1
          ENDIF
      ENDIF
      IF( renum .NE. 0 ) THEN
          IF( sort .EQ. 'shot' ) THEN
              lbuf(3) = renum
              lbuf(7) = 0
          ELSE
              lbuf(6) = renum
          ENDIF
      ENDIF
      IF( retrac .NE. 0 ) THEN
          IF( sort .EQ. 'shot' ) THEN
              lbuf(4) = retrac
              lbuf(7) = 0
          ELSE
              lbuf(7) = retrac
          ENDIF
          retrac = retrac + 1 
      ENDIF
      IF( renum + retrac .NE. 0 ) THEN                                  ! is either given?
          IF( sort .EQ. 'shot' ) THEN
!****    retrac other than 1 will cause the following to cause havic
              IF( lbuf(4)+1 .GT. jntrcs .OR.
     &            ( ltr .NE. 0 .AND. lbuf(4) .EQ. ltr )) THEN
                  IF( lretrac .NE. 0 ) retrac = lretrac
                  IF( renum .NE. 0 ) renum = renum + 1
              ENDIF
          ELSE
              IF( lbuf(51) .EQ. -1 .OR. sort .EQ. 'stack' ) THEN
                  IF( lretrac .NE. 0 ) retrac = lretrac
                  IF( renum .NE. 0 ) renum = renum + 1
              ENDIF
          ENDIF
      ENDIF
      IF( jform .EQ. NIUST ) THEN
          CALL swap32( lbuf(46),1)
          CALL swap16( ibuf(91),1)
          ltemp = ibuf(53) * 1000. + ibuf(91)    ! ltemp is number of microseconds
          IF( ltemp .GT. 32767 ) THEN
              ibuf(55) = ltemp / 100 
              ibuf(108) = -10
          ELSEIF( ltemp .GT. 3276 ) THEN
              ibuf(55) = ltemp / 10
              ibuf(108) = -100
          ELSE
              ibuf(55) = ltemp
              ibuf(108) = 0
          ENDIF
          CALL swap16( ibuf(91),1)
          CALL swap32( lbuf(46),1)
!          lbuf(4) = 1            ! just in case the user doesn't do it
          segyrev = 1.
!          ibinhdr(151) = 256     ! doex automatically sets it to 256 also so don't do it
          ifixed = 1
!          ibinhdr(152) = ifixed  ! doex does this also so don't do it
      ENDIF
      delay = REAL(ibuf(55)) / 1000.
      IF( segyrev .GE. 1 .AND. ibuf(108) .NE. 0 ) THEN
          scalar_mils = ibuf(108)
          IF( scalar_mils .GT. 0 ) THEN
              delay = delay * scalar_mils
          ELSE
              delay = delay / (-1 * scalar_mils)
          ENDIF
      ELSE
!****     disko will try to write rev 1 and if 108 is non-zero then the delay scalar is set!
!****        (Knudsen rev 0 uses 108.  Subsequent sioseis processing will make the output rev 1)
          ibuf(108) = 0
      ENDIF
!****  ah sh!@#$%,  Geometrics allows negative delays and Knudsen
!****  writes 16 bit unsigned delays and UTIG OBS use 32 bit signed.
!**** If the delay is negative, then it was an unsigned 16 bit integer
!**** UTIG has this clever scheme mm making the delay into a 32 bit integer
!**** by using bytes 107,108,109,110 for the delay.
!      IF( ibuf(54) .NE. 0 .AND. ibuf(55) .NE. 0 ) THEN
      IF( jform .EQ. 10 ) THEN                                          ! UTIG-OBS
          IF( icompt .EQ. 2 .OR. icompt .EQ. 4 ) CALL swap16(ibuf(54),2)
          ishort(1) = ibuf(54)
          ishort(2) = ibuf(55)
          IF( icompt .EQ. 2 .OR. icompt .EQ. 4 ) CALL swap32(llong,1)
!          delay = FLOAT(llong) / 1000.
      ENDIF
      buf(46) = delay
      IF( ifmt .EQ. 3 .OR. ifmt .EQ. 4 ) THEN                           ! 16 bit integer?
          nwrds = nwrds / 2
          nbytes = nbytes / 2
      ENDIF
      IF( secs .GT. 0 ) THEN
          numdat_orig = numdat
          numdat = secs/si + 1                                          ! is the user changing the amount of data to process?
!          ibuf(58) = numdat
           CALL long2ushort( numdat, ibuf(58) )
      ENDIF
      IF( jform .EQ. 2 ) nwrds = ((nsampss-1) / 128 + 1) * 128          ! SSC data is in 512 bytes blocks
      IF( jform .EQ. 7 ) nbytes = nbytes + 320
      IF( IAND(lprint,4) .NE. 0 ) THEN
          PRINT *,' no=',no,' itrno=',itrno,' numdat=',numdat,' ntrcs=',
     *            jntrcs,' ifmt=',ifmt,' nbytes=',nbytes
      ENDIF
  340 IF( numdat .GT. maxsamps ) THEN
          PRINT *,' ***  WARNING  ***  Too much data!  Want ',numdat,
     &        ' maximum is ',maxsamps
      ENDIF
      IF( ireal .NE. 0 ) CALL adrdisc64( idunit, iaddr64 )                  ! get the disk address before reading attempt
      IF( icompt .NE. 4 ) THEN                                          ! VMS VAX doesn't have byte oriented I/O
          CALL rddiscb( idunit, lbuf(numhdr+1), nbytes, istat )
      ELSE
          CALL rddisc( idunit, lbuf(numhdr+1), nwrds, istat )
          istat = istat * nbperw                                        ! kludge!!
      ENDIF
      IF( istat .NE. nbytes ) THEN
          IF( ireal .NE. 0 ) THEN
!              CALL adrdisc( idunit, iaddr )
              CALL frefil( -2, idunit, istat )
              CALL SLEEP(1)
              icount = icount + 1  
              IF( icount .EQ. 1 ) PRINT *,' Waiting for a shot.'
              CALL getfil( -4, idunit, path, istat )
              CALL podiscb64( idunit, 1, iaddr64 )
              IF( icount .LT. 100 ) GOTO 340
          ENDIF
          PRINT *,' ***  ERROR  ***  Error in file ',path
          PRINT *,' while reading shot/rp ',no,' trace ',itrno,
     *            ' nbytes = ',nbytes,' istat=',istat
          IF( mlists .NE. nlists ) PRINT *,
     &        ' Is DISKIN LNO missing from a list?'
          istop = -1
          ibuf(15) = 2                                                  ! set the dead trace flag
          CALL frefil( 2, idunit, istat )                               ! close the file, release the unit
          RETURN
      ENDIF
!****
!****  decimate the data if the user requested it
!****
  350 CONTINUE
      IF( decimf .GT. 1 ) THEN
          IF( ibuf(15) .NE. 2 ) THEN                                    ! is it a dead trace?
              IF( ifmt .NE. 3 ) THEN                                    ! is it 16 bit data?
                  i = 0
                  DO 400 j = 1, numdat, decimf
                     i = i + 1
                     lbuf(numhdr+i) = lbuf(numhdr+j)
  400             CONTINUE
              ELSE
                  i = 0
                  DO 410 j = 1, numdat, decimf
                     i = i + 1
                     ibuf(numhdr+numhdr+i) = ibuf(numhdr+numhdr+j)
  410             CONTINUE
              ENDIF
          ENDIF
          numdat = (numdat+1) / decimf
          buf(49) = buf(49) * decimf
          si = si * decimf
          ibuf(59) = ibuf(59) * decimf
          CALL long2ushort( numdat, ibuf(58) )
      ENDIF
!****
!****  Now reformat the data
!****
      IF( ifmt .EQ. 1 ) THEN                                            ! IBM floating point
          IF( icompt .EQ. vaxunix .OR. icompt .EQ. vaxvms ) THEN
              IF( jform .NE. 5 .AND. jform .NE. 7 )
     *            CALL swap32( buf(numhdr+1), numdat )
          ELSE
              IF( jform .EQ. 5 .OR. jform .EQ. 7 )
     *            CALL swap32( buf(numhdr+1), numdat )
          ENDIF
          IF( iuseap .EQ. 0 ) THEN
              in = 0
              CALL ibm2fp( buf(numhdr+1), numdat, buf(numhdr+1) )
          ELSE
              in = 1
              CALL inap( buf, numdat )                                  ! assign the ap
              CALL apput( buf(numhdr+1), in, numdat, 3 )                ! convert on the fly
          ENDIF
          GOTO 1000
      ENDIF
      IF( ifmt .EQ. 2 ) THEN                                            ! 32 bit integers?
          in = 0
          IF( icompt .EQ. vaxunix .OR. icompt .EQ. vaxvms ) THEN
              IF( jform .NE. 5 .AND. jform .NE. 7 .AND. jform .NE. 12 )
     *            CALL swap32( buf(numhdr+1), numdat )
          ELSE
              IF( jform .EQ. 5 .OR. jform .EQ. 7 .OR. jform .EQ. 12 )
     *            CALL swap32( buf(numhdr+1), numdat )
          ENDIF
          DO  i = 1, numdat
              buf(numhdr+i) = lbuf(numhdr+i)
          ENDDO
          GOTO 1000
      ENDIF
      IF( ifmt .EQ. 3 ) THEN                                            ! 16 bit integers?
          in = 0
          IF( icompt .EQ. vaxunix .OR. icompt .EQ. vaxvms ) THEN
              IF( jform .NE. 5 .AND. jform .NE. 7 .AND. jform .NE. 12 )
     *            CALL swap16( ibuf(numhdr+numhdr+1), numdat )
          ELSE
              IF( jform .EQ. 5 .OR. jform .EQ. 7 .OR. jform .EQ. 12 )
     *            CALL swap16( ibuf(numhdr+numhdr+1), numdat )
          ENDIF
          IF( jform .NE. 7 ) THEN
              DO i = 1, numdat
                 scr(i) = ibuf(numhdr+numhdr+i)
              ENDDO
          ELSE
!****         The ODEC data is unsigned 16 bits
              DO i = 1, numdat-80                                       ! don't do the trailer
                 ltemp = ibuf(numhdr+numhdr+i)
                 ltemp = IAND(ltemp,32767)
                 IF( ibuf(numhdr+numhdr+i) .LT. 0 ) ltemp = ltemp+32768
                 scr(i) = ltemp
              ENDDO
              DO i = 1, 80                                              ! now move the trailer
                 lscr(numhdr+numdat-80+i) = lbuf(numhdr+numdat/2-80+i)
              ENDDO
          ENDIF
!****     the next stupid statement gets around an HP compiler bug!
          temp = scr(numdat)
          DO i = 1, numdat
             buf(numhdr+i) = scr(i)
          ENDDO
          GOTO 1000
      ENDIF
      IF( ifmt .EQ. 4 ) THEN                                            ! 16 bit integers?
          in = 0
          IF( icompt .EQ. vaxunix .OR. icompt .EQ. vaxvms ) THEN
              IF( jform .NE. 5 .AND. jform .NE. 7 )
     *            CALL swap16( buf(numhdr+1), numdat )
          ELSE
              IF( jform .EQ. 5 .OR. jform .EQ. 7 )
     *            CALL swap16( buf(numhdr+1), numdat )
          ENDIF
          CALL sfp2fp( ibuf( numhdr+numhdr+1 ), numdat, scr )
          DO i = 1, numdat
             buf(numhdr+i) = scr(i)
          ENDDO
          GOTO 1000
      ENDIF
!****
!****   IEEE - Native word floating point
      IF( ifmt .GT. 4 ) THEN
          in = 0
          IF( icompt .EQ. vaxunix .OR. icompt .EQ. vaxvms ) THEN
              IF( jform .NE. 5 .AND. jform .NE. 7 )
     *            CALL swap32( buf(numhdr+1), numdat )
          ELSE
              IF( jform .EQ. 5 .OR. jform .EQ. 7 )
     *            CALL swap32( buf(numhdr+1), numdat )
          ENDIF
      ENDIF
!****
!****   Now discard the portions of the data that the user doesn't want
!****
 1000 CONTINUE
      istart = 0                                                        ! this is an additive to an index to the start of data
      IF( set(2) .NE. 0 ) THEN
          istart = (set(1)-delay) / si + .5
          IF( delay .GT. set(1) ) THEN
              n = NINT( (delay-set(1)) / si )
              IF( n+numdat .GT. 33000 ) THEN
                  PRINT *,' Bad delay or set. delay=',delay,' set=',set
                  n = 5000
                  ibuf(15) = 2
                  GOTO 1031
              ENDIF
              DO i = 1, n
                 scr(i) = 0.
              ENDDO
              DO i = 1, numdat
                 scr(n+i) = buf(numhdr+i)
              ENDDO
              numdat = numdat + n
              DO i = 1, numdat
                 buf(numhdr+i) = scr(i)
              ENDDO
 1031         istart = 0
          ENDIF
          buf(46) = set(1)
          delay = set(1)
          ibuf(55) = set(1) * 1000. + .5                            ! reset the delay in the header
          IF( set(2) .GT. delay + (numdat-1)*si ) THEN
              n = NINT( (set(2)-delay)/si  ) + 1 - numdat
              DO i = 1, n
                 buf(numhdr+numdat+i) = 0.
              ENDDO
          ENDIF
          numdat = NINT((set(2)-set(1))/si) + 1
!          ibuf(58) = numdat
          CALL long2ushort( numdat, ibuf(58) )
          IF( istart .GT. 0 ) THEN
              DO i = 1, numdat
                 buf(numhdr+i) = buf(numhdr+istart+i)
              ENDDO
          ENDIF
      ENDIF
!****
!****  Now zero the back of the trace if secs made it long
!****
      IF( secs .GT. numdat_orig * si ) THEN
          DO i = numdat_orig + 1, numdat
             buf(numhdr+i) = 0.
          ENDDO
      ENDIF
!****
!****  Take care of non SIOSEIS gather stuff here
!****
      IF( forgat .NE. 0. ) THEN
          IF( lbuf(51) .EQ. -1 ) THEN
              nextno = nextno + noinc
              nexttr = ftr
          ENDIF
          CALL rddisc( idunit, scr, 7, istat )                          ! read the next trace header
          IF( (icompt .EQ. vaxunix) .OR. (icompt .EQ. vaxvms) .OR.
     &        (jform .EQ. 5) .OR. (jform .EQ. 7) )
     &        CALL swap32( lscr, 7 )
          lbuf(51) = 0
          IF( istat .LT. 0 ) lbuf(51) = -1                              ! end of file?
          IF( ltr .NE. 0 .AND. itrno .EQ. ltr ) lbuf(51) = -1           ! ltr given and this is ltr
          IF( jform .EQ. 2 ) THEN                                       ! SSC gathers?
              IF( lbuf(6) .NE. lscr(1) ) iforcnt = iforcnt + 1
          ELSE
              IF( lbuf(6) .NE. lscr(6) ) iforcnt = iforcnt + 1
          ENDIF
          IF( iforcnt .EQ. forgat ) THEN
              lbuf(51) = -1
              iforcnt = 0
          ENDIF
          CALL podisc( idunit, 2, -7 )                                  ! position the disk back to the start
      ENDIF
!****
!****  If using two input files
!****
      IF( lun2 .NE. 0 ) THEN                                            ! one or two files simultaneously?
          IF( ialt .EQ. 2 ) THEN
              IF( idunit .EQ. lun1 ) THEN
                  idunit = lun2
              ELSE
                  idunit = lun1
              ENDIF
          ENDIF
          IF( itrno .EQ. ltr ) THEN
              IF( idunit .EQ. lun2 ) THEN                               ! are we on the second file?
                  idunit = lun1                                         ! switch back to the first file
              ELSE
                  nextno = no                                           ! look for the same shot number
                  nexttr = ftr
                  lbuf(51) = 0
                  ndeadtr = 0
                  idunit = lun2                                         ! read from the second file
              ENDIF
          ENDIF
      ENDIF
!      IF( lbuf(7) .NE. 0 ) ncdp = ncdp + 1                              ! count the traces in this rp
      IF( mintrs .NE. 0 ) THEN
          IF( sort .NE. 'shot' ) THEN
              IF( lbuf(51) .LT. 0 .AND. lbuf(7) .LT. mintrs ) THEN
                  lbuf(51) = 0
                  ndead = mintrs - lbuf(7)
              ENDIF
          ENDIF
      ENDIF           

      RETURN
      END
