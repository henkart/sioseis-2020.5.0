      subroutine ldeogeom(lhead,ihead,itype,ildgounit,xbinp,ioff0p,
     &                    nfoldp,ngpspp,offsetp,offset2p)
!*GMK
!*** Modifications:
!*** Mod 23 Oct 1992 by pch to start at the beginning of the file on each
!***            shot.  Also revert back to checking for shot number and
!***            shot time being equal,  Sepr line 27 has duplicate shot
!***            numbers.  The line also has dropped shots.
!*** Mod 16 Dec 1992 by GMK (ESPs distance d is actually in 1/10's of meters
!*              d = d/10. Also kill ghost shots THISFLAGSHOT = 3. Also kill 
!*               shots w/ bougus bearing if type = 5 (ESP).
!*   Mod 17 Dec 1992 Had a mix-up w/ OS bearing and MS bearing in code -- fixed.
!*   Mod 18 Dec 1992 Place LDGO line # into CDP slot in trace header and fake some CDP trace
!*              #'s to sort on such that the CDP (really ESP) can be gathered.
!*   Mod 15 Apr 1994 ifudgemin so as to allow 1 sec mismatches near minute mark
!*               Oh well, time to CHECK only SHOt #'s  for comparison, everything else is
!*               unreliable over an entire experiment.
!*   Mod 21 Dec 1994 Basically only shot # reliable so...check nav file for shot # only
!*               but will warn user if any mismatch occurs, Yikes!
!  mod July 95 by pch
! change ifudgehour = abs(tthr-ihr)  to ifudgehour = abs(tthr-ihour)
! and ifudgeminute to ifudgemin
!  18 Mar 1997 - add the line3 option so that we don't kill any traces
!                when line3 = 1.  Also set    thisshotflag = 1  when its 2!
!  30 Apr 97 - If type .EQ. -3, then force every shot in the nav file to
!              be 1 (thisshotflag = 1)
!  June 2000 - g77 didn't like float(ilogcdpdis) because ilogcdpdis is a real
!                              float(ioffset)
!  June 2000 - g77 didn't line line3 being a logical
!
!
!*    The following code was modified from a LDGO program that calculates offset(range)
!*    for various geometries, namely type = 3 (CDP), type = 4 (WAP) and type = 5 (ESP).
!*    These parameters can be calculated assuming one is given a LDGO log file and some
!*    other geometrical information (i.e. gun-antenna offset for both ships etc.....).
!*    To find the proper shot in the log file, it is neccessary to know the field shot
!*    number (bytes 245-248 BCD) for trace 0 header which is read during the input process.
!*    The following variables are returned from the log file (binary) and are used to 
!*    calculate offset and CDP number (type = 3 or 4):
!*    Navigation file:
!*        info header:   position (4byte words) total size = 1024 bytes
!*                       idirshot1 = 2        1st shot number
!*                       idirshotn = 3        last shot number
!*                       inumrec   = 7        total number of shots
!*                       idirangfact = 19     angle factor (i.e. angle = angle/idirangfact
!*                       ildgoline = 20       LDGO line number
!*
!*        shot headers:  position (4byte words) total size = 256 bytes * number of shots
!*                       msshot = 2            field shot number found in lbuf(3)
!*                       thisshotflag = 5      = 1 if nav ship (i.e Ewing), or = 2 if other ship 
!*                       msdeadshot = 6        bad trace, nav ship
!*                       osdeadshot = 7        bad trace, other ship
!*                       ttjd = 9              julian day >
!*                       tthr = 10             hour        >  >  >
!*                       ttmn = 11             min          >  >  use to verify shot
!*                       ttsec = 12            sec           >      an additional check!
!*                       finalrange = 16       ship to ship range (tenths of meteres)
!*                       mscourse = 19         nav ship course
!*                       oscourse = 20         other ship course
!*                       bearostoms = 21       bearing other ship to nav ship
!*                       distcdp1 = 22         distance along nav. line from 5000 m reference point
!*                       mslat = 23            latitude of nav ship (MS)  deg*10**6
!*                       mslon = 24            longitude of nav ship (MS) deg*10**6
!*                       mscdpno = 25          cdp number corresponding to this shot MS
!*                       spcdpbin = 30         cdp bin size meters (*100)
                                        
      integer*4 logbuf(256)
      real*4    buf(256)
      EQUIVALENCE ( buf(1), logbuf(1) )
      integer*4 ishot, ishotold, itrace, lhead(111),   
     &idirshot1,idirshotn, inumrec, idirangfact, ishotpos, numtrys,
     &maxtrys,ildgounit, ishotmis, msshot, thisshotflag, msdeadshot,
     &osdeadshot, finalrange, mscourse, oscourse, bearostoms,
     &mslat, mslon, msmaxcdp, spcdpbin, ttjd, tthr, ttmn, ttsec, type,
     &ilogssh, ilogrsh, ilogstorb, ilogrange, ixsdist , ktrace, istrv,
     &icdpno, nfold, ildgoline

      real*4 xbin, ioff0, ngpsp, dag1, dag2, d, theta, phi, ymd,
     &xmd, yda, xda, xxmt, xat, yat, angfact, radeg, range, dist,
     &rsshp, refdist, distcdp1, ngpspp, ilogcdpdis, ioffset

      integer*2 iday, imin, ihour, isec, ihead(111)

      logical first, first1
!****  line3 on tera has dead trace flag on everything, don't kill any
!****  shots when line3 is set to 1, only kill when line3 is 0
      data first /.true./, maxtrys /10/, first1 /.true./, line3/0/

      save
                    
      type = IABS(itype)
!*    Get variables from the lhead(lbuf) or ihead(ibuf) arrays regarding shot trace                    
      ishot  = lhead(3)
      itrace = lhead(4)             
!*    If same shot, different trace - then skip past log file to savetime
      if (ishot.eq.ishotold) goto 99
      iday   = ihead(80)
      ihour  = ihead(81)
      imin   = ihead(82)
      isec   = ihead(83)

!*    check to see if shot number is zero, if so kill it and exit subroutine
      if (ishot.eq.0 .AND. line3 .EQ. 0 ) then
        ihead(15) = 2
        goto 999
      endif

!*    search through log file to find shot and compare the GMT to assure
!*    that the proper trace is identified

!*    get line variables         
      if (first) then
         call podisc( ildgounit, 1, 0)                                  ! rewind
         call rddiscb( ildgounit, logbuf, 1024, istat)
         if (istat .NE. 1024 ) THEN
             PRINT *,' ***  ERROR  ***    Bad ldgo nav file.'
             stop
         ENDIF
!*       get variables form line log header 1024 bytes long
         idirshot1 = logbuf(2)
         idirshotn = logbuf(3)
         inumrec   = logbuf(7)    
         idirangfact = logbuf(19)
         ildgoline = logbuf(20)
!*       disable if loop since we don't need info for each shot
!         first = .false.
      endif
                        
!*    look for shot # in log file
!*    skip 1024 bytes (255 4byte words) and nshots*256 bytes (64 4 byte words)  --> ishotpos
!*    if shot = 0 i.e. lhead(3)=0, then go to mute trace end then exit - bad trace
   10 CONTINUE
!*    read nav shot values
      call rddiscb( ildgounit, logbuf, 256, istat )
      if (istat .NE. 256 ) THEN
         PRINT *,' ***  ERROR  *** Could not find shot ',lhead(3),
     &        ' trace ',lhead(4),' in the ldgo nav file.'
         ihead(15) = 2
         RETURN
      ENDIF
!*    set  match variables from nav header and compare to values read in from trace 0
      msshot = logbuf(2)                                                !nav file shot num
      ttjd   = logbuf(9)                                                !time julian day
      tthr   = logbuf(10)                                               !time hour
      ttmn   = logbuf(11)                                               !time minute
      ttsec  = logbuf(12)                                               !time sec
      ifudgesec = abs(ttsec-isec)                                       !sometimes shot tme trace 0 and logfile shots
                                                                        !are off by a second - hence its fudged to be within
                                                                        !1 second of shot time
      ifudgeday = abs(ttjd-iday)                                        !day are sometimes off by 1 too
      ifudgemin = abs(ttmn-imin)                                        ! minute will be off, if off by a sec. near min. mark
      ifudgehour = abs(tthr-ihour)                                        ! check hour mark also


!*    Did we find the proper shot in nav file

      if ( msshot.eq.lhead(3) ) then
!* Lets check for fun and see how off the DSS-240 vaules are:
!* Write vaules to screen for trace 1 only
      if ( itrace .eq. 1) then
       if (ifudgeday.ne.0) print*, '***Warning, day mismatch of ',
     &  ifudgeday, ' with navigation file for shot#', msshot
       if (ifudgehour.ne.0) print*, '***Warning, hour mismatch of ',
     &  ifudgehour, ' with navigation file for shot#', msshot
       if (ifudgemin.ne.0) print*, '***Warning, minute mismatch of ',
     &  ifudgemin, ' with navigation file for shot#', msshot     
       if (ifudgesec.ne.0) print*, '***Warning, second mismatch of ',
     &  ifudgesec, ' with navigation file for shot#', msshot
      endif     
     
!*   Lamont timing mark on DSS-240 sucks the big one -- not consistent
!*   we will risk for now, not using another check since its not robust enough
!*   will use fudge parameters for now
!*   if so, read variables for particular shot from 256 byte nav. header 
         thisshotflag = logbuf(5)                                       !=1, nav; = 2, os
         IF( itype .EQ. -3 ) thisshotflag = 1
         msdeadshot = logbuf(6)                                         != -1 bad nav shot
         osdeadshot = logbuf(7)                                         != -1 bad os shot
         finalrange = logbuf(16)                                        !ship to ship range
         mscourse = logbuf(19)                                          !nav ship course
         oscourse = logbuf(20)                                          !other ship course
         bearostoms = logbuf(21)                                        !bearing other ship to mother ship
!*GMK now read as real, distcdp         
         distcdp1 = buf(22)                                             !distance to initial cdp nav point
         mslat = logbuf(23)                                             !lat of shot nav ship *10**6
         mslon = logbuf(24)                                             !long of shot nav ship *10**6
         msmaxcdp = logbuf(27)                                          !max num cdps this shot nav ship
         spcdpbin = logbuf(30)                                          !cdp bin size *100
!         print *,(ii,logbuf(ii),ii=1,30)
      else
!*     Basically looks at each value in the LDGO log file until a match is achieved,
!*     or end of file is encountered
         goto 10
      endif

!*    check if type = 3 (cdp) and os ship shots were recorded, then kill os shots or..
!*    check if type = 5 (esp) and ms ship shots were recorded, then kill ms shots
!*    and place bougus range in header (-999999)
99    if (type.eq.3 .and. thisshotflag.eq.2 ) then
         IF( line3 .EQ. 0 ) THEN
             ihead(15) = 2
             lhead(10) = -999999
             goto 999
         ELSE
             thisshotflag = 1
         ENDIF
      else if (type.eq.5 .and. thisshotflag.eq.1.AND. line3 .EQ. 0) then
         ihead(15) = 2
         lhead(10) = -999999
         goto 999
      endif

!*    check to see if bogus bearing is placed in either logbuf(19,20,21)
!*    if so, kill trace and place bougus range in header (-999999)
       if (type .ne. 3 .AND. line3 .EQ. 0 ) then
        if (mscourse.eq.-999999 .or. oscourse.eq.-999999 .or.
     &  bearostoms.eq.-999999) then
          ihead(15) = 2
          lhead(10) = -999999
          goto 999
        endif
       endif

!*    check to see if shot is really DSS - 240 trick shot to get the record length
!*    beyond 19.8 seconds - Will kill traces from ghost shot! If you want 0-40 sec.
!*    of data use SIOSEIS process CAT to concatenate shots (remember there will be
!*    a 2 second gap 19.8-20.0 and 39.8-40.0 DSS - 240 did not record then!
!*    will kill traces and place bogus range in header for now!
      if (type.eq.5 .and. thisshotflag.eq.3 .AND. line3 .EQ. 0 ) then
          ihead(15) = 2
          lhead(10) = -999999
          goto 999
      endif
           
!*    again check to see if there is a bad trace in log file
      if (thisshotflag.eq.1 .and. msdeadshot.eq.-1.AND.line3.EQ.0) then
         ihead(15) = 2
         goto 999
      elseif (thisshotflag.eq.2.and.osdeadshot.eq.-1.AND.line3.EQ.0)then
         ihead(15) = 2
         goto 999
      endif                   

!*    we are now going to do the offset and cdp number calculations for type = 3,4 and 5
!*    1st, we define the navigation and shooting parameters from the ldgo geometry file 
!*    via sioseis (for both the shooting and navigation ships)
         
      if (first1) then
	xbin  =	abs(xbinp)                                                   !bin size in meters for cdp.
	ioff0 = abs(float(ioff0p))                                           !recording ship distance antenna to 1st receiver.)
        nfold = abs(nfoldp) 
!*GMK Group is now real                                                      !recording ship # of channels.
        ngpsp = abs(ngpspp)                                                  !recording ship group spacing.
!*      I've hardwired the istrv value to be 0 for the DSS 240 system
        istrv =	0                                                            !recording ship streamer channel 1, 0-far, 1-near
        dag1 = abs(offsetp)                                                  !navigation ship distance antenna to guns.
        dag2 = abs(offset2p)                                                 !other ship ship distance antenna to guns.
        first1 = .false.
      endif

!*    define the following from appropriate nav log record knowing the shooting and
!*    recording ships.
        ilogssh   = oscourse                                            !shooting ship course
        ilogrsh   = mscourse                                            !recording ship course
        ilogstorb = bearostoms                                          !bearing os to ms.
        ilogrange = finalrange                                          !ship to ship range this shot

      ixsdist = 0
!*    ixdist used for explosive ESP shots - not done today
	if(type.eq.5) then
           ixsdist = 0                                                   !shooting ship extra source distance, set = 0
	endif

!*    reverse bearing to ms to os ship - used for two ship, two streamer work!
!*    recent experiments do not include this geometry - nevertheless is included for future work
	if(thisshotflag.eq.1) then
           ilogstorb=ilogstorb-180                                      !we really want bearing shooting to recording.
	endif                                                

!*    distance from initial navigation point to shot location
          ilogcdpdis = distcdp1                                         !distance along line ms this shot.

!*    offset calculation.
        radeg=3.14159/180.
!*    angfact is used to retrieve more precise angle measurements (usually 100)
        angfact = float(idirangfact)
        if(angfact.eq.0) angfact=1.0
        angfact=1./angfact
        if (thisshotflag.eq.2) then
           d=float(ilogrange)/10.0                                      !ship ship distance (tenths of meters)
           theta=angfact*float(ilogstorb-ilogrsh)*radeg                 !bearing angle between shooting ship and nav. ship
           phi=angfact*float(ilogssh-ilogrsh)*radeg                     !feather angle of shooting ship w/ respect to nav. ship
           ymd=d*sin(theta)                                             !perpindicular projection w/ respect to nav. ship of s-s distance
           xmd=d*cos(theta)                                             !parallel projection w/ respect to nav. ship of s-s distance
           yda=(dag2+float(ixsdist))*sin(phi)                           !perpindicular projection w/ respect to nav. ship of source feathering
           xda=(dag2+float(ixsdist))*cos(phi)                           !parallel projection w/ respect to nav. ship of source feathering
       else
           ymd=0
           xmd=0
           yda=0
           xda=dag1                                                     !airgun to antenna distance nav ship = dag1
       endif
!*    calculate offset w/ respect to streamer channel
        ktrace=itrace
        if(istrv.eq.1) then                                             !streamer channel 1 is near trace, not the case however
           ktrace=nfold-ktrace+1
        endif
        xxmt=ioff0+float(nfold-ktrace)*ngpsp                            !distance from antenna to itrace channel on streamer
        xat=xxmt-xmd-xda
        yat=-ymd-yda
        range=sqrt(xat**2+yat**2)                                       !calculate range using X**2 + Y**2 = Z**2
        if(xat.gt.0) then
           range=-range
        endif
        if(type.ne.5) then                                              !ranges absolute for waps and cdps
           range=abs(range)
        endif

!*      place integer range in segy header
!*GMK Make ioffset real
        ioffset= range                                                  !make range an integer
        lhead(10) = nint(ioffset)

!*    calculate cdp number (LDGO line number) and psuedo cdp trace number for
!*    ESP binning of ranges i.e. 25 m range increments will have same trace #
!*    which is preferable for SIOSEIS process SORT

      if (type.eq.5) then
        lhead(6) = ildgoline
        lhead(7) = nint(10000 + (ioffset/25.0) )
      endif
        
!*    calculate cdp number corresponding to log (type 3 and 4)
          
      if (type.eq.3 .or. type.eq.4) then
                    
	rsshp=0.0

!*    if the shooting ship is not the navigation ship, i.e. thisshotflag = 2  
	if(thisshotflag.eq.2) then
	  rsshp= float(ilogrange)/10.0
	endif
                                                              
!*    if the shooting ship is not the lead ship, i.e. not 1 (2 source, 2 receiver geometry)
	if(thisshotflag.eq.2) then
	  rsshp=-rsshp
	endif

!*    calculate half offset variable
	refdist=ioffset*.5          

!*    if the shooting and recording ships are the same or the shooting ship
!*    is the lead ship(1)
	if(thisshotflag.eq.1) then
	  refdist=-refdist
	endif
!*	 print *,' refdist=',refdist,' ilogcdpdis=',ilogcdpdis,
!*     &   ' rsshp=',rsshp,' dag1=',dag1

!*    midpoint distance away from navigational initial point
        if (thisshotflag.eq.1) then
          dist=ilogcdpdis+rsshp-dag1+refdist
        else
          dist=ilogcdpdis+rsshp-dag2+refdist
        endif 

!*    calculate cdp number and place it in segy header
!*	print *,' dist=',dist,' ilog=',ilogcdpdis
	icdpno = int((dist+xbin/2.0)/xbin + 1)
        lhead(6) = icdpno

!*    endif to calculate cdp
      endif
 
!*    put lat, lon, dead trace into segy header,
!*    and ready it to be passed onto next process
  
      lhead(19) = mslon                                                  !longitude*10**6
      lhead(20) = mslat                                                  !latitude*10**6
      ihead(15) = 1                                                      !good trace

999   ishotold = ishot

      return
      end

