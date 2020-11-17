      subroutine FK2TED(scr, lscr )
!                              Process FK2TX
!                              ------- -----
!-----------------------------------------------------------------------------
!  Document date: 17th October 1988    version 1.0    a.j. harding
!
!     Process FK2TX transforms data back from the FK (frequency-wavenumber)
!  domain into the TX (time-space) domain.  There is, in general, no need to
!  specify any parameters to FKT2X since any FK domain processes are performed
!  by separate processes such as FKFILT and the parameters needed for the back
!  transformation, such as the number of wavenumbers, are stored in either the
!  binary header or the trace header.
!
!     However PATH1 & PATH2 parameters like those in TX2FK may be provided so
!  as to control the disk location of scratch datasets used by the program.
!  These are useful if there is insufficient space in the default locations
!  used by SIOSEIS.
!
!  The Parameter Dictionary
!  --- --------- ----------
!  PATH1  - This is where input FK traces are accumulated prior to back
!            back transformation. It is deleted after back transformation, and
!           before the output tx data is written to disk. Thus it is usually
!           safe to put this file in the same directory as the final output. If
!           procedure TX2FK is present in the procs list then this scratch
!           file is the same as the 2nd scratch file of TX2FK. Care should be
!           taken not to specify this file twice. Only the first filename wil
!           be used by SIOSEIS.
!           Default: Implementation dependent.
!
!  PATH2  - The name of the second scratch file to be used by SIOSEIS.
!           Default: Implementation dependent.
!
!  IHDRPATH - Filename containing a set of original TX trace headers that were
!           written by process TX2FK. These headers will be added to the TX
!           traces after transformation. On the assumption that a 2-D process
!           was performed in F-K all traces will be marked live and all mute
!           entries will be zeroed in the header.
!
!  OPAD   = YES. All padding both in time and range added to data prior to
!           2-D FFT will be output. Preset NO.
!
!  LPRINT -  The debug parameter.
!
!  END    - terminates each parameter list.
!
! Modifications
!    3/14/89 ajh : Converted to 6 letter externals. And spaces at start of print
! Last Modified
!    6/21/89 ajh : Added input parameters ihdrpath, opad. Changed call for
!                  fk2ted. Pass variables via disk file.
!    31 Aug 05 - Add parameter KILL
!    15 Aug 07 - g95 didn't like doing an internal read of ihdrpath
!-------------------------------------------------------------------------------
!
! INCLUDE FILES
! INCLUDE seisio.inc
! This include file include all (?) of the parameters needed to define SIOSEIS
! I/O routines. In particular it defines positions of objects in the SEGY
! headers & constants for these values
!
! EBCDIC HEADER
! Declare Header lengths. These are the length of the EXTERNAL disk image in
! host words.
      integer     CRAYEBC, NORMEBC
      parameter ( CRAYEBC  = 400)
      parameter ( NORMEBC  = 2*CRAYEBC)
!
! BINARY HEADER
      integer     CRAYBIN, NORMBIN
      parameter ( CRAYBIN  =  50)                                       ! length of external disk image
      parameter ( NORMBIN  = 2*CRAYBIN)
!
!  The position of objects in the Binary header
      integer    NTRCPTR,IFMTPTR, IDTYPPTR, INKPTR, ITSIPTR
      integer    ITDELPTR, IVERSPTR
      parameter (NTRCPTR  =  7)                                         ! No of traces / gather
      parameter (IFMTPTR  = 13)                                         ! The SEGY Format of the data
      parameter (IDTYPPTR = 31)                                         ! The Domain/ID of the data e.q. T-X data
      parameter (INKPTR   = 32)                                         ! The No. of wavenumbers in an f-k dataset
      parameter (ITSIPTR  = 33)                                         ! Record of tx sample interval in us
      parameter (ITDELPTR = 34)                                         ! Record of tx time delay in ms.
      parameter (IVERSPTR = 35)                                         ! The SIOSEIS version number
!
!..Constants for SEGY Format (IFMTPTR)
      integer IBMFP, INT16, INT32, HOSTFP
      parameter (IBMFP  = 1)                                            ! IBM floating point
      parameter (INT32  = 2)
      parameter (INT16  = 3)
      parameter (HOSTFP = 5)                                            ! Host FP. In actuality Host FP >= 5
!
!.. Data Domains/Types (IDTYPPTR)
      integer    IDTX, IDFKRCT, IDFKPLR, IDFKPLRU
      parameter (IDTX     = 1)
      parameter (IDFKRCT  = 2)
      parameter (IDFKPLR  = 3)
      parameter (IDFKPLRU = 8)
!
!.. Version Number
      integer  CURVERS
      parameter( CURVERS = 210)
!
! TRACE HEADER
!.. Declare header length on external disk files. Internal length is given by
! Numhdr in common block READT.
      integer  CRAYTHDR, NORMTHDR
      parameter ( CRAYTHDR =  30)                                       ! Length of external disk image
      parameter ( NORMTHDR = 2*CRAYTHDR)
!
!    The positions of elements of the trace header are defined by elements of
! the common block SEGYPTR. This block is initialized in routine SETPTR. The
! position of elements differs on the CRAY since it does not allow mixing of
! Integer*16 & Integer*32 words but only has Integer*64.
!
      integer LLSEQPTR                                                  ! Trace Sequence No. within line
      integer LRSEQPTR                                                  ! Trace Sequence No. within reel
      integer LSHOTPTR                                                  ! Shot number or Stacked Number
      integer LSHTRPTR                                                  ! Trace number within Shot
      integer LRPNPTR                                                   ! RP or CDP number
      integer LRPTRPTR                                                  ! Trace No. within CDP
      integer ITRIDPTR                                                  ! Trace ID Live(1)/Dead(2)
      integer LDISTPTR                                                  ! Source to Receiver Distance
      integer LWBDPTR                                                   ! Water Bottom depth at source
      integer LSXCOPTR                                                  ! Source X co-ordinate
      integer LRXCOPTR                                                  ! Receiver X co-ordingate
      integer IDELMPTR                                                  ! Deep water delay in ms.
      integer ISTMPTR                                                   ! Start Mute time in ms.
      integer IENDMPTR                                                  !   End Mute time in ms
      integer ISAMPPTR                                                  ! Number of data samples in trace
      integer ISIPTR                                                    ! Sample interval of trace in us
      integer IYRPTR                                                    ! Year data was recorded
      integer IDAYPTR                                                   ! Day of year
      integer IHRPTR                                                    ! Hour of day
      integer IMINPTR                                                   ! Minute of hour
      integer ISECPTR                                                   ! Second of Minute
      integer IGMTPTR                                                   ! Time zone of header Local(1)/GMT(2)
      integer LDELSPTR                                                  ! Deep water delay in secs.
      integer LSMUSPTR                                                  ! Start of mute in seconds
      integer LEMUSPTR                                                  ! End of mute in seconds
      integer LSISPTR                                                   ! Sample interval in seconds
      integer LWBTSPTR                                                  ! Water bottom time in seconds
      integer LGATPTR                                                   ! No. of traces in stacked data(>0) or End of Gather(<0)
      integer LSSMSPTR                                                  ! Start of surgical mute in secs.
      integer LESMSPTR                                                  ! End of surgical mute in secs
      integer LSBPTR                                                    ! SeaBeam slant range
!
      common /SEGYPTR/ LLSEQPTR, LRSEQPTR, LSHOTPTR, LSHTRPTR, LRPNPTR,
     *                 LRPTRPTR, ITRIDPTR, LDISTPTR, LWBDPTR,  LSXCOPTR,
     *                 LRXCOPTR, IDELMPTR, ISTMPTR,  IENDMPTR, ISAMPPTR,
     *                 ISIPTR,   IYRPTR,   IDAYPTR,  IHRPTR,   IMINPTR,
     *                 ISECPTR,  IGMTPTR,  LDELSPTR, LSMUSPTR, LEMUSPTR,
     *                 LSISPTR,  LWBTSPTR, LGATPTR,  LSSMSPTR, LESMSPTR,
     *                 LSBPTR
!
!.. Define Constants for Trace Header
!
!.. Constants for Trace ID (ITRIDPTR)
      integer    LIVETR, DEADTR                                         ! Parameters for trace Type
      parameter (LIVETR  = 1)
      parameter (DEADTR  = 2)
!
!.. Define constants for use with the DISKIO subroutines getfil, frefil,
! podisc
! GETFIL constants
      integer    CREATTMP, CREATNEW, OPENOLD
      parameter (CREATTMP   = 1)                                        ! Create a temporary file
      parameter (CREATNEW   = 3)                                        ! Create a new named file
      parameter (OPENOLD    = 4)                                        ! Open an existing file
! FREFIL constants
      integer    RLDLCLS, RLCLS
      parameter (RLCLS   = 2)                                           ! Release and Close the file
      parameter (RLDLCLS = 3)                                           ! Release, delete and close
! PODISC constants
      integer    POSABS, POSREL, DSKSTRT
      parameter (POSABS     = 1)
      parameter (POSREL     = 2)
      parameter (DSKSTRT    = 0)                                        ! The starting position on the disk
! INCLUDE txfkgbl.inc
!.. This is the global common block for the FK routines. At present it is used
!.. by TX2FEX, FK2TEX & their associated edit routines to indicate whether they
!.. have been called & to pass file unit numbers.
!
!  Last Modified 10/18/88
!
      logical txinit, fkinit                                            ! Execution initialization ?
      logical txed, fked                                                ! Edit initialization ?
      integer ltxunt1, ltxunt2                                          ! File Stream ID for tx2fex scratch file
      integer lfkunt1, lfkunt2                                          ! File Stream ID for fk2tex scratch file
      integer txunit                                                    ! Used Stream ID for tx2fex scratch file
      integer ohdrtxfk                                                  ! File unit holding tx headers
      logical tmptxhdr                                                  ! Indicates if tx header file is temporary
      integer txprestk
      common /TXFKE/ txinit, fkinit, txed, fked, txunit,
     $               ltxunt1, ltxunt2, lfkunt1, lfkunt2,
     $               ohdrtxfk, tmptxhdr, txprestk, ntx2fk, nfk2tx, range
! INCLUDE txfk.inc
! Include file TXFK.INC
! include file for the sioseis 2D FFT and fk migration routines
!
!
! Parameters for coord type
      integer CDRECT, CDPOLAR, CDPOLARU
      parameter ( CDRECT   = 1)
      parameter ( CDPOLAR  = 2)
      parameter ( CDPOLARU = 3)
!
! parameters for Windows
      integer BoxCar
      parameter (BoxCar = 5)
! INCLUDE siogbl.inc
!... Computer Type
      integer PRIME, VAXUNIX, APOLLO, VAXVMS, CRAY
      integer CONVEX, IEEE
      parameter (PRIME   = 1)
      parameter (VAXUNIX = 2)
      parameter (APOLLO  = 3)
      parameter (VAXVMS  = 4)
      parameter (CRAY    = 5)
      parameter (CONVEX  = 6)
      parameter (IEEE    = 7)
!
! PROGRAM
      real     scr(111)
      integer lscr(111)
!
      parameter (NPARS  = 6)                                            ! the number of user parameters
!
!
      integer parmunit, nwrds
      common /FK2TXC/ parmunit, nwrds                                   ! Pass parameters to execution phase.
!
      common /EDITS/ ierror,iwarn,irun,now,icompt
!
      character*8  names(NPARS)
      character*1  types(NPARS)
      integer     length(NPARS)
      character*100 token
!
      character*100 ihdrpath
      integer status, opad, kill
!
      data names/'LPRINT','OPAD','PATH1','PATH2','IHDRPATH','KILL'/
      data types/'L', 5*'A'/
      data length/6, 4, 2*5, 8, 4/
!
      lprint   = 0
      opad     = 0
      ihdrpath = ' '
      dead = 0
!..        set the presets in FK global common
      call txfkinit
      lfkunt1 = 0
      lfkunt2 = 0
      fked  =  .true.                                                   ! Let everybody else know we have been called
!..
!       The current command line in the system buffer may have the parameters.
!...                                Get a parameter list from the user.
  100 continue
      call getoke(token,nchars)                                         ! Get a token from the user parameter line
      call upcase(token,nchars)                                         ! Convert the token to uppercase
!
      if(nchars.le.0) then                                              ! Get another user parameter line
        if(now.eq.1) print 140
  140   format(' <  Enter Parameters  >')
        call rdline
        go to 100
      endif
!
      do 190 i = 1, npars                                               ! See if it is a parameter name
        len    = length(i)                                              ! Get the legal parameter name length
        iparam = i                                                      ! Save the index
        if(token(1:nchars).eq.names(i)(1:len).and.nchars.eq.len)
     *  go to 200
  190   continue
!                                                                       ! Still looking for the name
      if(token(1:nchars).eq.'END'.and.nchars.eq.3) go to 1000           ! end of list?
      print 191, token(1:nchars)
  191 format(' ***  Error  *** FK2TX does not have a parameter ',
     *  'named ',A10)
      ierror = ierror + 1
      go to 100
!
!..                                Found the parameter name, now find the value
  200 continue
      nparam  = iparam
  210   continue
        call getoke(token,nchars)
        if(nchars.le.0) then                                            ! End of line?
          if(now.eq.1) print 140                                        ! Allow parameters on a different line
          call rdline                                                   ! Get another line
          go to 210
        endif
!
  230 continue
      if( types(nparam).eq.'A') then                                    ! An alpha Parameter
!
        if ( names(nparam).eq.'PATH1') then
          if ( txed.and.(ltxunt2.ne.0) ) then
            lfkunt1 = ltxunt2
            print *,' *** Warning FK2TX **** First scratch file already'
     $             ,' specified by TX2FK. Using that file'
          else
            call getfil(CREATNEW,lfkunt1,token,status)
            if (status.ne.0) then
              print *,' ***ERROR*** Error opening PATH1: ',token
              ierror = ierror + 1
            endif
          endif
!
        else if ( names(nparam).eq.'PATH2') then
          call getfil(CREATNEW,lfkunt2,token,status)
          if (status.ne.0) then
            print *,' ***ERROR*** Error opening PATH2: ',token
            ierror = ierror + 1
          endif
!
        else if ( names(nparam).eq.'IHDRPATH') then
          ihdrpath = token(1:nchars)
          if ( ( icompt.eq.VAXUNIX ) .or. (icompt .eq. IEEE) )          ! UNIX machines
     $      ihdrpath(nchars+1:nchars+1) = char(0)                       ! Add null terminator
          if ( (icompt .eq. CRAY) .and. (nchars .gt. 8) ) then
             print *,' ***  ERROR  ***  Cray pathnames must',
     $               'be 8 or fewer characters long.'
             ierror = ierror + 1
          endif
        endif                                                           ! of if names == PATH1
!
        if (names(nparam) .eq. 'OPAD') then
          call upcase(token,nchars)
          if ( token(1:nchars) .eq. 'YES') then
            opad = 1
          else if (token(1:nchars) .eq. 'NO') then
            opad = 0
          else
            print *,' *** ERROR *** Illegal OPAD'
            ierror = ierror + 1
          endif
        endif
        IF( names(nparam) .EQ. 'KILL' ) THEN
            CALL upcase(token,nchars)
            IF( token(1:nchars) .EQ. 'YES' ) THEN
                kill = 1
            ELSEIF ( token(1:nchars) .EQ. 'NO' ) THEN
                kill = 0
            ELSE
                PRINT *,' ***  ERROR  ***  Illegal KILL.'
                ierror = ierror + 1
            ENDIF
         ENDIF
      endif                                                             ! of if  type == A
!
!
      if (types(nparam).eq.'L') then
        call dcode(token,nchars,areal,istat)                            ! try and decode it
        if(istat.ne.2) then                                             ! =2 means it is a numeric
          ierror = ierror + 1                                           ! dcode printed an error
          go to 100
        endif
        lprint = areal                                                  ! lazy method. Replace with lvals when more parameters
      endif
      go to 100
 1000 continue
!
      if (lprint.ne.0) then
        print *,' File stream 1 (lfkunt1): ',lfkunt1
        print *,' File stream 2 (lfkunt2): ',lfkunt2
      endif
!                                                                       ! Look for the end of the parameter list
 2020 call getoke(token,nchars)
      call upcase(token,nchars)
      if(nchars.eq.0) then                                              ! End of a line?. Then get another
        if (now.eq.1) print 140
        call rdline
        go to 2020
      endif
      if(token(1:nchars).ne.'END'.or.nchars.ne.3) then
        print *,' ***  ERROR  ***  FK2TX only permits one list to be ',
     *    ' given.'
        ierror = ierror + 1
      endif
!
!.. Write the parameter names to a disk file
      call getfil (CREATTMP, parmunit, token, istat)
      if (istat.ne.0) then
        print *,' *** TX2FK ERROR *** Opening parameter file'
        STOP
      endif
      call podisc (parmunit, POSABS, DSKSTRT)
!
!      read (ihdrpath, '(25A4)') (lscr(i), i = 1,25)
      CALL wrdisc( parmunit, ihdrpath, 25 )
      lscr(26) = opad
      lscr(27) = lprint
      lscr(28) = kill
      nwrds    = 28
      call wrdisc(parmunit, lscr(26), nwrds-25)
      if (lprint.ne.0) print *,' Exiting fk2ted'
      return
      end
