      logical function chkprc()
!-------------------------------------------------------------------------------
! chkprc checks to see if there if the user has specified an FK process in the
! procs list.
!   Outputs:
!      .TRUE.  if an FK process is specified.
!      .FALSE. if no FK process specified
!
!   Last Modified:
!      10/7/88
!      6 Sep 93 by pch - increase pnames and fkproc to be character*7
!-------------------------------------------------------------------------------
!
!.. The following are needed to check the user processes (set in GETPRO)
      parameter ( NAMESS  = 47)
      integer pnum                                                      ! The number of proceeses called
      character*7 pnames(NAMESS)
      common /porder/ pnum, iorder(NAMESS)
      common /pnames/ pnames                                            ! The names of the processes called
!
      parameter ( NFKPROC = 4)
      character*7 fkproc(NFKPROC)                                       ! List of known FK processes
      logical    found
!
      data  fkproc / 'FK2TX ','FKMIGR','FKFILT','TX2FK'/
!
! Check the Processes in the procs list for an FK process
      i = 1
      found = .FALSE.
10    if ( (i.le.pnum) .AND. .not.found) then
        j = 1
15      if ( (j.le.NFKPROC) .AND. .not.found) then
          if (pnames(i).eq.fkproc(j) ) found = .TRUE.
          j = j + 1
          GO TO 15
        endif
        i = i + 1
        GO TO 10
      endif
!
      chkprc = found
      return
      end
