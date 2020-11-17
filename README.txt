ver 2020.5.0 (14 July 2020)
   1) geom - Allow type 9 to use UTM coordinates (SEGY trace header integer 45)
   2) diskin - Use unsigned integer for numdat from the binary header in died.f
             - Print a warning if the SEGY trace id (ihdr(15)) is 0 (not set).
   3) diskio.c - put the tmp files back in the current directory since /tmp is NOT
                 cleared on boot as it used to be
   4) gather - Bad change in 2020.3.* when eliminating DO loop ending in label statement
               rather than CONTINUE.
   5) filter - Get nsamps from common rather than the trace header
             - Abort if more than 32767 samples
   6) contro, filters - Max MAXLEN=32767 so filters can hold 32767 complex numbers
   7) geom - Do not clobber shot coordinates on types 9, 18, 21
           - Add type 21, an ASCII nav shot number, decimal lat, decimal long
   8) diskin - Allow nsamps > 32767 for random access of fno
   9) plot - Make nsecs preset to what the doc says - the length of the first trace
   10) contro - Bad common/pltcom/ size
   11) plot - Changed COMMON /porder/ iorder(4) to same as contro iorder(100).  OSX cares.
   12) mute - take nsamps from common numdat rather than the trace header (>32767)
   13) geom - Detect bad NAVFIL under type 21
   14) velan - Variable ifin needed to be set 0 at the beginning
   15) plot - Add warning if ANN and ANN2 are both HEADER.
   16) geom - type 21 zero out the receiver lat/long since they aren' specified
   17) stack - Allow decimal degrees for lat/long (SEGY REV2 type 3)
             - Don't use receiver coordinates for midpoint calculation if they are 0
   18) prout - Check for too many INDICES - 10 max!
   19) gather - Check for too much data in the edit phase.
   20) geom - ABORT if negative CMP number detected
lsd - Honor coordinate type 1 (length or UTM)
