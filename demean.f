      SUBROUTINE DEMEAN(IN,IOUT,N,MEAN)
!     DEMEAN IS AN AP ROUTINE TO REMOVE THE MEAN VALUE OF A TRACE.  MEAN REMOVAL
! IS THE SAME AS REMOVING THE DC SHIFT OR DC BIAS.
!
! ARGUMENTS:
!  IN     - THE AP ADDRESS OF THE INPUT ARRAY (THE BIASED ARRAY).
!  IOUT   - THE AP ADDRESS OF THE OUTPUT ARRAY TO RECEIVE THE UNBIASED DATA.
!           (IOUT MAY BE THE SAME AS IN).
!  N      - THE NUMBER OF ELEMENTS IN THE IN ARRAY.
!  MEAN   - THE AP ADDRESS WHERE THE MEAN VALUE CAN BE STORED.
!
!
!  WRITTEN AND COPYRIGHTED BY:
!  PAUL HENKART, SCRIPPS INSTITUTION OF OCEANOGRAPHY, 27 OCTOBER 1981
!  ALL RIGHTS ARE RESERVED BY THE AUTHOR.  PERMISSION TO COPY OR REPRODUCE THIS
!  SUBROUTINE, BY COMPUTER OR OTHER MEANS, MAY BE OBTAINED ONLY FROM THE AUTHOR.
!
      CALL MEANV(IN,1,MEAN,N)                                           ! FIND THE MEAN VALUE
      CALL VNEG(MEAN,1,MEAN,1,1)                                        ! NEGATE THE MEAN VALUE
      CALL VSADD(IN,1,MEAN,IOUT,1,N)                                    ! SUBTRACT THE MEAN (ADD THE NEGATIVE)
      CALL VNEG(MEAN,1,MEAN,1,1)                                        ! MAKE THE MEAN AS IT WAS
      RETURN
      END
