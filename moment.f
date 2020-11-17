      SUBROUTINE moment(data,n,ave,adev,sdev,var,skew,curt)
!****  From Numerical Recipes
!  Given an array of DATA of length N, this routine returns its mean AVE, 
!  average deviation ADEV, standard deviation SDEV, variance VAR,
!  skewness SKEW, and kurtosis CURT.
      INTEGER n
      REAL adev,ave,curt,sdev,skew,var,data(n)
      INTEGER j
      REAL p,s,ep

      IF(n.le.1) PRINT*,' n must be at least 2 in moment'
      s=0.
      do 11 j=1,n
        s=s+data(j)
11    continue
      ave=s/n
      adev=0.
      var=0.
      skew=0.
      curt=0.
      ep=0.
      do 12 j=1,n
        s=data(j)-ave
        ep=ep+s
        adev=adev+abs(s)
        p=s*s
        var=var+p
        p=p*s
        skew=skew+p
        p=p*s
        curt=curt+p
12    continue
      adev=adev/n
      var=(var-ep**2/n)/(n-1)
      sdev=sqrt(var)
      if(var.ne.0.)then
        skew=skew/(n*sdev**3)
        curt=curt/(n*var**2)-3.
      else
         PRINT*,' no skew or kurtosis when zero variance in moment'
      endif
      return
      END
