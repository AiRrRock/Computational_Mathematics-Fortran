
      REAL FUNCTION FUN(X)
      REAL X 
      FUN= 1/(1+25*X*X)
      END function fun 

      use Spline
      EXTERNAL FUN
     
      INTEGER NOFUN,i
      REAL FUN,A,B,RELERR,ABSERR,RESULT,ERREST,FLAG
      DATA A/-1./,B/1.0/,RELERR/1.E-06/,ABSERR/0.0/
      integer, parameter :: N=21 
      REAL FXSTUPID 
      REAL XAR(N) , FAR(N)
      REAL B1(N),C1(N),D1(N), integr, h
      do i=1,N
      XAR(i)=-1+(i-1)*0.1
      FAR(i)=FUN(XAR(I))
      end do
      
      CALL SPLIN(N,XAR,FAR,B1,C1,D1) 
      integr = 0
      DO i=1,N-1
         h = XAR(i+1)-XAR(i)
         integr = integr + FAR(i)*H - B1(i)/2*h**2 +C1(i)/6*h**3 - D1(i)/24*h**4
      END do
            
      FXSTUPID = 2*atan(5.0)/5  
      PRINT*, 'TYPE OF VALUE               ', ' Value  ','   Difference between computated value and oiginal value'
      PRINT*,  'Original value             ',    FXSTUPID     
      CALL QUANC8(FUN,A,B,ABSERR,RELERR,RESULT,ERREST,NOFUN,FLAG)
      Print*, 'Quanc8 value of function   ',    RESULT,'   ',   FXSTUPID-RESULT 
      !PRINT 1,RESULT,ERREST,NOFUN,FLAG
      CALL QUANC8(LAN,A,B,ABSERR,RELERR,RESULT,ERREST,NOFUN,FLAG)
      
      !PRINT 1,RESULT,ERREST,NOFUN,FLAG
      print*, 'Quanc8 val of Lagrunge     ',    RESULT,'   ',   FXSTUPID-result
      print*, 'Spline aproximation        ',    integr,'   ',   FXSTUPID-integr
      STOP
   
    contains 
    REAL FUNCTION LAN(X)  
    REAL X, pup, pdown
    integer j,i
    i=0
    LAN = 0
    do i=1, N
        pup=1
        pdown=1
        do j =1,N
           if (i /= j) then 
             pup = pup* (X -XAR(j))
             pdown = pdown* (XAR(i)-XAR(j))
           end if 
        end do
        LAn = LAn + FAR(i)*pup/pdown
    end do
    END function
    
 END
