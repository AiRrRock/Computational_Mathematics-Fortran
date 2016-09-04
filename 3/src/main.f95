SUBROUTINE fun (T,Y,YP)
      REAL T,Y(2),YP(2)
       YP(1) = -40*Y(1)+260*Y(2)+1/(10*t*t+1)
       YP(2) = 30*Y(1)-270*Y(2) +exp(-2*t)
      RETURN
END
program Lab3      
      EXTERNAL fun
      INTEGER :: NEQN=2
      REAL :: T=0.0,Y(2),TOUT,RELERR=0.1E-03,ABSERR=0.0
      REAL :: TFINAL=0.4,TPRINT=0.02,WORK(27)
      INTEGER IWORK(5),IFLAG
      Y(1)=0
      Y(2)=1
      IFLAG=1
      TOUT=T
      call rungekutta(0.01)
      call rungekutta(0.001)
      print* ,"RUNGE-KUTTA-FEllBERG 45"
   10 CALL RKF45(FUN,NEQN,Y,T,TOUT,RELERR,ABSERR,IFLAG,WORK,IWORK)
      PRINT 11,T,Y(1),Y(2)
      GO TO (80,20,30,40,50,60,70,80),IFLAG
   20 TOUT=TPRINT+T
      IF(T.LT.TFINAL) GO TO 10
      STOP
   30 PRINT 31,RELERR,ABSERR
      GO TO 10
   40 PRINT 41
      GO TO 10
   50 ABSERR=0.1E-07
      PRINT 31,RELERR,ABSERR
      GO TO 10
   60 RELERR=RELERR*10.0
      PRINT 31,RELERR,ABSERR
      IFLAG=2
      GO TO 10
   70 PRINT 71
      IFLAG=2
      GO TO 10
   80 PRINT 81
      STOP
   11 FORMAT(' T=',F7.4,2X,'Y1=',F10.6,2X,'Y2=',F10.6)
   31 FORMAT(' ÃPAHÈÖÛ ÏOÃPEØHOCTEÉ ÈÇMEHEHÛ  ',' RELERR=',E10.3,2X,'ABSERR=',E10.3)
   41 FORMAT(' MHOÃO ØAÃOB ')
   71 FORMAT(' MHOÃO BÛXOÄOB ')
   81 FORMAT(' HEÏPABÈËÜHÛÉ BÛÇOB ')
contains
subroutine RungeKutta(h)
REAL    K0(NEQN),K1(NEQN),K2(NEQN),YN(NEQN),Ytmp(NEQN),Y(NEQN)
REAL    :: h, t=0.0
integer :: counter=0
      T=0.0
      Y(1) =0
      Y(2) =1.0
print*, "Runke-Kutta the 3rd with the step=",h
157   Ytmp=Y
      call fun(t, Ytmp,YN)
      K0 = h*YN
      Ytmp = Y +(K0/2)
      call fun (t+(h/2.0),Ytmp,YN)
      K1= h*YN
      Ytmp= Y +((2*K1)-K0)
      call fun (t+h,Ytmp,YN)
      K2 = h*YN
      YN = Y +1.0/6.0*(K0+4*K1+K2)
      Y =YN
      T=h+T
      counter = counter +1 
      IF ((h.LE.0.001) .and. mod(counter,20).eq.0) PRINT 11 ,T,YN(1),YN(2)
      IF ((h.GE.0.01).and. mod(counter,2).eq.0) PRINT 12, T,YN(1),YN(2)
      IF(T.LT.0.4-h) GO TO 157

11 FORMAT(' T=',F7.4,2X,'Y1=',F15.6,2X,'Y2=',F15.6)
12 FORMAT(' T=',F7.4,2X,'Y1=',ES15.6,2X,'Y2=',ES15.6)

end subroutine RungeKutta      
end program Lab3
