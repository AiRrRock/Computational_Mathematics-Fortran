program Lab2      
      INTEGER i,j,z
      integer, parameter ::NDIM=5, N=5 ,STEPS=3, REALNESS=4
      REAL(REALNESS) :: R(N,N),Original(N,N),NOTA(N,N),A(N,N),U(N,N),L(N,N),X,COND,IPVT(N),WORK(N),DET,B(N),E(N,N)
      do z=1,STEPS
         print*, "STEP",Z
         if (z==1) x=1.1
         if (z==2) x=1.001
         if (z==3) x=1.000001
         print*,"Value of X", X
         forall(i=1:N,j=1:N) A(i,j)=j
         forall (i=1:N,j=1:N,i==j) A(i,j)=1+x
         A(1,1)= 1
         Original=A
         print*, "Original Matrix"
         print 101, ((A(i,j),j=1,N),i=1,N)
         E=0
         forall(i=1:N,j=1:N,i==j)E(i,j)=1
         call DECOMP(NDIM,N,A,COND,IPVT,WORK)
         print*, "Conditional number of original Matrix", Cond
         U=0
         forall(i=1:N,j=1:N,j>=i) U(i,j)=A(i,j)
         L=0
         forall(i=1:N,j=1:N,i==j)L(i,j)=1
         forall(i=1:N,j=1:N,J<i) L(i,j)=A(i,j)
         forall(i=1:N,j=1:N,j<i) L(i,J)=-L(i,j)      
         !print 101, ((U(i,j),j=1,N),i=1,N)
         !print 101, ((L(i,j),j=1,N),i=1,N)
         !A= matmul(L,U)
         DET=IPVT(N)
         do i=1,N
            DET=DET*A(i,i)
         end do
         do i=1,N
            B=0
            B(i)=1
            call SOLVE(NDIM,N,A,B,IPVT)
            NOTA(1:N,i)=B
         end do
         print*, "Reciprocal Matrix"
         print 101,((NOTA(i,j),j=1,N),i=1,N)   
         R=matmul(Original,NOTA) 
         R=R-E
         print*, "Desired Matrix"
         print 101,((R(i,j),j=1,N),i=1,N)
         X=norm2(R)
         print*, "Norm of Matrix R"
         print 102,X 
         call DECOMP(NDIM,N,R,COND,IPVT,WORK) 
      end do 
         
 101  FORMAT((5ES14.7))
 102  FORMAT(1ES12.6) 
end program Lab2
