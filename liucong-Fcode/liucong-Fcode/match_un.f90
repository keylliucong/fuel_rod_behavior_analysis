SUBROUTINE GAUSS(NT,A,X)
! Performing Gaussian elimination to solve a set of N simultaneous linear equations
IMPLICIT NONE 
real A,SAV,R,X 
INTEGER NT,IM1,K,L,I,J,M,MM 
DIMENSION A(NT,NT+1),X(NT) 
DO 20 I = 2,NT 
 DO 20 J = I,NT 
   IF (ABS(A(I-1,I-1)).LT.1.D-12) A(I-1,I-1)=0.D0 
         IF (A(I-1,I-1)) 1,2,1 
    2    IM1= I-1 
         DO 21 M =I,NT 
          IF (A(M,IM1)) 3,21,3  
    3     DO 22 MM = IM1,NT+1 
           SAV= A(M,MM) 
           A(M,MM)= A(IM1,MM) 
   22     A(IM1,MM)= SAV 
   21    CONTINUE 
    1   R= A(J,I-1)/A(I-1,I-1) 
        DO 20 K= I,NT+1 
   20  A(J,K)= A(J,K)-R*A(I-1,K)
    
DO 30 I= 2,NT 
      K= NT-I+2 
      R= A(K,NT+1)/A(K,K) 
  DO 30 J=I,NT 
      L= NT-J+1 
   30 A(L,NT+1) = A(L,NT+1)-R*A(L,K) 
      DO 40 I= 1,NT 
   40  X(I)= A(I,NT+1)/A(I,I) 
END SUBROUTINE GAUSS


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



subroutine agaus(a,b,x)
 implicit none
 integer,parameter:: n=4            !!!!!!!!!!!计算8阶的方程组
 integer::i,j,l,k,is
 real(8)::a(n,n),x(n),b(n),js(n)

 real(8)::d,t

   l=1      !逻辑变量,判断是否有解
   do k=1,n-1
     d=0.0
     do i=k,n
      do j=k,n
       if (abs(a(i,j))>d) then
         d=abs(a(i,j))
         js(k)=j
         is=i
       end if
   end do
     end do      !把行绝对值最大的元素换到主元位置
     if (d+1.0==1.0) then   
       l=0
     else     !最大元素为0无解
   if(js(k)/=k) then

    do i=1,n
      t=a(i,k)
   a(i,k)=a(i,js(k))
   a(i,js(k))=t
       end do       !最大元素不在K行，K行
   end if
   if(is/=k) then
    do j=k,n
     t=a(k,j)
  a(k,j)=a(is,j)
  a(is,j)=t     !交换到K列
       end do
    t=b(k)
    b(k)=b(is)
    b(is)=t    
   end if          !最大元素在主对角线上
     end if     !消去
  if (l==0) then
    write(*,100)
    write(*,*) "l==0"
    pause
    return
  end if
    do j=k+1,n
     a(k,j)=a(k,j)/a(k,k)
       end do
       b(k)=b(k)/a(k,k)     !求三角矩阵
    do i=k+1,n
     do j=k+1,n
    a(i,j)=a(i,j)-a(i,k)*a(k,j)
        end do
  b(i)=b(i)-a(i,k)*b(k)
       end do
     end do
     if (abs(a(n,n))+1.0==1.0) then
    l=0
    write(*,100)
    write(*,*) "l==0"
    pause
    return
  end if
  x(n)=b(n)/a(n,n)
  do i=n-1,1,-1
    t=0.0
    do j=i+1,n
      t=t+a(i,j)*x(j)
       end do
       x(i)=b(i)-t
     end do
100 format(1x,'fail')
    js(n)=n
 do k=n,1,-1
   if (js(k)/=k) then
    t=x(k)
    x(k)=x(js(k))
    x(js(k))=t
   end if
    end do
    return 
 end 




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







SUBROUTINE INVERT(A,AINV) 
! Performing invert a 3x3 matrix
IMPLICIT NONE 
real A,AINV,DETA 
INTEGER J,K 
DIMENSION A(3,3),AINV(3,3) 
  AINV(1,1)= A(2,2)*A(3,3)-A(3,2)*A(2,3) 
  AINV(1,2)= -(A(1,2)*A(3,3)-A(3,2)*A(1,3)) 
  AINV(1,3)= A(1,2)*A(2,3)-A(2,2)*A(1,3) 
  AINV(2,1)= -(A(2,1)*A(3,3)-A(3,1)*A(2,3)) 
  AINV(2,2)= A(1,1)*A(3,3)-A(3,1)*A(1,3) 
  AINV(2,3)= -(A(1,1)*A(2,3)-A(2,1)*A(1,3)) 
  AINV(3,1)= A(2,1)*A(3,2)-A(3,1)*A(2,2) 
  AINV(3,2)= -(A(1,1)*A(3,2)-A(3,1)*A(1,2)) 
  AINV(3,3)= A(1,1)*A(2,2)-A(2,1)*A(1,2) 
  DETA= A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))-A(1,2)*(A(2,1)*A(3,3)-A(3,1)*A(2,3)) &
      + A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2)) 
DO 10 J = 1,3 
  DO 20 K = 1,3 
     20 AINV(J,K)=AINV(J,K)/DETA 
   10  CONTINUE 
END SUBROUTINE INVERT 



!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC!


subroutine sanjiao(N,A_MESH,G,X)
INTEGER N
REAL L,U,R,Y,A_MESH,X,G	
DIMENSION A_MESH(N,N),X(N),U(N),L(N),R(N),Y(N),G(N)
U(1)=A_MESH(1,1)
R(1)=A_MESH(1,2)
Y(1)=G(1)

DO I=2,N-1
L(I)=A_MESH(I,I-1)/U(I-1)
U(I)=A_MESH(I,I)-L(I)*R(I-1)
R(I)=A_MESH(I,I+1)
Y(I)=G(I)-L(I)*Y(I-1)
END DO

L(N)=A_MESH(N,N-1)/U(N-1)
U(N)=A_MESH(N,N)-L(N)*R(N-1)
Y(N)=G(N)-L(N)*Y(N-1)

X(N)=Y(N)/U(N)
DO I=N-1,1,-1
X(I)=1./U(I)*(Y(I)-R(I)*X(I+1))
END DO 

END SUBROUTINE sanjiao






subroutine matrix_multiply(A,B,C)
implicit none
integer,parameter:: n=3            

real(8)::a(n,n),b(n,n),c(n,n)
integer::i,j,k

do i=1,n
do j=1,n

  c(i,j)=0.
  do k=1,n
  c(i,j)=c(i,j)+a(i,k)*b(k,j)
  end do

end do
end do

end subroutine matrix_multiply




subroutine series_multiply(A,B,C)
implicit none
integer,parameter:: n=3            

real(8)::a(n,n),b(n),c(n)
integer::j,k


do j=1,n

  c(j)=0.
  do k=1,n
  c(j)=c(j)+a(j,k)*b(k)
  end do

end do


end subroutine series_multiply




subroutine agaus2(a,b,x)
 implicit none
 integer,parameter:: n=3            !!!!!!!!!!!计算8阶的方程组
 integer::i,j,l,k,is
 real(8)::a(n,n),x(n),b(n),js(n)

 real(8)::d,t

   l=1      !逻辑变量,判断是否有解
   do k=1,n-1
     d=0.0
     do i=k,n
      do j=k,n
       if (abs(a(i,j))>d) then
         d=abs(a(i,j))
         js(k)=j
         is=i
       end if
   end do
     end do      !把行绝对值最大的元素换到主元位置
     if (d+1.0==1.0) then   
       l=0
     else     !最大元素为0无解
   if(js(k)/=k) then

    do i=1,n
      t=a(i,k)
   a(i,k)=a(i,js(k))
   a(i,js(k))=t
       end do       !最大元素不在K行，K行
   end if
   if(is/=k) then
    do j=k,n
     t=a(k,j)
  a(k,j)=a(is,j)
  a(is,j)=t     !交换到K列
       end do
    t=b(k)
    b(k)=b(is)
    b(is)=t    
   end if          !最大元素在主对角线上
     end if     !消去
  if (l==0) then
    write(*,100)
    write(*,*) "l==0"
    pause
    return
  end if
    do j=k+1,n
     a(k,j)=a(k,j)/a(k,k)
       end do
       b(k)=b(k)/a(k,k)     !求三角矩阵
    do i=k+1,n
     do j=k+1,n
    a(i,j)=a(i,j)-a(i,k)*a(k,j)
        end do
  b(i)=b(i)-a(i,k)*b(k)
       end do
     end do
     if (abs(a(n,n))+1.0==1.0) then
    l=0
    write(*,100)
    write(*,*) "l==0"
    pause
    return
  end if
  x(n)=b(n)/a(n,n)
  do i=n-1,1,-1
    t=0.0
    do j=i+1,n
      t=t+a(i,j)*x(j)
       end do
       x(i)=b(i)-t
     end do
100 format(1x,'fail')
    js(n)=n
 do k=n,1,-1
   if (js(k)/=k) then
    t=x(k)
    x(k)=x(js(k))
    x(js(k))=t
   end if
    end do
    return 
 end 