module match
	use variables 
	use hcore
	use physical_properties
    implicit none
    contains
	
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TDMA算法!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TDMA算法!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TDMA算法!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!形参说明：                                                     !!!!!!!!!!!!!!
	!!!!!!!!!!!!!!          A—————双精度二维数组,体积为N*N,输入参数.方程组的!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!          系数矩阵,要求主对角线占绝对优势.                     !!!!!!!!!!!!!!
	!!!!!!!!!!!!!!          B—————双精度实型一维数组,长度为N,输入参数.方程组!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!          右端的常数向量.                                      !!!!!!!!!!!!!!
	!!!!!!!!!!!!!!          N—————整型变量,输入参数.方程组阶数.             !!!!!!!!!!!!!!
	!!!!!!!!!!!!!!          X—————双精度实型一维数组,长度为N,输出参数.返回方!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!          程组的解向量.                                        !!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine TDMA(n,A,B,X)
	integer n,i
	real(8)::L,U,R,Y
	dimension A(n,n),X(n),U(n),L(n),R(n),Y(n),B(n)
	double precision A,B,X
	U(1)=A(1,1)
	R(1)=A(1,2)
	Y(1)=B(1)
	do i=2,n-1
		L(i)=A(i,i-1)/U(i-1)
		U(i)=A(i,i)-L(i)*R(i-1)
		R(i)=A(i,i+1)
		Y(i)=B(i)-L(i)*Y(i-1)
	end do
	L(n)=A(n,n-1)/U(n-1)
	U(n)=A(n,n)-L(n)*R(n-1)
	Y(n)=B(n)-L(n)*Y(n-1)
	X(n)=Y(n)/U(n)
	do i=n-1,1,-1
		X(i)=1./U(i)*(Y(i)-R(i)*X(i+1))
	end do 
	end subroutine
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TDMA算法!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TDMA算法!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TDMA算法!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
	
	
	
end module
	