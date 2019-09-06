program newtonmd

implicit none
real:: error
integer::istep,n,i,z,pointer,j, steps, q, points,outputfile,num, kappa
real :: mu1, mu2, stepsize, mu 
real:: func
real, dimension (:), allocatable :: w0, muarray
real, dimension (:,:), allocatable :: jarray,c,fx,x0, x1,xf, bin0
character (LEN = 60) :: inputval, filename1, filename2 
!CHARACTER*29 :: path ='C:\Users\17029\Desktop\gfortran\datfiles\'

error=1.0e-06
inputval = 'S'
kappa = 0

! initialising n 
print*,"Enter the dimension"
READ(*,*)  n

allocate ( x0(n,1) )
allocate ( bin0(n,1) )
allocate ( w0(n) )
allocate ( fx(n,1))
allocate ( x1(n,1))
allocate ( xf(n,1))
allocate ( jarray(n,n) )
allocate ( c(n,n))


print*,"Enter mu1 value"
READ(*,*)  mu1

print*,"Enter mu2 value"
READ(*,*)  mu2

print*,"Enter points value"
READ(*,*)  points

!print*,"Enter the Guess Vector"
!DO i = 1, n
!   READ(*,*)  x0(i,1)
!END DO

print*,"Enter the aeta Vector"
DO i = 1, n
   READ(*,*)  w0(i)
END DO

!print*,"Enter kappa value, Enter -1 for all"
!READ(*,*)  kappa


stepsize = (mu2-mu1)/(points+1)
allocate ( muarray(points + 2))


muarray(1) = mu1
do i = 2,(points+1)
muarray(i) = muarray(i-1) + stepsize
end do
muarray(points+2) = mu2
outputfile = 1  



do num = 0, (2**(n) -1)


write(filename1, '("Output",I0,".dat")') outputfile 
filename2=  filename1

print*, "this work 2"
call dectobin_gen(n,bin0,num)

print*, 'this works'
print*, bin0


!if ( sum(bin0) /=  kappa) CYCLE


x0 = bin0
print*, "this works 3"
open(90, file=filename2,action="write",POSITION="APPEND")

q = 1
do while ( q /= (points + 3))

mu = muarray(q)  


pointer= 1

xf = x0

!do while(   ((abs(NORM2(xf)-NORM2(x1))/abs(NORM2(xf))) > 0.00001) )


!do while(   (pointer /= 70) .or. ((abs(NORM2(xf)-NORM2(x0))/abs(NORM2(x0))) == 0 )   )


do while( pointer /= 70 )


z = 1 

print*, x0


!! generating jacobian 

call Jac(x0, mu, n, w0, jarray)

do i = 1,n
do j = 1,n
   print*, i,j,jarray(i,j)
end do
end do


!! generating it's inverse
print*, "generating inverse"


call inverse(jarray,c,n)

print*,'c is', c


!! generating fx(x0)

do while (z /= n+1  )
call eqn(x0, w0, z, n, mu, func)
fx(z,1) = func
z = z+1
end do

print*, fx

!! computing x1

x1 = x0 - matmul(c,fx)


!!do i = 1,n
!!do j = 1,n
!!   print*, i,j,c(i,j)
!!end do
!!end do

!!print*,'fx is', fx


!!print*,'c*fx is', matmul(c, fx) 

!!print*, size(c), size(fx)


!!print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!!print*, x1
!!print*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

x0 = x1
pointer = pointer + 1
end do

   write(90, '(F8.5,X)', advance='no') muarray(q)

   do j=1,n
      write(90, '(F8.5,X)', advance='no') x1(j,1)
   end do

   write(90, *) ''  ! this gives you the line break


q = q+1


end do

outputfile = outputfile + 1
close(90)
end do 

end program

!!!!!!!!!! SUBROUTINES !!!!!!!!!!!!!



subroutine eqn(x0, w0, i, n, mu, f)

integer, intent(in):: n
real :: mu 
real, DIMENSION(1:n) :: x0
real, DIMENSION(1:n) :: w0 
real, intent(out):: f
integer:: i
real :: sumn

call sum(x0, w0, i, n, sumn)

print*,'sum = ', sumn
print*, 'mu =', mu

f = ((x0(i))**2) - (x0(i)) - (mu*sumn)

return
END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sum(x0, w0, i, n, sumn)

integer, intent(in):: n
real, DIMENSION(n,1) :: x0 
real, DIMENSION(1:n) :: w0 
real, intent(out):: sumn 
integer:: i, p,k

k = i

p = k

sumn = 0

do while (k /= n)
sumn = sumn + (x0(p,1) - x0((k + 1),1))/(w0(p) - w0((k+1)) )
k = k + 1
end do

k = p

do while (k /= 1)
sumn = sumn + (x0(p,1) - x0((k - 1),1))/(w0(p) - w0((k-1)) )
k = k -1
end do

return 
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Jac(x0, mu, n, w0, jind)

integer :: n, j , i, k
real, DIMENSION(n,1) :: x0 
real, DIMENSION(1:n) :: w0 
real, DIMENSION(1:n, 1:n) :: jind
real :: mu, sumnaeta

print*, 'n value is ', n

j = 1

i = 1



do while( i /= n+1 )

do j = 1,n

print*, i, j

if (i == j) then
k = i
call sumaeta(w0, k, n, sumnaeta)
print*, 'once'
jind(i,j) =  2*(x0(i,1)) - 1 - (mu*sumnaeta) 
end if

if(i /= j) then
print*, 'twice'
jind(i,j) =  mu/(w0(i)-w0(j))
end if

end do 

i = i +1
end do


end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dectobin_gen(n,bin0,num)

integer :: n,num,i,temp
real, DIMENSION(n,1) :: bin0

do i = 1,n

bin0(i,1) = 0

end do


temp= num

do i=1,n

if (mod(temp,2)==1) bin0(i,1)=1


temp   =   temp/2

if (temp==0) exit


end do

end





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sumaeta( w0, i, n, sumnaeta)

integer :: n
real, DIMENSION(1:n) :: w0 
real :: sumnaeta, sumn 
integer:: i, p

p = i

sumn = 0

do while (i /= n)
sumn = sumn + (1/(w0(p)-w0(i+1)))
i = i + 1
end do

i = p

do while (i /= 1)
sumn = sumn + (1/(w0(p) - w0(i-1)))
i = i -1
end do

sumnaeta = sumn

return 
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
real :: a(n,n), c(n,n)
real :: L(n,n), U(n,n), b(n), d(n), x(n)
real :: coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

