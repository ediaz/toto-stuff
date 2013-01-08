module dsp
!$ use omp_lib
implicit none

contains



!
!
!  Recursive exponential dsp_smoothing.
!  This modules contains functions for 2d and 3d
!  processing.
!
!  It also contains interfaces so the user calls
!  one function name, no matter if the input is 
!  1D or 2D
!
!
!  Author : Esteban D\'{i}az 
!  date   : 10/2012
!
!

subroutine dsp_smooth1d(a,x,y,n1)
real :: a
real ,dimension(:)  :: x,y
integer :: n1

call dsp_smooth1_1d(a,x,y,n1)

end subroutine


subroutine dsp_smooth2d(a,x,y,n1,n2)
real :: a
real ,dimension(:,:)  :: x,y
integer :: n1,n2

call dsp_smooth1_2d(a,x,y,n1,n2)
call dsp_smooth2(a,y,y,n1,n2)

end subroutine



!
! Recursive exponential dsp_smoothing
! for a 1D vector
!
! Author: Esteban D\'{i}az 
! Date  : 15.10.2012

subroutine dsp_smooth1_1d(a,x,y,n1)
real  :: a
real :: b,yi
real ,dimension(:)  :: x,y
integer  :: i 
integer n1,n2


b=1.0-a

yi =x(1)
y(1) = yi

do i=2,n1-1
    yi = a*yi +b*x(i)
    y(i) = yi
enddo

yi = (a*yi +x(n1))/(1.0+a) 
y(n1) = yi

do i=n1-1,1,-1
    yi = a*yi +b*y(i)
    y(i) = yi
enddo 

end subroutine 


!
! Recursive exponential dsp_smoothing
! along first dimension for a 2D 
! array
!
! Author: Esteban D\'{i}az 
! Date  : 15.10.2012
subroutine dsp_smooth1_2d(a,x,y,n1,n2)
real :: a
real, dimension(:,:)  :: x
real, dimension(:,:)  :: y
integer:: n1,n2
integer :: i



y=0.0
do i=1,n2
  call dsp_smooth1_1d(a,x(:,i),y(:,i),n1)
enddo
end subroutine 












!
! Recursive exponential dsp_smoothing
! along first dimension for a 2D 
! array
!
! Author: Esteban D\'{i}az 
! Date  : 15.10.2012
subroutine dsp_smooth1p_2d(a,x,y,n1,n2)
real :: a
real, dimension(:,:)  :: x
real, dimension(:,:)  :: y
real, dimension(:), allocatable :: xaux,yaux
integer:: n1,n2
integer :: i

y=0.0
!$omp parallel do private (i) shared(a,x,y,n1) schedule(dynamic)
do i=1,n2
  call dsp_smooth1_1d(a, x(:,i) ,y(:,i),n1)
enddo
!$omp end parallel do
end subroutine 


!
! Recursive exponential dsp_smoothing
! along second dimension for a 2D 
! array, by transposing the array
! (I did this just for QC of Smooth2)
!
! It looks very clean in Fortran
!
! Author: Esteban D\'{i}az 
! Date  : 15.10.2012

subroutine dsp_smooth2t(a,x,y,n1,n2)
real :: a
real, dimension(:,:)  :: x
real, dimension(:,:)  :: y
integer:: n1,n2
integer :: i
y=0.0
do i=1,n1
  call dsp_smooth1_1d(a,x(i,:),y(i,:),n2)
enddo
end subroutine 

!
! Recursive exponential dsp_smoothing
! along second dimension for a 2D 
! array
!
! Author: Esteban D\'{i}az 
! Date  : 15.10.2012
subroutine dsp_smooth2 (a,x,y,n1,n2)
real    :: a,b
integer :: i1,i2
integer n1,n2
real, dimension(:,:) :: x,y 

b = 1.0 -a

do i1=1,n1
  y(i1,1) = x(i1,1) !first
enddo

do i2=2,n2-1
  do i1=1,n1
    y(i1,i2) = a*y(i1,i2-1)+b*x(i1,i2) !forward
  enddo
enddo

do i1=1,n1
  y(i1,n2) = (a*y(i1,n2-1) +x(i1,n2))/(1.0+a) !last
enddo


do i2=n2-1,1,-1
  do i1=1,n1
    y(i1,i2) = a*y(i1,i2+1)+b*y(i1,i2) !forward
  enddo
enddo


end subroutine


subroutine dsp_smooth2p(a,x,y,n1,n2)
real :: a
real, dimension(:,:)  :: x
real, dimension(:,:)  :: y
integer:: n1,n2
integer :: i
integer :: chunk,nthreads 
INTEGER, EXTERNAL :: omp_get_max_threads



nthreads =omp_get_max_threads()

chunk = int(n1/nthreads)

!$omp parallel do private (i) shared(a,x,y,chunk) schedule(dynamic)
do i = 1,n1,chunk
  call dsp_smooth2_chunk(a,x,y,i,chunk)
enddo
!$omp end parallel do

endsubroutine


! Subroutine for chunky processing
subroutine dsp_smooth2_chunk(a,x,y,lb,chunk)
real :: a,b
real, dimension(:,:)  :: x
real, dimension(:,:)  :: y
integer :: n1,n2, chunk,i1,i2
integer :: lb,ub 

n2 = size(x,2)
n1 = size(x,1)

ub = lb +chunk - 1 
ub = min(ub,n1)
b = 1.0 -a

do i1=lb,ub
  y(i1,1) = x(i1,1) !first
enddo

do i2=2,n2-1
  do i1=lb,ub
    y(i1,i2) = a*y(i1,i2-1)+b*x(i1,i2) !forward
  enddo
enddo

do i1=lb,ub
  y(i1,n2) = (a*y(i1,n2-1) +x(i1,n2))/(1.0+a) !last
enddo

do i2=n2-1,1,-1
  do i1=lb,ub
    y(i1,i2) = a*y(i1,i2+1)+b*y(i1,i2) !forward
  enddo
enddo

end subroutine



function mean1(x,n) result(m)
    real :: m
    real, dimension(:),intent(in) :: x
    integer :: n 
    
    m=0.0 

    m=sum(x)/n
end function


function mean2(x,n) result(m)
    real :: m
    real, dimension(:,:),intent(in) :: x
    integer :: n 
    
    m=0.0 

    m=sum(x)/n
end function


function benchmark(f,fname,a,x,y,n1,n2,maxtime) result(flops)
  use stopwatch
  type (sw_struct) :: sw
  real :: maxtime, flops
  integer :: nsmooth 
  character(10) :: fname


  real :: a
  real , dimension(:,:) :: x,y
  integer :: n1,n2

  interface 
    subroutine f(a,x,y,n1,n2)
      real :: a
      real , dimension(:,:) :: x,y
      integer :: n1,n2
      
    end subroutine
  end interface


  call stopwatch_init(sw)
  nsmooth=0
  
  call start_sw(sw)
  do while (time_sw(sw) < maxtime ) 
      call f(a,x,y,n1,n2)
      nsmooth =nsmooth +1
  enddo
  call stop_sw(sw)
  
  write(*,'(A,I5,A,F13.10)')fname//': mflops=', &
      int(6.0d-6*n1*n2*nsmooth/time_sw(sw)), &
      '  mean= ',mean2(y,n1*n2)

  flops = 6.0d-6*n1*n2*nsmooth/time_sw(sw)

endfunction

end module dsp
