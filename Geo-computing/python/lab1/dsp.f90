module dsp
implicit none

contains


subroutine smooth_dsp(a,x,y)
real  :: a
real :: b,yi
real ,dimension(:)  :: x,y
integer  :: i , n 


n=SIZE(x) 
b=1.0-a

yi =x(1)
y(1) = yi

do i=2,n-1
    yi = a*yi +b*x(i)
    y(i) = yi
enddo

yi = (a*yi +x(n))/(1.0+a) 
y(n) = yi

do i=n-1,1,-1
    yi = a*yi +b*y(i)
    y(i) = yi
enddo 

end subroutine 

  subroutine smooth_dsp2(a,x,y)
  real :: a,b,yi,sx,sy
  real ,dimension(:)  :: x,y
  integer  :: i , n 
  
  n=SIZE(x) 
  b=1.0-a
  sx = b; sy = a
  yi = sx*x(1)
  y(1) = yi
  
  do i=2,n-1
    yi = a*yi +b*x(i)
    y(i) = yi
  enddo
  sx = sx/(1.0+a) ; sy = sy/(1.0+a)
  yi = sy*yi +sx*x(n) 
  y(n) = yi
  
  do i=n-1,1,-1
    yi = a*yi +b*y(i)
    y(i) = yi
  enddo 
  end subroutine smooth_dsp2 



function mean_dsp(x) result(m)
    real :: m
    real, dimension(:),intent(in) :: x
    integer :: n 
    
    n= SIZE(x) 
    m=0.0 

    m=sum(x)/n
end function





end module dsp
