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

do i=2,n
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


function mean_dsp(x) result(m)
    real :: m
    real, dimension(:),intent(in) :: x
    integer :: i,n 
    
    n= SIZE(x) 
    m=0.0 

    m=sum(x)/n
    !do i=1,n
    !    m = m +x(i) 
    !enddo

    !m = m/n

end function





end module dsp
