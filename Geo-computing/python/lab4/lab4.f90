program lab4
use dsp
use stopwatch
use random
implicit none


real, allocatable, dimension(:,:) :: x,y
integer:: i,n1,n2



n1 =10000; n2=1000


allocate(x(n1,n2),y(n1,n2))
call get_random(x,398793)






do i=1,100
call dsp_smooth1p_2d(0.9,x,y,n1,n2)
print *, mean2(y,n1*n2)
enddo


deallocate(x,y)
end program lab4
