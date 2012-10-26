program lab4
use dsp
use random
implicit none

real, allocatable, dimension(:,:) :: x,y
integer :: n1,n2,nsmooth,i
real    :: maxtime, a ,flops1,flops2
real, external :: benchmark2
integer ::  omp_get_thread_num



n1 =1001; n2=1003
allocate(x(n1,n2),y(n1,n2))


a = 0.99  ! smoothing parameter.


call get_random(x,398793)


maxtime=4.0

flops1 = benchmark(dsp_smooth1_2d,'smooth1   ',a,x,y,n1,n2,maxtime)
flops2 = benchmark(dsp_smooth1p_2d,'smooth1P  ',a,x,y,n1,n2,maxtime)
write(*,*) 'rate= ',flops2/flops1


flops1 = benchmark(dsp_smooth2,'smooth2   ',a,x,y,n1,n2,maxtime)
flops2 = benchmark(dsp_smooth2p,'smooth2p   ',a,x,y,n1,n2,maxtime)
write(*,*) 'rate= ',flops2/flops1



deallocate(x,y)
end program lab4
