program lab0
use stopwatch2
use dsp 


implicit none


type (sw_struct) :: sw

real,allocatable,dimension(:) :: x,y
real :: maxtime, a
integer  :: i,n,nsmooth 


n = 1001 ! number of samples

allocate ( x(n),y(n)) ! memory allocation

x =0.0 ! array initialization.
y = x 

! ---------------------
! input vector:
x(n)=1.0 
x(int(n/2)+1) = x(n)
x(1) = x(n)
! ---------------------

a = 0.99  ! smoothing parameter.

call stopwatch_init(sw)
maxtime=2.0 ! time for benchmark
nsmooth=0

call start_sw(sw)
do while (time_sw(sw) < maxtime ) 
    call smooth_dsp(a,x,y)
    nsmooth =nsmooth +1
enddo

call stop_sw(sw)



write (*,'(A)') 'Fortran'
write (*,'(A, I10)') 'nsmooth =', nsmooth
write (*,'(A, F11.8)') '   mean =', mean_dsp(y)
write (*,'(A, F12.9)') '   time =', time_sw(sw) 
write (*,'(A, I10)') ' mflops =', int(6.0d-6*n*nsmooth/time_sw(sw))



! close everything a deallocate: 
deallocate (y,x) 
call stopwatch_close()

end program lab0
