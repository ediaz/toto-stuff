module stopwatch
implicit none
   
logical, private :: running 
real    :: time_ , start_
!integer, private, allocatable, dimension(:) :: time_array

contains



subroutine stopwatch_init()

running = .FALSE.
time_ = 0.0 
start_ = 0.0 
!allocate(time_array(8))

end subroutine

subroutine stopwatch_close()
  !  deallocate(time_array)
end subroutine



function get_time_sw() result (t) 
    real :: t


    ! The way I measure time is not error-save.
    ! This won't work if I use the sw for long runs
    ! which include a change in the day, month, or year. 
    
    ! Beware of using this sofware at midnight, or even
    ! worse, at new year's evening.
    
    call cpu_time(t)
!    t  = (time_array(5)*60*60 + time_array(6)*60 &
!             + time_array(7) +time_array(8)*1d-03) 
end function

subroutine start_sw()
    if (.not.running) then
        start_ = get_time_sw() 
        running = .TRUE.
    endif
end subroutine 

subroutine stop_sw()
    if (running) then 
        time_ = time_ +get_time_sw()  -start_ 
        running=.FALSE.
    end if
end subroutine stop_sw

subroutine reset_sw()
    call stop_sw() 
    time_ = 0.0
end subroutine reset_sw

subroutine restart_sw()
    call stop_sw()
    call start_sw()
end subroutine restart_sw

function time_sw() result (t )
real :: t
    if (running ) then
        t = time_ -start_ + get_time_sw() 
    else
        t = time_
    end if
end function time_sw
        
end module stopwatch
