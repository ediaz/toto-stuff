module stopwatch2

implicit none

type sw_struct
  logical, private :: running 
  real    :: time_ , start_
end type sw_struct








   








contains



subroutine stopwatch_init(sw)
type (sw_struct) :: sw

sw%running = .FALSE.
sw%time_ = 0.0 
sw%start_ = 0.0 
!allocate(time_array(8))

end subroutine

subroutine stopwatch_close()
  !  deallocate(time_array)
end subroutine



function get_time_sw() result (t) 
  real :: t


  call cpu_time(t)
end function

subroutine start_sw(sw)
  type (sw_struct) :: sw

    if (.not.sw%running) then
        sw%start_ = get_time_sw() 
        sw%running = .TRUE.
    endif
end subroutine 

subroutine stop_sw(sw)
  type (sw_struct) :: sw

    if (sw%running) then 
        sw%time_ = sw%time_ +get_time_sw()  -sw%start_ 
        sw%running=.FALSE.
    end if
end subroutine stop_sw

subroutine reset_sw(sw)
  type (sw_struct) :: sw

    call stop_sw(sw) 
    sw%time_ = 0.0
end subroutine reset_sw

subroutine restart_sw(sw)
  type (sw_struct) :: sw
  
    call stop_sw(sw)
    call start_sw(sw)
end subroutine restart_sw

function time_sw(sw) result (t )
  type (sw_struct) :: sw
  real :: t

    if (sw%running ) then
        t = sw%time_ -sw%start_ + get_time_sw() 
    else
        t = sw%time_
    end if
end function time_sw
        
end module stopwatch2
