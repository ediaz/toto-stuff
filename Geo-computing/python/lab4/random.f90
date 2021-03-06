module random 
implicit none


interface get_random
  module procedure get_random_1d, get_random_2d    
end interface



contains

subroutine init_seed(seed)
  integer seed
call srand(seed)
end subroutine


subroutine get_random_1d(vector,iseed)

  real, dimension(:),intent(inout) :: vector
  integer            :: iseed
  integer i

  call init_seed(iseed)

  do i=1,size(vector)
    vector(i) = rand()
  enddo


end subroutine get_random_1d

subroutine get_random_2d(vector,iseed)

  real, dimension(:,:),intent(inout) :: vector
  integer            :: iseed
  integer  :: i2,i1

  call init_seed(iseed)

  do i2=1,size(vector,2)
    do i1=1,size(vector,1)
      vector(i1,i2) = rand()
    enddo
  enddo


end subroutine get_random_2d







end module random
