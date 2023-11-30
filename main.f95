program main
use modd

implicit none

integer :: N, M
real(8) :: L_x, L_y
real(8) :: sigma
real(8) :: T_x, dT_y, T_0
real(8) :: D
real(8) :: dx, dy, dt
real(8), allocatable :: x(:), y(:)

! Read the following parameters from file 'INPUT':
call InitializeParameters(L_x, L_y, N, M, sigma, T_x, dT_y, T_0, D)

! Get grid step sizes dx and dt, as well as the array of coordinates "x":
allocate(x(0:N-1), y(0:M-1))
call InitializeGrid(N, M, D, L_x, L_y, x, dx, y, dy, dt)

! Set initial condition (t = 0):
!call SetIC(N, M, T_0, T_old)

deallocate(x, y)

end
