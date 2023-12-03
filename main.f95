program main
use modd

implicit none

integer :: N, M
real(8) :: L_x, L_y
real(8) :: T_x, dT_x, T_0
real(8) :: D
real(8) :: C
real(8) :: dx, dy, dt, t, t_stop
real(8), allocatable :: x(:), y(:)
real(8), allocatable :: u_old(:,:), u_new(:,:)

! Read the following parameters from file 'INPUT':
call InitializeParameters(L_x, L_y, N, M, T_x, dT_x, T_0, D, t_stop)

! Calculate grid step sizes dx, dy and dt, as well as the arrays of coordinates "x" and "y":
allocate(x(0:N-1), y(0:M-1))
call InitializeGrid(N, M, D, L_x, L_y, x, dx, y, dy, dt)

! Allocate temperature map:
allocate(u_old(0:N-1,0:M-1), u_new(0:N-1,0:M-1))

! Set initial condition (t = 0):
call SetIC(N, M, T_0, u_old)

! Start timer:
t = 0.d0

do while (t <= t_stop)
	! Set boundary condition:
	call SetBC(N, M, T_x, dT_x, dx, u_old)
	! Perform one time step:
    call Step(N, M, D, dx, dy, dt, T_x, dT_x, u_old, u_new)
    ! Update initial conditions:
    call UpdateIC(N, M, u_old, u_new)
    ! Update timer:
	t = t + dt
end do

! Save data to file 'RESULT':
call SaveData(N, M, dt, t_stop, u_new)

deallocate(x, y)

end
