module modd
use matrix_solver
implicit none
real(8) :: pi = 4 * atan(1.d0)

contains

subroutine InitializeParameters(L_x, L_y, N, M, sigma, T_x, dT_y, T_0, D)
! Subroutine for setting model parameters
integer :: N, M
real(8) :: L_x, L_y, sigma, T_x, dT_y, T_0, D, C
open(unit = 1, file = 'INPUT')
read(1, *) L_x
read(1, *) L_y
read(1, *) N
read(1, *) M
read(1, *) sigma
read(1, *) T_x
read(1, *) dT_y
read(1, *) T_0
read(1, *) D
end subroutine InitializeParameters

subroutine InitializeGrid(N, M, D, L_x, L_y, x, dx, y, dy, dt)   !dt = ???
! Subroutine for initialising the grid
real(8) :: L_x, L_y, D
integer :: N, M, i
real(8) :: dx, dt, x(0:N-1), dy, y(0:M-1)
dx = L_x / (N - 1)
dy = L_y / (M - 1)
dt = 
x(0) = 0; x(N-1) = L_x
do i = 1, N-2
	x(i) = x(i-1) + dx
enddo
y(0) = 0; y(N-1) = L_y
do i = 1, M-2
	y(i) = y(i-1) + dy
enddo
end subroutine InitializeGrid


subroutine SetIC(N, M, T_0, T_old)
! Subroutine for setting the initial condition
integer :: N, M, i, j
real(8) :: T_old(0:N-1,0:M-1), T_0
do i = 0, N-1
    do j = 0, M-1
    	T_old(i,j) = T_0
	enddo
enddo
end subroutine SetIC

end module
