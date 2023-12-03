module modd
use matrix_solver
implicit none
real(8) :: pi = 4 * atan(1.d0)

contains

subroutine InitializeParameters(L_x, L_y, N, M, T_x, dT_x, T_0, D, t_stop)
! Subroutine for setting model parameters
integer :: N, M
real(8) :: L_x, L_y, T_x, dT_x, T_0, D, t_stop
open(unit = 1, file = 'INPUT')
read(1, *) L_x
read(1, *) L_y
read(1, *) N
read(1, *) M
read(1, *) T_x
read(1, *) dT_x
read(1, *) T_0
read(1, *) D
read(1, *) t_stop
end subroutine InitializeParameters

subroutine InitializeGrid(N, M, D, L_x, L_y, x, dx, y, dy, dt)
! Subroutine for initialising the grid
integer :: N, M, i
real(8) :: L_x, L_y, D
real(8) :: dx, dt, x(0:N-1), dy, y(0:M-1)
dx = L_x / (N - 1)
dy = L_y / (M - 1)
dt = 0.1d0
x(0) = 0; x(N-1) = L_x
do i = 1, N-2
	x(i) = x(i-1) + dx
enddo
y(0) = 0; y(M-1) = L_y
do i = 1, M-2
	y(i) = y(i-1) + dy
enddo
end subroutine InitializeGrid


subroutine SetIC(N, M, T_0, u)
! Subroutine for setting the initial condition
integer :: N, M, i, j
real(8) :: u(0:N-1, 0:M-1), T_0
do i = 0, N-1
    do j = 0, M-1
    	u(i, j) = T_0
	enddo
enddo
end subroutine SetIC

subroutine SetBC(N, M, T_x, dT_x, dx, u)
!Subroutine for setting the boundary condition
integer :: N, M, i
real(8) :: T_x, dT_x
real(8) :: dx
real(8) :: u(0:N-1,0:M-1)
do i = 0, N-1
    u(i, 0) = T_x
    u(i, M-1) = T_x
enddo
do i = 0, M-1
    u(1, i) = u(0, i) + dT_x * dx
    u(N-2, i) = u(N-1, i) + dT_x * dx
enddo
end subroutine SetBC

subroutine Step(N, M, D, dx, dy, dt, T_x, dT_x, u_old_2d, u_new_2d)
! Performing one time step:
integer :: N, M
integer :: i, j
real(8) :: D, dx, dy, dt
real(8) :: T_x, dT_x
real(8) :: u_old_2d(0:N-1, 0:M-1), u_new_2d(0:N-1, 0:M-1)
real(8), allocatable :: u_new(:)
real(8), allocatable :: alpha_vector(:), beta_vector(:), gamma_vector(:), b_vector(:)
! First half-step (alongside X):
allocate(alpha_vector(0:N-2), beta_vector(0:N-1), gamma_vector(1:N-1), b_vector(0:N-1), u_new(0:N-1))
! For BC:
beta_vector(0) = - 1.d0 / dx
alpha_vector(0) = 1.d0 / dx
b_vector(0) = dT_x
gamma_vector(N-1) = 1.d0 / dx
beta_vector(N-1) = - 1.d0 / dx
b_vector(N-1) = dT_x
do i = 1, N-2
	beta_vector(i) = - (2 * D / dx**2 + 2 / dt) ! the main diagonal
	alpha_vector(i) = D / dx**2 ! the upper subdiagonal
	gamma_vector(i) = D / dx**2 ! the lower subdiagonal
enddo
do j = 1, M-2
	do i = 1, N-2
		b_vector(i) = - D / dy**2 * (u_old_2d(i, j+1) + u_old_2d(i, j-1)) + (2 * D / dy**2 - 2 / dt) * u_old_2d(i, j) ! right parts
	enddo
	call tridiagonal(N, alpha_vector, beta_vector, gamma_vector, b_vector, u_new)
	u_new_2d(:, j) = u_new
enddo
deallocate(alpha_vector, beta_vector, gamma_vector, b_vector, u_new)
! Second half-step (alongside Y):
allocate(alpha_vector(0:M-2), beta_vector(0:M-1), gamma_vector(1:M-1), b_vector(0:M-1), u_new(0:M-1))
! For BC:
beta_vector(0) = 1.d0
alpha_vector(0) = 0.d0
b_vector(0) = T_x
gamma_vector(M-1) = 0.d0
beta_vector(M-1) = 1.d0
b_vector(M-1) = T_x
do j = 1, M-2
	beta_vector(j) = - (2 * D / dy**2 + 2 / dt) ! the main diagonal
	alpha_vector(j) = D / dy**2 ! the upper subdiagonal
	gamma_vector(j) = D / dy**2 ! the lower subdiagonal
enddo
do i = 1, N-2
	do j = 1, M-2
		b_vector(j) = - D / dx**2 * (u_old_2d(i+1, j) + u_old_2d(i-1, j)) + (2 * D / dx**2 - 2 / dt) * u_old_2d(i, j) ! right parts
	enddo
	call tridiagonal(M, alpha_vector, beta_vector, gamma_vector, b_vector, u_new)
	u_new_2d(i, :) = u_new
enddo
deallocate(alpha_vector, beta_vector, gamma_vector, b_vector)
end subroutine Step

subroutine UpdateIC(N, M, u_old, u_new)
!Subroutine for updating the initial condition
integer :: N, M
real(8) :: u_old(0:N-1, 0:M-1), u_new(0:N-1, 0:M-1)
u_old = u_new
end subroutine UpdateIC

subroutine SaveData(N, M, dt, t_stop, u)
!Subroutine to save data
integer :: N, M, j
real(8) :: u(0:N-1, 0:M-1)
real(8) :: dt, t_stop
open(unit = 2, file = 'RESULT')
write(2, *) N, M, dt, t_stop
do j = 0, M-1
	write(2, *) u(:, j)
enddo
end subroutine SaveData


end module
