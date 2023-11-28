module modd
use matrix_solver
implicit none
real(8) :: pi = 4 * atan(1.d0)

contains

subroutine InitializeParameters(L_x, L_y, N, M, sigma, T_x, dT_y, T_0, D)
! Subroutine for setting model parameters
integer :: N, M
real(8) :: L_x, L_y, sigma, T_x, dT_y, T_0, D
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

end module
