program main
use modd

implicit none

integer :: N, M
real(8) :: L_x, L_y
real(8) :: sigma
real(8) :: T_x, dT_y, T_0
real(8) :: D

! Read the following parameters from file 'INPUT':
call InitializeParameters(L_x, L_y, N, M, sigma, T_x, dT_y, T_0, D)


end
