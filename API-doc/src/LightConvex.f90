module LightConvex
   use stdlib_linalg_constants, only: ilp, dp, lk
   use assert_m, only: assert
   implicit none(external)
   private

   public :: simplex

   !-----------------------------
   !-----     Constants     -----
   !-----------------------------

   real(dp), parameter :: eps = epsilon(1.0_dp)

   !--------------------------------------
   !-----     LINEAR PROGRAMMING     -----
   !--------------------------------------

   interface simplex
      pure module subroutine dense_standard_simplex(A, nleq, ngeq, neq, iposv, info)
         implicit none(external)
         real(dp), intent(inout) :: A(:, :)
        !! Simplex tableau of dimension n+2 x m
         integer(ilp), intent(in) :: nleq, ngeq, neq
        !! Number of constraints of each type.
         integer(ilp), intent(out) :: info
        !! Return flag:
        !!  - info = 1  : Optimal solution has been found.
        !!  - info = 0  : Objective function is unbounded.
        !!  - info = -1 : Problem is infeasible.
         integer(ilp), intent(out) :: iposv(:)
        !! Book-keeping for the primal and slack variables.
      end subroutine dense_standard_simplex
   end interface simplex

end module LightConvex
