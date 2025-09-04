module LightConvex
   use stdlib_linalg_constants, only: ilp, dp, lk
   use assert_m, only: assert
   implicit none(external)
   private

   public :: simplex
   public :: unbounded_status, optimal_status, infeasible_status, &
             maxiter_exceeded

   !-----------------------------
   !-----     Constants     -----
   !-----------------------------

   real(dp), parameter :: eps = epsilon(1.0_dp)
    !! Machine precision
   real(dp), parameter :: tol = sqrt(eps)
   integer(ilp), parameter :: unbounded_status = 1
    !! Return flag for an unbounded problem.
   integer(ilp), parameter :: optimal_status = 0
    !! Return flag for when an optimal solution has been computed.
   integer(ilp), parameter :: infeasible_status = -1
    !! Return flag for an infeasible problem.
   integer(ilp), parameter :: maxiter_exceeded = -2
    !! Return flag for excessive number of iterations.

   !======================================
   !======================================
   !=====                            =====
   !=====     LINEAR PROGRAMMING     =====
   !=====                            =====
   !======================================
   !======================================

   ! ----- Abstract types -----

   !> Base type for the pivoting rule used in the simplex algorith.
   type, abstract :: abstract_pivot_rule
   end type abstract_pivot_rule

   !> Base type to determine which method is used to find an
   !  initial feasible point.
   type, abstract :: abstract_feasible_initialization
   end type abstract_feasible_initialization

   !--------------------------------------------
   !-----     PRIMAL SIMPLEX ALGORITHM     -----
   !--------------------------------------------

   !----- Low-level algorithm -----

   interface simplex
      module subroutine dense_standard_simplex(A, nleq, ngeq, neq, iposv, &
                                               maxiter, info, pivot, &
                                               initialization)
         implicit none(external)
         real(dp), intent(inout) :: A(:, :)
        !! Simplex tableau of dimension n+2 x m
         integer(ilp), intent(in) :: nleq, ngeq, neq
        !! Number of constraints of each type.
         integer(ilp), intent(in) :: maxiter
        !! Maximum number of iterations.
         integer(ilp), intent(out) :: info
        !! Return flag:
        !!  - info = -2 : Maximum number of iterations has been exceeded.
        !!  - info = -1 : Problem is infeasible.
        !!  - info = 0  : Optimal solution has been found.
        !!  - info = 1  : Objective function is unbounded.
         integer(ilp), intent(out) :: iposv(:)
        !! Book-keeping for the primal and slack variables.
         class(abstract_pivot_rule), intent(in) :: pivot
        !! Which pivoting rule is being used.
         class(abstract_feasible_initialization), intent(in) :: initialization
      end subroutine dense_standard_simplex
   end interface simplex

   !----- Pivot selection rules -----

   !> Default pivot selection rule.
   type, extends(abstract_pivot_rule), public :: Dantzig
   end type Dantzig

   !----- Initialization methods -----

   !> Auxiliary function maximization (default).
   type, extends(abstract_feasible_initialization), public :: auxiliary_function
   end type auxiliary_function

end module LightConvex
