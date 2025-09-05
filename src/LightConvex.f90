module LightConvex
   use stdlib_linalg_constants, only: ilp, dp, lk
   use assert_m, only: assert
   implicit none(external)
   private

   public :: simplex
   public :: unbounded_status, optimal_status, infeasible_status, maxiter_exceeded
   public :: is_optimal, is_feasible, is_unbounded

   !-----------------------------
   !-----     Constants     -----
   !-----------------------------

   real(dp), parameter :: eps = epsilon(1.0_dp)
    !! Machine precision
   real(dp), parameter, public :: tol = sqrt(eps)
    !! Tolerance used inside the different solvers.
   integer(ilp), parameter :: unbounded_status = 1
    !! Return flag for an unbounded problem.
   integer(ilp), parameter :: optimal_status = 0
    !! Return flag for when an optimal solution has been computed.
   integer(ilp), parameter :: infeasible_status = -1
    !! Return flag for an infeasible problem.
   integer(ilp), parameter :: maxiter_exceeded = -2
    !! Return flag for excessive number of iterations.

   !> Base type for defining convex problems.
   type, abstract :: abstract_cvx_problem
      private
      integer(ilp) :: status
   end type abstract_cvx_problem

   interface
      pure logical(lk) module function is_optimal(problem) result(bool)
         implicit none(external)
         class(abstract_cvx_problem), intent(in) :: problem
      end function is_optimal

      pure logical(lk) module function is_feasible(problem) result(bool)
         implicit none(external)
         class(abstract_cvx_problem), intent(in) :: problem
      end function is_feasible

      pure logical(lk) module function is_unbounded(problem) result(bool)
         implicit none(external)
         class(abstract_cvx_problem), intent(in) :: problem
      end function is_unbounded
   end interface

   !> Base type for defining convex solvers.
   type, abstract :: abstract_cvx_solver
   end type abstract_cvx_solver

   !> Base type for defining solutions to convex problems.
   type, abstract :: abstract_cvx_solution
   end type abstract_cvx_solution

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

   !----- Linear programming problem -----

   !> Definition of a linear program.
   type, extends(abstract_cvx_problem), public :: dense_lp_type
      real(dp), allocatable :: A(:, :)
        !! Simplex tableau.
      integer(ilp) :: nleq, ngeq, neq
        !! Number of constraints of each type.
   end type dense_lp_type

   interface linear_program
      type(dense_lp_type) module function assemble_dense_lp(c, Aleq, bleq, Ageq, bgeq, Aeq, beq) result(problem)
         implicit none(external)
         real(dp), intent(in) :: c(:)
         real(dp), intent(in), optional :: Aleq(:, :), bleq(:)
         real(dp), intent(in), optional :: Ageq(:, :), bgeq(:)
         real(dp), intent(in), optional :: Aeq(:, :), beq(:)
      end function assemble_dense_lp
   end interface
   public :: linear_program

   !> Solution of a linear program.
   type, extends(abstract_cvx_solution), public :: lp_solution
      real(dp), allocatable :: x(:)
        !! Solution of the primal problem.
      real(dp), allocatable :: y(:)
        !! Solution of the dual problem.
      real(dp), allocatable :: s(:)
        !! Slack variables for the primal problem.
      real(dp), allocatable :: t(:)
        !! Slack variables for the dual problem.
      real(dp) :: objective_value
        !! Objective value at the optimum.
   end type lp_solution

   !--------------------------------------------
   !-----     PRIMAL SIMPLEX ALGORITHM     -----
   !--------------------------------------------

   type, extends(abstract_cvx_solver), public :: PrimalSimplex
      class(abstract_pivot_rule), allocatable :: pivot
      class(abstract_feasible_initialization), allocatable :: initialization
      integer(ilp) :: maxiter
   end type PrimalSimplex

   interface PrimalSimplex
      type(PrimalSimplex) module function initialize_primal_simplex_alg(pivot, initialization, maxiter) result(alg)
         implicit none(external)
         class(abstract_pivot_rule), intent(in), optional :: pivot
         class(abstract_feasible_initialization), intent(in), optional :: initialization
         integer(ilp), intent(in), optional :: maxiter
      end function initialize_primal_simplex_alg
   end interface PrimalSimplex

   !----- High-level interface -----

   interface solve
      type(lp_solution) module function solve_with_dense_simplex(problem, alg) result(solution)
         implicit none(external)
         type(dense_lp_type), intent(inout) :: problem
         type(PrimalSimplex), intent(in) :: alg
      end function solve_with_dense_simplex
   end interface
   public :: solve

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

contains

   module procedure is_optimal
   bool = problem%status == optimal_status
   end procedure is_optimal

   module procedure is_feasible
   bool = problem%status /= infeasible_status
   end procedure is_feasible

   module procedure is_unbounded
   bool = problem%status == unbounded_status
   end procedure

end module LightConvex
