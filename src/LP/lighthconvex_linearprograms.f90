module lightconvex_lp
   use assert_m, only: assert => assert_always
   use stdlib_optval, only: optval
   use lightconvex_constants, only: ilp, dp, lk, eps, tol, &
                                    optimal_status, infeasible_status, &
                                    unbounded_status, maxiter_exceeded
   use lightconvex_abstract, only: abstract_cvx_problem, &
                                   abstract_cvx_solver, &
                                   abstract_cvx_solution, &
                                   is_optimal, is_feasible, is_unbounded
   use lightconvex_linalg, only: qr_type, QR
   implicit none(external)
   private

   !==================================================
   !=====                                        =====
   !=====     DEFINITION OF A LINEAR PROGRAM     =====
   !=====                                        =====
   !==================================================

   !----- Linear programming problem -----

   !> Definition of a linear program.
   type, extends(abstract_cvx_problem), public :: dense_lp_type
      real(dp), allocatable :: c(:)
        !! Linear cost function
      real(dp), allocatable :: Aleq(:, :), bleq(:)
        !! <= inequality constraints.
      real(dp), allocatable :: Ageq(:, :), bgeq(:)
        !! >= inequality constraints.
      real(dp), allocatable :: Aeq(:, :), beq(:)
        !! == equality constraints.
   end type dense_lp_type

   interface linear_program
      type(dense_lp_type) module function create_dense_linear_program(c, Aleq, bleq, Ageq, bgeq, Aeq, beq) result(problem)
         implicit none(external)
         real(dp), intent(in) :: c(:)
            !! Linear cost function.
         real(dp), intent(in), optional :: Aleq(:, :), bleq(:)
            !! <= inequality constraints.
         real(dp), intent(in), optional :: Ageq(:, :), bgeq(:)
            !! >= inequality constraints.
         real(dp), intent(in), optional :: Aeq(:, :), beq(:)
            !! == equality constraints.
      end function create_dense_linear_program
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

   !==============================================
   !=====                                    =====
   !=====     LINEAR PROGRAMMING SOLVERS     =====
   !=====                                    =====
   !===============================================

   ! ----- Abstract types -----

   !> Base type for the pivoting rule used in the simplex algorith.
   type, abstract :: abstract_pivot_rule
   end type abstract_pivot_rule

   !> Base type to determine which method is used to find an
   !  initial feasible point.
   type, abstract :: abstract_feasible_initialization
   end type abstract_feasible_initialization

   !> Auxiliary function maximization (default).
   type, extends(abstract_feasible_initialization), public :: auxiliary_function
   end type auxiliary_function

   !--------------------------------------------
   !-----     PRIMAL SIMPLEX ALGORITHM     -----
   !--------------------------------------------

   type, extends(abstract_cvx_solver), public :: PrimalSimplex
      class(abstract_pivot_rule), allocatable :: pivot
        !! Rule to choose the pivoting column.
      class(abstract_feasible_initialization), allocatable :: initialization
        !! Method to find an initial feasible point.
      integer(ilp) :: maxiter
        !! Maximum number of iterations in the simplex method.
   end type PrimalSimplex

   interface PrimalSimplex
      type(PrimalSimplex) module function initialize_primal_simplex_alg(pivot, initialization, maxiter) result(alg)
         implicit none(external)
         class(abstract_pivot_rule), intent(in), optional :: pivot
            !! Rule to choose the pivoting column.
         class(abstract_feasible_initialization), intent(in), optional :: initialization
            !! Method to find an initial feasible point.
         integer(ilp), intent(in), optional :: maxiter
            !! Maximum number of iterations in the simplex method.
      end function initialize_primal_simplex_alg
   end interface PrimalSimplex

   !----- Low-level algorithm -----

   interface simplex
      pure module subroutine dense_standard_simplex(A, nleq, ngeq, neq, iposv, &
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
   public :: simplex

   !----- Pivot selection rules -----

   !> Default pivot selection rule.
   type, extends(abstract_pivot_rule), public :: Dantzig
   end type Dantzig

   !-----------------------------------------
   !-----     PRIMAL AFFINE SCALING     -----
   !-----------------------------------------

   type, extends(abstract_cvx_solver), public :: PrimalAffineScaling
      class(abstract_feasible_initialization), allocatable :: initialization
        !! Method to find an initial feasible point.
      integer(ilp) :: maxiter
        !! Maximum number of iterations in the PAS method.
   end type PrimalAffineScaling

   interface PrimalAffineScaling
      type(PrimalAffineScaling) module function initialize_primal_affine_scaling_alg(initialization, maxiter) result(alg)
         implicit none(external)
         class(abstract_feasible_initialization), intent(in), optional :: initialization
            !! Method to find an initial feasible point.
         integer(ilp), intent(in), optional :: maxiter
            !! Maximum number of iterations.
      end function initialize_primal_affine_scaling_alg
   end interface PrimalAffineScaling

   !----- High-level interface -----

   interface solve
      type(lp_solution) module function solve_with_dense_simplex(problem, alg) result(solution)
         implicit none(external)
         type(dense_lp_type), intent(inout) :: problem
            !! Problem to solve.
         type(PrimalSimplex), intent(in) :: alg
            !! Algorithm used.
      end function solve_with_dense_simplex

      type(lp_solution) module function solve_with_dense_PAS(problem, x0, alg) result(solution)
         implicit none(external)
         type(dense_lp_type), intent(inout) :: problem
            !! Problem to solve.
         real(dp), optional, intent(in) :: x0
            !! Initial guess for warm starting.
         type(PrimalAffineScaling), intent(in) :: alg
            !! Algorithm used.
      end function solve_with_dense_PAS
   end interface
   public :: solve

contains

   !==================================================
   !=====                                        =====
   !=====     DEFINITION OF A LINEAR PROGRAM     =====
   !=====                                        =====
   !==================================================

   !-----------------------------------------
   !-----     DENSE LINEAR PROGRAMS     -----
   !-----------------------------------------

   module procedure create_dense_linear_program
   integer(ilp) :: n
    !! Number of variables.

   n = size(c); problem%c = c

   !> Sanity check.
   call assert(assertion=(present(Aleq) .and. present(bleq)) .or. (.not. present(Aleq) .and. .not. present(bleq)), &
               description="Specification of <= constraints incomplete. Either Aleq or bleq is missing.")
   call assert(assertion=(present(Ageq) .and. present(bgeq)) .or. (.not. present(Ageq) .and. .not. present(bgeq)), &
               description="Specification of >= constraints incomplete. Either Aleq or bleq is missing.")
   call assert(assertion=(present(Aeq) .and. present(beq)) .or. (.not. present(Aeq) .and. .not. present(beq)), &
               description="Specification of == constraints incomplete. Either Aleq or bleq is missing.")
   call assert(assertion=present(Aleq) .or. present(Ageq) .or. present(Aeq), &
               description="No constraint has been provided. Ill-posed problem.")

   !> Consistency of the <= inequalities.
   if (present(Aleq)) then
      !> Check dimensions.
      call assert(assertion=size(Aleq, 1) == size(bleq), &
                  description="Aleq and bleq have an inconsistent number of rows.")
      call assert(assertion=size(Aleq, 2) == n, &
                  description="Number of columns of Aleq is inconsistent with the number of variables (size(c)).")
      call assert(assertion=all(bleq >= -eps), &
                  description="Right-hand side vector bleq needs to be non-negative.")

      !> If all good, allocate arrays.
      problem%Aleq = Aleq; problem%bleq = bleq
   end if

   !> Consistency of the >= inequalities.
   if (present(Ageq)) then
      !> Check dimensions.
      call assert(assertion=size(Ageq, 1) == size(bgeq), &
                  description="Ageq and bgeq have an inconsistent number of rows.")
      call assert(assertion=size(Ageq, 2) == n, &
                  description="Number of columns of Ageq is inconsistent with the number of variables (size(c)).")
      call assert(assertion=all(bgeq >= -eps), &
                  description="Right-hand side vector bgeq needs to be non-negative.")

      !> If all good, allocate arrays.
      problem%Ageq = Ageq; problem%bgeq = bgeq
   end if

   !> Consistency of the == inequalities.
   if (present(Aeq)) then
      !> Check dimensions.
      call assert(assertion=size(Aeq, 1) == size(beq), &
                  description="Aeq and beq have an inconsistent number of rows.")
      call assert(assertion=size(Aeq, 2) == n, &
                  description="Number of columns of Aeq is inconsistent with the number of variables (size(c)).")
      call assert(assertion=all(beq >= -eps), &
                  description="Right-hand side vector beq needs to be non-negative.")

      !> If all good, allocate arrays.
      problem%Aeq = Aeq; problem%beq = beq
   end if

   end procedure create_dense_linear_program
end module lightconvex_lp
