module lightconvex_lp
   use assert_m, only: assert => assert_always
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
   public :: simplex

   !----- Pivot selection rules -----

   !> Default pivot selection rule.
   type, extends(abstract_pivot_rule), public :: Dantzig
   end type Dantzig

   !----- Initialization methods -----

   !> Auxiliary function maximization (default).
   type, extends(abstract_feasible_initialization), public :: auxiliary_function
   end type auxiliary_function

contains

   !==================================================
   !=====                                        =====
   !=====     DEFINITION OF A LINEAR PROGRAM     =====
   !=====                                        =====
   !==================================================

   !-----------------------------------------
   !-----     DENSE LINEAR PROGRAMS     -----
   !-----------------------------------------

   module procedure assemble_dense_lp
   integer(ilp) :: m, n
    !! Number of constraints and number of variables.
   integer(ilp) :: offset

   n = size(c)  ! Number of variables.

   !> Sanity checks.
   call assert(assertion=(present(Aleq) .and. present(bleq)) .or. (.not. present(Aleq) .and. .not. present(bleq)), &
               description="Specification of <= constraints incomplete. Either Aleq or bleq is missing.")
   call assert(assertion=(present(Ageq) .and. present(bgeq)) .or. (.not. present(Ageq) .and. .not. present(bgeq)), &
               description="Specification of >= constraints incomplete. Either Aleq or bleq is missing.")
   call assert(assertion=(present(Aeq) .and. present(beq)) .or. (.not. present(Aeq) .and. .not. present(beq)), &
               description="Specification of == constraints incomplete. Either Aleq or bleq is missing.")

   associate (nleq => problem%nleq, ngeq => problem%ngeq, neq => problem%neq)
      nleq = 0; ngeq = 0; neq = 0

      !> Consistency of the <= inequalities.
      if (present(Aleq)) then
         !> Check dimensions.
         call assert(assertion=size(Aleq, 1) == size(bleq), &
                     description="Aleq and bleq have an inconsistent number of rows.")
         call assert(assertion=size(Aleq, 2) == n, &
                     description="Number of columns of Aleq is inconsistent with the number of variables (size(c)).")
         call assert(assertion=all(bleq >= -eps), &
                     description="Right-hand side vector bleq needs to be non-negative.")
         !> Number of <= constraints.
         nleq = size(Aleq, 1)
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
         !> Number of >= constraints.
         ngeq = size(Ageq, 1)
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
         !> Number of <= constraints.
         neq = size(Aeq, 1)
      end if

      !> Total number of constraints.
      m = nleq + ngeq + neq

      !----- Construct the Simplex tableau -----
      allocate (problem%A(m + 2, n + 1), source=0.0_dp); problem%A(1, 2:) = c   ! Cost function.
      offset = 1

      !> Add the <= inequalities.
      if (present(Aleq)) then
         problem%A(offset + 1:nleq + offset, 2:) = -Aleq
         problem%A(offset + 1:nleq + offset, 1) = bleq
         ! offset = offset + 1
      end if

      !> Add the >= inequalities.
      if (present(Ageq)) then
         problem%A(nleq + offset + 1:nleq + ngeq + offset, 2:) = -Ageq
         problem%A(nleq + offset + 1:nleq + ngeq + offset, 1) = bgeq
         ! offset = offset + 1
      end if

      !> Add the == constraints.
      if (present(Aeq)) then
         problem%A(nleq + ngeq + offset + 1:m + 1, 2:) = -Aeq
         problem%A(nleq + ngeq + offset + 1:m + 1, 1) = beq
      end if
   end associate
   end procedure assemble_dense_lp
end module lightconvex_lp
