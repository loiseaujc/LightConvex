module lightconvex_kkt
   use stdlib_optval, only: optval
   use lightconvex_constants, only: dp, ilp, lk
   use lightconvex_abstract, only: abstract_cvx_vector, AbstractMatrix, AbstractSymMatrix
   implicit none(type, external)
   private

   integer(ilp), parameter, public :: kkt_not_initialized = -99_ilp
   integer(ilp), parameter, public :: kkt_success = 0_ilp
   integer(ilp), parameter, public :: kkt_not_converged = -1_ilp
   integer(ilp), parameter, public :: kkt_invalid_regularization = -2_ilp
   integer(ilp), parameter, public :: kkt_numerical_error = -3_ilp

   !> Derived-type returned by all KKT solvers.
   type, public :: kkt_info
      !> Status of the KKT solver.
      integer(ilp) :: status = kkt_not_initialized
      !> Number of iterations for iterative solvers (0 for direct solvers).
      integer(ilp) :: n_iter = 0_ilp
      !> Relative residual of the unregularized problem (if known).
      real(dp) :: residual = huge(1.0_dp)
   end type kkt_info

   interface
      logical(lk) pure module function is_successful(info) result(bool)
         implicit none(type, external)
         type(kkt_info), intent(in) :: info
      end function is_successful
   end interface
   public :: is_successful

   !> Abstract type for defining KKT solvers.
   type, public, abstract :: abstract_kkt_solver
   contains
      procedure(update_iface), pass(self), deferred :: update
      procedure(solve_iface), pass(self), deferred :: solve
   end type abstract_kkt_solver

   abstract interface
      subroutine update_iface(self, d1, d2, info, reg1, reg2)
         import abstract_cvx_vector, abstract_kkt_solver, dp, kkt_info
         implicit none(type, external)
         class(abstract_kkt_solver), intent(inout) :: self
         class(abstract_cvx_vector), intent(in) :: d1
         class(abstract_cvx_vector), intent(in) :: d2
         type(kkt_info), intent(out) :: info
         real(dp), optional, intent(in) :: reg1, reg2
      end subroutine update_iface

      subroutine solve_iface(self, rhs_x, rhs_y, sol_x, sol_y, info, tol)
         import abstract_cvx_vector, abstract_kkt_solver, dp, kkt_info
         implicit none(type, external)
         class(abstract_kkt_solver), intent(inout) :: self
         class(abstract_cvx_vector), intent(in) :: rhs_x, rhs_y
         class(abstract_cvx_vector), intent(inout) :: sol_x, sol_y
         type(kkt_info), intent(out) :: info
         real(dp), optional, intent(in) :: tol
      end subroutine solve_iface
   end interface

   !-------------------------------------------------
   !-----     KKT SOLVER FOR DENSE PROBLEMS     -----
   !-------------------------------------------------

   type, public, extends(abstract_kkt_solver) :: dense_kkt_solver
      private
      integer(ilp) :: n = 0_ilp, m = 0_ilp  ! Number of variables (n) and constraints (m).
      real(dp), allocatable :: P(:, :)      ! Quadratic form. n x n matrix.
      real(dp), allocatable :: A(:, :)      ! Constraint matrix. m x n matrix.
      real(dp), allocatable :: d1(:), d2(:) ! Unregularized diagonals.
      real(dp), allocatable :: e1(:), e2(:) ! Regularized diagonals.
      real(dp), allocatable :: K(:, :)      ! LDLT factors of regularized KKT matrix.
      integer(ilp), allocatable :: ipiv(:)  ! Pivots from LDLT.
      real(dp), allocatable :: workspace(:) ! LAPACK workspace (queried only once).
      real(dp), allocatable :: z(:, :), r(:), dz(:, :)    ! Length n+m, allocated once. Working set.
      logical(lk) :: factorized = .false.
      integer(ilp) :: max_refine = 20_ilp
      real(dp) :: tol = 1.0e-08_dp
   contains
      procedure, pass(self) :: update => dense_update
      procedure, pass(self) :: solve => dense_solve
   end type dense_kkt_solver

   interface
      module subroutine dense_update(self, d1, d2, info, reg1, reg2)
         implicit none(type, external)
         class(dense_kkt_solver), intent(inout) :: self
         class(abstract_cvx_vector), intent(in) :: d1
         class(abstract_cvx_vector), intent(in) :: d2
         type(kkt_info), intent(out) :: info
         real(dp), optional, intent(in) :: reg1, reg2
      end subroutine dense_update

      module subroutine dense_solve(self, rhs_x, rhs_y, sol_x, sol_y, info, tol)
         implicit none(type, external)
         class(dense_kkt_solver), intent(inout) :: self
         class(abstract_cvx_vector), intent(in) :: rhs_x, rhs_y
         class(abstract_cvx_vector), intent(inout) :: sol_x, sol_y
         type(kkt_info), intent(out) :: info
         real(dp), optional, intent(in) :: tol
      end subroutine dense_solve
   end interface

   !-------------------------------------------
   !-----   FACTORIES FOR KKT SOLVERS     -----
   !-------------------------------------------

   interface kkt_solver
      module function create_dense_kkt_solver(P, A, max_refine, tol) result(solver)
         implicit none(type, external)
         class(AbstractSymMatrix), intent(in) :: P
         class(AbstractMatrix), intent(in) :: A
         integer(ilp), optional, intent(in) :: max_refine
         real(dp), optional, intent(in) :: tol
         type(dense_kkt_solver), allocatable :: solver
      end function create_dense_kkt_solver
   end interface
   public :: kkt_solver

contains

   module procedure is_successful
   bool = info%status == kkt_success
   end procedure is_successful

end module lightconvex_kkt
