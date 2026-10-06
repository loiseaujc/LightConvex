module lightconvex_kkt
   use stdlib_optval, only: optval
   use stdlib_linalg_lapack, only: sytrf, sytrs
   use lightconvex_constants, only: dp, ilp, lk
   use lightconvex_abstract, only: abstract_cvx_vector
   use lightconvex_dense_vectors, only: dense_vector
   implicit none(type, external)
   private

   integer(ilp), parameter, public :: kkt_not_initialized = -99_ilp
   integer(ilp), parameter, public :: kkt_success = 0_ilp
   integer(ilp), parameter, public :: kkt_not_converged = -1_ilp
   integer(ilp), parameter, public :: kkt_not_quasidefinite = -2_ilp
   integer(ilp), parameter, public :: kkt_numerical_error = -3_ilp

   !> Derived-type returned by all KKT solvers.
   type, public :: kkt_info
      private
      !> Status of the KKT solver.
      integer(ilp) :: status = kkt_not_initialized
      !> Number of iterations for iterative solvers (0 for direct solvers).
      integer(ilp) :: n_iter = 0_ilp
      !> Number of iterative refinement steps taken.
      integer(ilp) :: n_refine = 0_ilp
      !> Number of matrix-vector products (P, A', or A).
      integer(ilp) :: n_spmv = 0_ilp
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
      real(dp), allocatable :: z(:, :), r(:), dz(:)    ! Length n+m, allocated once. Working set.
      logical(lk) :: factorized = .false.
      integer(ilp) :: max_refine = 10_ilp
      real(dp) :: refine_tol = 1.0e-13_dp
   contains
      procedure, pass(self) :: update => dense_update
      procedure, pass(self) :: solve => dense_solve
   end type dense_kkt_solver

   !-------------------------------------------
   !-----   FACTORIES FOR KKT SOLVERS     -----
   !-------------------------------------------

   interface kkt_solver
      module function create_dense_kkt_solver(P, A, max_refine, refine_tol) result(solver)
         implicit none(type, external)
         real(dp), intent(in) :: P(:, :), A(:, :)
         integer(ilp), optional, intent(in) :: max_refine
         real(dp), optional, intent(in) :: refine_tol
         type(dense_kkt_solver), allocatable :: solver
      end function create_dense_kkt_solver
   end interface
   public :: kkt_solver

contains

   module procedure is_successful
   bool = info%status == kkt_success
   end procedure is_successful

   !------------------------------------
   !-----     DENSE KKT SOLVER     -----
   !------------------------------------

   module procedure create_dense_kkt_solver
   !> Allocate solver.
   allocate (solver)
   associate (n => size(P, 1), m => size(A, 1))
      !> Sanity checks.
      if (size(P, 2) /= n) error stop "dense_kkt_solver: P needs to be a square matrix."
      if (size(A, 2) /= n) error stop "dense_kkt_solver: A needs to have the same number of columns as P."
      solver%n = n; solver%m = m

      solver%max_refine = optval(solver%max_refine, max_refine)
      if (solver%max_refine < 0) error stop "dense_kkt_solver: max_refine needs to be positive."

      solver%refine_tol = optval(solver%refine_tol, refine_tol)
      if (solver%refine_tol < 0.0_dp) error stop "dense_kkt_solver: refine_tol needs to be positive."

      ! --------------------

      !> Allocate matrices.
      allocate (solver%P, source=P)
      allocate (solver%A, source=A)
      allocate (solver%K(n + m, n + m), source=0.0_dp)
      allocate (solver%z(n + m, 1), solver%r(n + m), solver%dz(n + m), source=0.0_dp)
      allocate (solver%d1(n), solver%d2(m), source=0.0_dp)
      allocate (solver%e1(n), solver%e2(m), source=0.0_dp)
      allocate (solver%ipiv(n + m), source=0_ilp)

      !> Workspace query.
      block
         character(1), parameter :: uplo = "L"
         integer(ilp), parameter :: lwork = -1_ilp
         integer(ilp) :: info
         real(dp) :: dummy_work(1)
         !> Workspace query.
         call sytrf(uplo, n + m, solver%K, n + m, solver%ipiv, dummy_work, lwork, info)
         if (info /= 0) error stop "dense_kkt_solver: error in sytrf."
         allocate (solver%workspace(int(dummy_work(1), kind=ilp)), source=0.0_dp)
      end block
   end associate
   end procedure create_dense_kkt_solver

   subroutine dense_update(self, d1, d2, info, reg1, reg2)
      implicit none(type, external)
      class(dense_kkt_solver), intent(inout) :: self
      class(abstract_cvx_vector), intent(in) :: d1
      class(abstract_cvx_vector), intent(in) :: d2
      type(kkt_info), intent(out) :: info
      real(dp), optional, intent(in) :: reg1, reg2

      self%factorized = .false.

      select type (d1)
      type is (dense_vector)
         select type (d2)
         type is (dense_vector)
            associate (n => d1%get_size(), m => d2%get_size(), uplo => "L")
               !> Sanity check.
               if (n /= self%n) error stop "kkt%update: d1 has inconsistent dimensions."
               if (m /= self%m) error stop "kkt%update: d2 has inconsistent dimensions."

               !> Regularized diagonals.
               self%e1 = d1%data + optval(reg1, 0.0_dp)
               self%e2 = d2%data + optval(reg2, 0.0_dp)

               if (any(self%e1 < 0.0_dp) .or. any(self%e2 < 0.0_dp)) then
                  info%status = kkt_not_quasidefinite
                  self%factorized = .false.
                  return
               end if

               !> Assembling the regularized KKT matrix.
               block
                  integer(ilp) :: i, j
                  do j = 1, n
                     !> Fill the (1, 1) block of K.
                     do i = 1, j
                        self%K(i, j) = self%P(i, j) ! Store only the lower triangular part.
                     end do
                     self%K(j, j) = self%K(j, j) + self%e1(j)   ! Regularized diagonal.

                     !> Fill the (2, 1) block of K.
                     do i = 1, m
                        self%K(n + i, j) = self%A(i, j)
                     end do
                  end do

                  !> Fill the (2, 2) diagonal block of K.
                  do j = 1, m
                     self%K(n + j, n + j) = -self%e2(j)
                  end do
               end block

               !> LDLT factorization of K.
               block
                  integer(ilp) :: lwork, lapack_info
                  self%ipiv = 0_ilp
                  lwork = size(self%workspace, kind=ilp)
                  call sytrf(uplo, n + m, self%K, n + m, self%ipiv, &
                             self%workspace, lwork, lapack_info)
                  if (lapack_info < 0) then
                     error stop "kkt%update: error in sytrf."
                  else if (lapack_info > 0) then
                     info%status = kkt_numerical_error
                     return
                  end if
               end block

               !> Store unregularized diagonal for iterative refinement.
               self%d1 = d1%data
               self%d2 = d2%data

               !> Update solver status.
               info%status = kkt_success
               self%factorized = .true.
            end associate
         class default
            error stop "kkt%update: Type of d2 should be dense_vector."
         end select
      class default
         error stop "kkt%update: Type of d1 should be dense_vector."
      end select
   end subroutine dense_update

   subroutine dense_solve(self, rhs_x, rhs_y, sol_x, sol_y, info, tol)
      implicit none(type, external)
      class(dense_kkt_solver), intent(inout) :: self
      class(abstract_cvx_vector), intent(in) :: rhs_x, rhs_y
      class(abstract_cvx_vector), intent(inout) :: sol_x, sol_y
      type(kkt_info), intent(out) :: info
      real(dp), optional, intent(in) :: tol

      if (.not. self%factorized) then
         error stop "kkt%solve: KKT matrix has not been factorized. Aborting."
      end if

      select type (rhs_x)
      type is (dense_vector)
         select type (rhs_y)
         type is (dense_vector)
            select type (sol_x)
            type is (dense_vector)
               select type (sol_y)
               type is (dense_vector)

                  associate (n => self%n, m => self%m, uplo => "L")
                     if (rhs_x%get_size() /= n) error stop "kkt%solve: rhs_x has inconsistent dimensions."
                     if (sol_x%get_size() /= n) error stop "kkt%solve: sol_x has inconsistent dimensions."
                     if (rhs_y%get_size() /= m) error stop "kkt%solve: rhs_y has inconsistent dimensions."
                     if (sol_y%get_size() /= m) error stop "kkt%solve: sol_y has inconsistent dimensions."

                     !> Working vector.
                     self%z(:n, 1) = rhs_x%data; self%z(n + 1:, 1) = rhs_y%data

                     !> Solve the linear system.
                     block
                        integer(ilp) :: lapack_info
                        call sytrs(uplo, n + m, 1, self%K, n + m, self%ipiv, self%z, n + m, lapack_info)
                        if (lapack_info /= 0) error stop "kkt%solve: Error in sytrs."
                     end block

                     !> Approximated solutions.
                     sol_x%data = self%z(:n, 1)
                     sol_y%data = self%z(n + 1:, 1)

                     block
                        integer(ilp) :: iter
                        !> Iterative refinement.
                        do iter = 1, self%max_refine
                           exit
                        end do

                        !> Book-keeping.
                        info%status = merge(kkt_success, kkt_not_converged, iter < self%max_refine)
                        info%n_iter = 0
                     end block
                  end associate

               class default
                  error stop "kkt%solve: Type of sol_y should be dense vector."
               end select
            class default
               error stop "kkt%solve: Type of sol_x should be dense vector."
            end select
         class default
            error stop "kkt%solve: Type of rhs_y should be dense_vector."
         end select
      class default
         error stop "kkt%solve: Type of rhs_x should be dense_vector."
      end select
   end subroutine dense_solve

end module lightconvex_kkt
