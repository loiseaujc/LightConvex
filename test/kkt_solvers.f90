module TestKKTSolvers
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use stdlib_math, only: all_close, is_close
   use stdlib_linalg, only: norm, eye, solve
   use lightconvex_constants, only: ilp, dp
   use lightconvex, only: dense_vector, dense_kkt_solver, kkt_info, kkt_solver, is_successful
   implicit none(external)
   private

   public :: collect_dense_kkt_solvers_tests
contains
   subroutine collect_dense_kkt_solvers_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("Dense KKT solver testsuite", test_dense_kkt_solver_testsuite)]
   end subroutine collect_dense_kkt_solvers_tests

   subroutine test_dense_kkt_solver_testsuite(error)
      type(error_type), allocatable, intent(out) :: error

      !--------------------------------------------------
      !-----     SMALL SCALE LEAST-NORM PROBLEM     -----
      !--------------------------------------------------
      block
         integer(ilp), parameter :: m = 2, n = 4
         real(dp) :: P(n, n), A(m, n), b(m)
         real(dp), parameter :: xref(n) = [1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp] ! Reference primal solution.
         real(dp), parameter :: yref(m) = [-0.5_dp, -0.5_dp]                ! Reference dual solution.
         type(dense_vector), allocatable :: x, y, rhs_x, rhs_y, d1, d2
         type(dense_kkt_solver), allocatable :: kkt
         type(kkt_info) :: info

         !> Problem's matrices.
         P = eye(n, mold=1.0_dp)
         A(1, :) = [2.0_dp, -1.0_dp, 1.0_dp, -1.0_dp]
         A(2, :) = [0.0_dp, 1.0_dp, -1.0_dp, -1.0_dp]
         b = [3.0_dp, 1.0_dp]

         !> Problem's vectors.
         x = dense_vector(n); d1 = dense_vector(n)
         y = dense_vector(m); d2 = dense_vector(m)
         rhs_x = dense_vector(n); rhs_y = dense_vector(b)

         !> Create KKT solver.
         kkt = kkt_solver(P, A)
         call kkt%update(d1, d2, info)
         call check(error, is_successful(info)) ! Successfully initialized the KKT solver.
         if (allocated(error)) return

         !> Solve the problem.
         call kkt%solve(rhs_x, rhs_y, x, y, info)
         call check(error, is_successful(info)) ! Successfully computed the solution.
         call check(error, info%n_refine == 0)  ! No iterative refinement (no regularization was used).
         if (allocated(error)) return

         !> Check primal solution.
         call check(error, all_close(xref, x%data))
         if (allocated(error)) return

         !> Check dual solution.
         call check(error, all_close(yref, y%data))
         if (allocated(error)) return
      end block

      !-----------------------------------------------------------------------
      !-----     SMALL SCALE EQUALITY-CONSTRAINED STRICTLY CONVEX QP     -----
      !-----------------------------------------------------------------------
      block
         integer(ilp), parameter :: m = 3, n = 5
         real(dp) :: P(n, n), q(n)                              ! Quadratic cost.
         real(dp) :: A(m, n), b(m)                              ! Equality constraints.
         real(dp) :: xref(n), yref(m)                           ! Reference primal and dual solutions.
         real(dp), parameter :: obj_ref = 0.0393038880729888_dp ! Reference objective value.
         real(dp) :: obj
         type(dense_vector), allocatable :: x, y, rhs_x, rhs_y, d1, d2
         type(dense_kkt_solver), allocatable :: kkt
         type(kkt_info) :: info

         !> Problem's matrices.
         P = eye(n, mold=1.0_dp)
         q = [0.73727161_dp, 0.75526241_dp, 0.04741426_dp, -0.11260887_dp, -0.11260887_dp]
         A(:, 1) = [3.6_dp, 0.0_dp, -9.72_dp]
         A(:, 2) = [-3.4_dp, -1.9_dp, -8.67_dp]
         A(:, 3) = [-3.8_dp, -1.7_dp, 0.0_dp]
         A(:, 4) = [1.6_dp, -4.0_dp, 0.0_dp]
         A(:, 5) = [1.6_dp, -4.0_dp, 0.0_dp]
         b = [1.02_dp, 0.03_dp, 0.081_dp]

         !> Problem's vectors.
         x = dense_vector(n); d1 = dense_vector(n)
         y = dense_vector(m); d2 = dense_vector(m)
         rhs_x = dense_vector(q); rhs_y = dense_vector(b)

         !> Create KKT solver.
         kkt = kkt_solver(P, A)
         call kkt%update(d1, d2, info)
         call check(error, is_successful(info))
         if (allocated(error)) return

         !> Solve the problem.
         call kkt%solve(rhs_x, rhs_y, x, y, info)
         call check(error, is_successful(info))
         call check(error, info%n_refine == 0)  ! No need for regularization.
         if (allocated(error)) return

         !> Check primal solution.
         xref = [0.07313507_dp, -0.09133482_dp, -0.08677699_dp, 0.03638213_dp, 0.03638213_dp]
         call check(error, all_close(xref, x%data, abs_tol=1.0e-6_dp))
         if (allocated(error)) return

         !> Check objective function.
         obj = 0.5_dp*dot_product(x%data, matmul(P, x%data)) - dot_product(x%data, q)
         call check(error, is_close(obj, obj_ref, abs_tol=1.0e-6_dp))
         if (allocated(error)) return

         !> Check dual solution.
         yref = [-0.0440876_dp, 0.01961271_dp, -0.08465554_dp]
         call check(error, all_close(yref, y%data, abs_tol=1.0e-6_dp))
         if (allocated(error)) return
      end block

      !---------------------------------------------------------------------------
      !-----     SMALL SCALE LEAST-NORM PROBLEM WITH DUPLICATED CONSTRAINT    -----
      !----------------------------------------------------------------------------
      block
         integer(ilp), parameter :: m = 3, n = 4
         real(dp) :: P(n, n), A(m, n), b(m)
         real(dp), parameter :: xref(n) = [1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp] ! Reference primal solution.
         real(dp), parameter :: yref(m - 1) = [-0.5_dp, -0.5_dp]                ! Reference dual solution.
         type(dense_vector), allocatable :: x, y, rhs_x, rhs_y, d1, d2
         type(dense_kkt_solver), allocatable :: kkt
         type(kkt_info) :: info

         !> Problem's matrices.
         P = eye(n, mold=1.0_dp)
         A(1, :) = [2.0_dp, -1.0_dp, 1.0_dp, -1.0_dp]
         A(2, :) = [0.0_dp, 1.0_dp, -1.0_dp, -1.0_dp]
         A(3, :) = [0.0_dp, 1.0_dp, -1.0_dp, -1.0_dp]
         b = [3.0_dp, 1.0_dp, 1.0_dp]

         !> Problem's vectors.
         x = dense_vector(n); d1 = dense_vector(n)
         y = dense_vector(m); d2 = dense_vector(m)
         rhs_x = dense_vector(n); rhs_y = dense_vector(b)

         !> Create KKT solver.
         kkt = kkt_solver(P, A)
         call kkt%update(d1, d2, info, reg1=1e-6_dp, reg2=1e-6_dp)
         call check(error, is_successful(info)) ! Successfully initialized the KKT solver.
         if (allocated(error)) return

         !> Solve the problem.
         call kkt%solve(rhs_x, rhs_y, x, y, info)
         call check(error, is_successful(info)) ! Successfully computed the solution.
         call check(error, info%n_refine > 0)   ! Iterative refinement used to handle the
         if (allocated(error)) return           ! redundant constraint.

         !> Check primal solution.
         call check(error, all_close(xref, x%data, abs_tol=epsilon(1.0_dp)))
         if (allocated(error)) return

         ! !> Check dual solution.
         ! print *, "yref :", yref
         ! print *, "y    :", y%data
         ! call check(error, all_close(yref, y%data(1:m - 1), abs_tol=epsilon(1.0_dp)))
         ! if (allocated(error)) return
      end block

   end subroutine test_dense_kkt_solver_testsuite
end module TestKKTSolvers
