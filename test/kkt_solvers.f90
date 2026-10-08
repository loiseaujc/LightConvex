module TestKKTSolvers
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use stdlib_stats_distribution_normal, only: rvs_normal
   use stdlib_math, only: all_close, is_close
   use stdlib_linalg, only: eye, solve, diag, norm, mnorm
   use lightconvex_constants, only: ilp, dp
   use lightconvex, only: dense_vector, dense_matrix, dense_sym_matrix, &
                          dense_kkt_solver, kkt_info, kkt_solver, is_successful
   implicit none(external)
   private

   public :: collect_dense_kkt_solvers_tests

   real(dp), parameter :: atol = 1.0e-12_dp
   real(dp), parameter :: rtol = sqrt(atol)
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
         kkt = kkt_solver(dense_sym_matrix(P), dense_matrix(A))
         call kkt%update(d1, d2, info)
         call check(error, is_successful(info)) ! Successfully initialized the KKT solver.
         if (allocated(error)) return

         !> Solve the problem.
         call kkt%solve(rhs_x, rhs_y, x, y, info)
         call check(error, is_successful(info)) ! Successfully computed the solution.
         if (allocated(error)) return

         !> Check primal solution.
         call check(error, all_close(xref, x%data, abs_tol=atol))
         if (allocated(error)) return

         !> Check dual solution.
         call check(error, all_close(yref, y%data, abs_tol=atol))
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
         kkt = kkt_solver(dense_sym_matrix(P), dense_matrix(A))
         call kkt%update(d1, d2, info)
         call check(error, is_successful(info))
         if (allocated(error)) return

         !> Solve the problem.
         call kkt%solve(rhs_x, rhs_y, x, y, info)
         call check(error, is_successful(info))
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
         real(dp), parameter :: yref(m - 1) = [-0.5_dp, -0.5_dp]            ! Reference dual solution.
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
         kkt = kkt_solver(dense_sym_matrix(P), dense_matrix(A))
         call kkt%update(d1, d2, info, reg1=rtol, reg2=rtol)
         call check(error, is_successful(info)) ! Successfully initialized the KKT solver.
         if (allocated(error)) return

         !> Solve the problem.
         call kkt%solve(rhs_x, rhs_y, x, y, info)
         call check(error, is_successful(info)) ! Successfully computed the solution.
         if (allocated(error)) return           ! redundant constraint.

         !> Check primal solution.
         call check(error, all_close(xref, x%data, abs_tol=atol))
         if (allocated(error)) return

         !> Check dual solution.
         call check(error, &
                    all_close(matmul(transpose(A(:2, :)), yref), matmul(transpose(A), y%data), &
                              abs_tol=atol))
         if (allocated(error)) return
      end block

      !-----------------------------------
      !-----     RANDOM PROBLEMS     -----
      !-----------------------------------
      block
         integer(ilp), parameter :: nproblems = 1024, max_n = 128
         integer(ilp) :: i, j, m, n
         real(dp) :: u
         real(dp), allocatable :: P(:, :), q(:), A(:, :), b(:)
         type(dense_vector), allocatable :: x, y, rhs_x, rhs_y, d1, d2
         type(dense_kkt_solver), allocatable :: kkt
         type(kkt_info) :: info

         do i = 1, nproblems
            !> Random problem size.
            call random_number(u)
            n = 1_ilp + floor(max_n*u, kind=ilp)
            call random_number(u)
            m = 1_ilp + floor((n - 1)*u, kind=ilp) ! Ensure max(m) = n - 1.

            !> Allocate data.
            allocate (P(n, n), q(n), A(m, n), b(m), source=0.0_dp)

            !> Random problem.

            do j = 1, n
               P(:, j) = rvs_normal(array_size=n)
               A(:, j) = rvs_normal(array_size=m)
            end do
            P = matmul(P, transpose(P))
            q = rvs_normal(array_size=n)
            b = rvs_normal(array_size=m)

            x = dense_vector(n); rhs_x = dense_vector(q); d1 = dense_vector(n)
            y = dense_vector(m); rhs_y = dense_vector(b); d2 = dense_vector(m)

            !> Create KKT solver.
            kkt = kkt_solver(dense_sym_matrix(P), dense_matrix(A))
            call kkt%update(d1, d2, info, reg1=rtol, reg2=rtol)
            call check(error, is_successful(info))
            if (allocated(error)) return

            !> Solve the problem.
            call kkt%solve(rhs_x, rhs_y, x, y, info)
            call check(error, is_successful(info) .or. info%residual <= rtol)
            if (allocated(error)) then
               print *, "Problem's dimensions:"
               print *, "     - # of variables   :", n
               print *, "     - # of constraints :", m
               print *, "Residual of KKT solver  :", info%residual
               return
            end if

            !> Deallocate data.
            deallocate (P, q, A, b)
         end do
      end block

      !----------------------------------------------------------------------------
      !------      EDGE CASES : UNCONSTRAINED SYSTEM WITH SINGLE VARIABLE     -----
      !----------------------------------------------------------------------------
      block
         integer(ilp), parameter :: n = 1, m = 0
         integer(ilp) :: i, j
         real(dp) :: u
         real(dp), allocatable :: P(:, :), q(:), A(:, :), b(:), xref(:)
         type(dense_vector), allocatable :: x, y, rhs_x, rhs_y, d1, d2
         type(dense_kkt_solver), allocatable :: kkt
         type(kkt_info) :: info

         !> Allocate data.
         allocate (P(n, n), q(n), A(m, n), b(m), source=0.0_dp)

         !> Random problem.

         do j = 1, n
            P(:, j) = rvs_normal(array_size=n)
         end do
         P = matmul(P, transpose(P))
         q = rvs_normal(array_size=n)

         x = dense_vector(n); rhs_x = dense_vector(q); d1 = dense_vector(n)
         y = dense_vector(m); rhs_y = dense_vector(b); d2 = dense_vector(m)

         !> Create KKT solver.
         kkt = kkt_solver(dense_sym_matrix(P), dense_matrix(A))
         call kkt%update(d1, d2, info, reg1=rtol, reg2=rtol)
         call check(error, is_successful(info))
         if (allocated(error)) return

         !> Solve the problem.
         call kkt%solve(rhs_x, rhs_y, x, y, info)
         call check(error, is_successful(info) .or. info%residual <= rtol)
         if (allocated(error)) then
            print *, "Problem's dimensions:"
            print *, "     - # of variables   :", n
            print *, "     - # of constraints :", m
            print *, "Residual of KKT solver  :", info%residual
            return
         end if

         !> Check solution.
         xref = solve(P, q)
         call check(error, all_close(xref, x%data, abs_tol=atol))
         if (allocated(error)) return
      end block

      !-------------------------------------------------
      !-----     EDGE CASE : SINGULAR P MATRIX     -----
      !-------------------------------------------------
      block
         integer(ilp), parameter :: n = 128, rk = 64, m = n - rk + 1
         integer(ilp) :: i, j
         real(dp) :: u
         real(dp), allocatable :: P(:, :), q(:), A(:, :), b(:)
         type(dense_vector), allocatable :: x, y, rhs_x, rhs_y, d1, d2
         type(dense_kkt_solver), allocatable :: kkt
         type(kkt_info) :: info

         !> Allocate data.
         allocate (P(n, n), q(n), A(m, n), b(m), source=0.0_dp)

         !> Random problem.

         do j = 1, rk
            P(:, j) = rvs_normal(array_size=n)
         end do
         do j = 1, m
            A(j, :) = rvs_normal(array_size=n)
         end do
         P = matmul(P, transpose(P))
         q = rvs_normal(array_size=n)
         b = rvs_normal(array_size=m)

         x = dense_vector(n); rhs_x = dense_vector(q); d1 = dense_vector(n)
         y = dense_vector(m); rhs_y = dense_vector(b); d2 = dense_vector(m)

         !> Create KKT solver.
         kkt = kkt_solver(dense_sym_matrix(P), dense_matrix(A))
         call kkt%update(d1, d2, info)
         call check(error, is_successful(info))
         if (allocated(error)) return

         !> Solve the problem.
         call kkt%solve(rhs_x, rhs_y, x, y, info)
         call check(error, is_successful(info) .or. info%residual <= rtol)
         if (allocated(error)) then
            print *, "Problem's dimensions:"
            print *, "     - # of variables   :", n
            print *, "     - # of constraints :", m
            print *, "Residual of KKT solver  :", info%residual
            return
         end if
      end block

      !-------------------------------------------------------
      !------      EDGE CASE : SINGULAR KKT MATRIX      ------
      !-------------------------------------------------------
      block
         integer(ilp), parameter :: n = 128, m = 10
         integer(ilp) :: i, j
         real(dp) :: u
         real(dp), allocatable :: P(:, :), q(:), A(:, :), b(:)
         type(dense_vector), allocatable :: x, y, rhs_x, rhs_y, d1, d2
         type(dense_kkt_solver), allocatable :: kkt
         type(kkt_info) :: info

         !> Allocate data.
         allocate (P(n, n), q(n), A(m, n), b(m), source=0.0_dp)

         !> Random problem.
         call random_number(q); q(n) = 0.0_dp
         P = diag(q)

         do j = 1, m
            A(j, :) = rvs_normal(array_size=n)
         end do
         A(:, n) = 0.0_dp
         q = rvs_normal(array_size=n); q(n) = 0.0_dp
         b = rvs_normal(array_size=m)

         x = dense_vector(n); rhs_x = dense_vector(q); d1 = dense_vector(n)
         y = dense_vector(m); rhs_y = dense_vector(b); d2 = dense_vector(m)

         !> Create KKT solver.
         kkt = kkt_solver(dense_sym_matrix(P), dense_matrix(A))
         call kkt%update(d1, d2, info, reg1=1e-6_dp, reg2=1e-6_dp)
         call check(error, is_successful(info))
         if (allocated(error)) return

         !> Solve the problem.
         call kkt%solve(rhs_x, rhs_y, x, y, info)
         call check(error, is_successful(info) .or. info%residual <= rtol)
         if (allocated(error)) then
            print *, "Problem's dimensions:"
            print *, "     - # of variables   :", n
            print *, "     - # of constraints :", m
            print *, "Residual of KKT solver  :", info%residual
            return
         end if
      end block

      !-------------------------------------------------------------------
      !------     EDGE CASE : ADMM/IPM LIKE DIAGONALS WITH M > N     -----
      !-------------------------------------------------------------------
      block
         integer(ilp), parameter :: n = 10, m = 100
         integer(ilp) :: i, j
         real(dp) :: u
         real(dp), allocatable :: P(:, :), q(:), A(:, :), b(:)
         type(dense_vector), allocatable :: x, y, rhs_x, rhs_y, d1, d2
         type(dense_kkt_solver), allocatable :: kkt
         type(kkt_info) :: info

         !> Allocate data.
         allocate (P(n, n), q(n), A(m, n), b(m), source=0.0_dp)

         !> Random problem.

         do j = 1, n
            P(:, j) = rvs_normal(array_size=n)
            A(:, j) = rvs_normal(array_size=m)
         end do

         P = matmul(P, transpose(P)); P = P/mnorm(P, 2); A = A/mnorm(A, 2)
         q = rvs_normal(array_size=n); q = q/norm(q, 2)
         b = rvs_normal(array_size=m); b = b/norm(b, 2)

         x = dense_vector(n); rhs_x = dense_vector(q); d1 = dense_vector(n)
         y = dense_vector(m); rhs_y = dense_vector(b); d2 = dense_vector(m)

         call random_number(d1%data); d1%data = 10.0_dp**(8.0_dp*d1%data - 4.0_dp) ! 1e-4 .. 1e4
         call random_number(d2%data); d2%data = 10.0_dp**(8.0_dp*d2%data - 4.0_dp)

         !> Create KKT solver.
         kkt = kkt_solver(dense_sym_matrix(P), dense_matrix(A))
         call kkt%update(d1, d2, info, reg1=1.0e-8_dp, reg2=1.0e-8_dp)
         call check(error, is_successful(info))
         if (allocated(error)) return

         !> Solve the problem.
         call kkt%solve(rhs_x, rhs_y, x, y, info)
         call check(error, is_successful(info) .or. info%residual <= rtol)
         if (allocated(error)) then
            print *, "Problem's dimensions:"
            print *, "     - # of variables   :", n
            print *, "     - # of constraints :", m
            print *, "Residual of KKT solver  :", info%residual
            return
         end if
      end block

   end subroutine test_dense_kkt_solver_testsuite
end module TestKKTSolvers
