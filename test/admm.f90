module TestADMM
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use stdlib_stats_distribution_normal, only: rvs_normal
   use stdlib_math, only: all_close, is_close
   use stdlib_linalg, only: eye, solve, diag, norm, mnorm
   use lightconvex_constants, only: ilp, dp
   use lightconvex, only: dense_vector, dense_matrix, dense_sym_matrix, &
                          dense_kkt_solver, kkt_info, kkt_solver, is_successful, &
                          ADMM, admm_solver, qp_problem, qp_solution
   implicit none(external)
   private

   public :: collect_admm_tests

   real(dp), parameter :: atol = 1.0e-12_dp
   real(dp), parameter :: rtol = sqrt(atol)
contains
   subroutine collect_admm_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("ADMM dense solver testsuite", test_dense_admm_testsuite)]
   end subroutine collect_admm_tests

   subroutine test_dense_admm_testsuite(error)
      type(error_type), allocatable, intent(out) :: error

      !--------------------------------------------------
      !-----     SMALL SCALE LEAST-NORM PROBLEM     -----
      !--------------------------------------------------
      block
         integer(ilp), parameter :: m = 2, n = 4
         real(dp) :: P(n, n), q(n), A(m, n), b(m)
         real(dp), parameter :: xref(n) = [1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp]
         real(dp), parameter :: yref(m) = [0.5_dp, 0.5_dp]
         real(dp), allocatable :: prim_res(:), dual_res(:)
         type(dense_vector), allocatable :: x, y, z
         type(dense_kkt_solver), allocatable :: kkt

         !> Problem's matrices.
         P = eye(n, mold=1.0_dp)
         A(1, :) = [2.0_dp, -1.0_dp, 1.0_dp, -1.0_dp]
         A(2, :) = [0.0_dp, 1.0_dp, -1.0_dp, -1.0_dp]
         b = [3.0_dp, 1.0_dp]

         !> Workspace.
         x = dense_vector(xref)
         y = dense_vector(m)
         z = dense_vector(m)

         !> KKT solver.
         kkt = kkt_solver(dense_sym_matrix(P), dense_matrix(A))

         !> Solve the problem with ADMM.
         call admm_solver(ADMM(maxiter=100_ilp, check_every=1_ilp), &
                          dense_sym_matrix(P), dense_vector(q), &
                          dense_matrix(A), dense_vector(b), dense_vector(b), kkt, &
                          x, y, z, prim_res, dual_res)

         !> Check primal solution.
         print *, maxval(abs(xref - x%data))
         call check(error, all_close(xref, x%data, abs_tol=1.0e-2_dp))
         if (allocated(error)) return

      end block
   end subroutine test_dense_admm_testsuite
end module TestADMM
