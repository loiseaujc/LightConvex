module TestLinearPrograms
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use stdlib_linalg_constants, only: ilp, dp, lk
   use LightConvex, only: simplex
   implicit none(external)
   private

   public collect_dense_simplex_problems
contains
   subroutine collect_dense_simplex_problems(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("Num. recipes example", test_num_recipes_problem)]
   end subroutine collect_dense_simplex_problems

   subroutine test_num_recipes_problem(error)
      type(error_type), allocatable, intent(out) :: error
      integer(ilp), parameter :: m = 4, n = 4
      real(dp), dimension(m + 2, n + 1) :: A
      integer(ilp), parameter :: nleq = 2, ngeq = 1, neq = 1
      real(dp) :: x(n), s(m), cost, cost_ref
      real(dp) :: xref(n), sref(m)
      integer(ilp) :: iposv(m), info, i, j

      !> Initialize the simplex tableau.
      A(1, :) = [0.0_dp, 1.0_dp, 1.0_dp, 3.0_dp, -0.5_dp]
      A(2, :) = [740.0_dp, -1.0_dp, 0.0_dp, -2.0_dp, 0.0_dp]
      A(3, :) = [0.0_dp, 0.0_dp, -2.0_dp, 0.0_dp, 7.0_dp]
      A(4, :) = [0.5_dp, 0.0_dp, -1.0_dp, 1.0_dp, -2.0_dp]
      A(5, :) = [9.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp]

      !> Solve the problem using the simplex method.
      block
         real(dp) :: start_time, end_time
         call cpu_time(start_time)
         call simplex(A, nleq, ngeq, neq, info, iposv)
         call cpu_time(end_time)
         print *, "     - Running time :", (end_time - start_time)*1000.0_dp, "milliseconds."
      end block

      !> Extract the primal and slack variables.
      x = 0.0_dp; s = 0.0_dp; cost = A(1, 1)
      do i = 1, m
         if (iposv(i) <= n) x(iposv(i)) = A(i + 1, 1)
         if (iposv(i) > n) s(iposv(i) - n) = A(i + 1, 1)
      end do

      !> Reference solution.
      cost_ref = 17.03_dp
      xref = [0.0_dp, 3.33_dp, 4.74_dp, 0.95_dp]
      sref = [730.55_dp, 0.0_dp, 0.0_dp, 0.0_dp]

      call check(error, info == 0)                          ! Optimal solution found.
      if (allocated(error)) return
      call check(error, maxval(abs(x - xref)) < 0.05_dp)    ! Matching primal.
      if (allocated(error)) return
      call check(error, maxval(abs(s - sref)) < 0.05_dp)    ! Matching slack.
      if (allocated(error)) return
      call check(error, abs(cost - cost_ref) < 0.05_dp)     ! Matching cost.
      if (allocated(error)) return

   end subroutine test_num_recipes_problem
end module TestLinearPrograms
