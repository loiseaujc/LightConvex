module TestLinearPrograms
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use stdlib_linalg_constants, only: ilp, dp, lk
   use LightConvex, only: simplex, &
                          Dantzig, &
                          auxiliary_function, &
                          infeasible_status, optimal_status, &
                          unbounded_status, maxiter_exceeded, &
                          tol
   implicit none(external)
   private

   public collect_dense_standard_simplex_problems
contains
   subroutine collect_dense_standard_simplex_problems(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("Num. recipes example", test_num_recipes_problem)]
      testsuite = [testsuite, new_unittest("Wiki example", test_wikipedia_example)]
      testsuite = [testsuite, new_unittest("Infeasible example", test_infeasible_lp)]
      testsuite = [testsuite, new_unittest("Unbounded example", test_unbounded_lp)]
      testsuite = [testsuite, new_unittest("Caltech examples", test_caltech_examples)]
   end subroutine collect_dense_standard_simplex_problems

   ! Test problem (10.8.6)-(10.8.7) from Numerical Recipes.
   ! The original LP reads
   !
   !   maximize    x1 + x2 + x3 - 0.5 x4
   !   subject to           x1 + 2 x3 <= 740
   !                     2 x2 - 7 x4 <= 0
   !                  x2 - x3 + 2 x4 >= 0.5
   !               x1 + x2 + x3 + x4  = 9
   !                  x1, x2, x3, x4 >= 0
   !
   ! The solution reported in the reference is given by
   ! x = [0.0, 3.33, 4.74, 0.95] (up to two significant digits)
   ! along with the slacks s = [730.55, 0, 0, 0]. The corresponding
   ! optimal cost is c = 17.03
   subroutine test_num_recipes_problem(error)
      type(error_type), allocatable, intent(out) :: error
      integer(ilp), parameter :: m = 4, n = 4, maxiter = 10
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
      call simplex(A, nleq, ngeq, neq, iposv, maxiter, info, &
                   Dantzig(), auxiliary_function())

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

   ! Test problem taken from wikipedia:
   !    https://en.wikipedia.org/wiki/Simplex_algorithm
   !
   ! The original LP reads
   !
   !   maximize    2x + 3y + 4z
   !   subject to  3x + 2y + z  <= 10
   !               2x + 5y + 3z <= 15
   !                    x, y, z >= 0
   !
   ! The solution is given by x = [0, 0, 5] with slack variables
   ! given by s = [5, 0]. Corresponding optimal cost is c = 20.
   subroutine test_wikipedia_example(error)
      type(error_type), allocatable, intent(out) :: error
      integer(ilp), parameter :: m = 2, n = 3, maxiter = 10
      real(dp), dimension(m + 2, n + 1) :: A
      integer(ilp), parameter :: nleq = 2, ngeq = 0, neq = 0
      real(dp) :: x(n), s(m), cost, cost_ref
      real(dp) :: xref(n), sref(m)
      integer(ilp) :: iposv(m), info, i, j

      !> Initialize the simplex tableau.
      A(1, :) = [0.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
      A(2, :) = [10.0_dp, -3.0_dp, -2.0_dp, -1.0_dp]
      A(3, :) = [15.0_dp, -2.0_dp, -5.0_dp, -3.0_dp]

      !> Solve the problem using the simplex method.
      call simplex(A, nleq, ngeq, neq, iposv, maxiter, info, &
                   Dantzig(), auxiliary_function())

      !> Extract the primal and slack variables.
      x = 0.0_dp; s = 0.0_dp; cost = A(1, 1)
      do i = 1, m
         if (iposv(i) <= n) x(iposv(i)) = A(i + 1, 1)
         if (iposv(i) > n) s(iposv(i) - n) = A(i + 1, 1)
      end do

      !> Reference solution.
      cost_ref = 20.0_dp
      xref = [0.0_dp, 0.0_dp, 5.0_dp]
      sref = [5.0_dp, 0.0_dp]

      call check(error, info == 0)              ! Optimal solution found.
      if (allocated(error)) return
      call check(error, maxval(abs(x - xref)) <= tol)  ! Matching primal.
      if (allocated(error)) return
      call check(error, maxval(abs(s - sref)) <= tol)  ! Matching slack.
      if (allocated(error)) return
      call check(error, abs(cost - cost_ref) <= tol)   ! Matching cost.
      if (allocated(error)) return

   end subroutine test_wikipedia_example

   ! This is a simple example illustrating an infeasible problem.
   ! The corresponding LP reads
   !
   !   maximize    x + y
   !   subject to  -x -y >= 1
   !                x, y >= 0
   !
   ! which obviously has no solution.
   subroutine test_infeasible_lp(error)
      type(error_type), allocatable, intent(out) :: error
      integer(ilp), parameter :: m = 1, n = 2, maxiter = 10
      real(dp), dimension(m + 2, n + 1) :: A
      integer(ilp), parameter :: nleq = 0, ngeq = 1, neq = 0
      integer(ilp) :: iposv(m), info, i, j

      !> Initialize the simplex tableau.
      A(1, :) = [0.0_dp, 1.0_dp, 1.0_dp]
      A(2, :) = [1.0_dp, 1.0_dp, 1.0_dp]

      !> Solve the problem using the simplex method.
      call simplex(A, nleq, ngeq, neq, iposv, maxiter, info, &
                   Dantzig(), auxiliary_function())

      call check(error, info == -1)              ! Problem is infeasible.
      if (allocated(error)) return

   end subroutine test_infeasible_lp

   ! This is a simple example illustrating an unbounded problem.
   ! The corresponding LP reads
   !
   !   maximize     x + y
   !   subject to   x + y >= 5
   !                x, y >= 0
   !
   ! which clearly has no finite solution.
   subroutine test_unbounded_lp(error)
      type(error_type), allocatable, intent(out) :: error
      integer(ilp), parameter :: m = 1, n = 2, maxiter = 10
      real(dp), dimension(m + 2, n + 1) :: A
      integer(ilp), parameter :: nleq = 0, ngeq = 1, neq = 0
      integer(ilp) :: iposv(m), info, i, j

      !> Initialize the simplex tableau.
      A(1, :) = [0.0_dp, 1.0_dp, 1.0_dp]
      A(2, :) = [5.0_dp, -1.0_dp, -1.0_dp]

      !> Solve the problem using the simplex method.
      call simplex(A, nleq, ngeq, neq, iposv, maxiter, info, &
                   Dantzig(), auxiliary_function())

      call check(error, info == 1)              ! Objective is unbounded.
      if (allocated(error)) return

   end subroutine test_unbounded_lp

   subroutine test_caltech_examples(error)
      type(error_type), allocatable, intent(out) :: error

      ! Lists of problems taken from
      !     K. C. Border, "The Gauss-Jordan and Simplex Algorithms",
      !     Caltech Division of the Humanities and Social Sciences, 2004

      ! The problem reads
      !
      !   maximize    2 x1 + 4 x2 + x3 + x4
      !   subject to  2 x1 + x2 <= 3
      !               x2 + 4 x3 + x4 <= 3
      !               x1 + 3 x2 + x4 <= 4
      !               x1, x2, x3, x4 >= 0
      !
      ! The primal solution is given by x = [1, 1, 0.5, 0] with cost c.T @ x = 6.5
      block
         integer(ilp), parameter :: m = 3, n = 4, maxiter = 10
         real(dp), dimension(m + 2, n + 1) :: A
         integer(ilp), parameter :: nleq = 3, ngeq = 0, neq = 0
         integer(ilp) :: iposv(m), info, i, j
         real(dp) :: x(n), cost
         real(dp) :: xref(n), cost_ref

         !> Simplex tableau.
         A(1, :) = [0.0_dp, 2.0_dp, 4.0_dp, 1.0_dp, 0.0_dp]
         A(2, :) = [3.0_dp, -2.0_dp, -1.0_dp, 0.0_dp, 0.0_dp]
         A(3, :) = [3.0_dp, 0.0_dp, -1.0_dp, -4.0_dp, -1.0_dp]
         A(4, :) = [4.0_dp, -1.0_dp, -3.0_dp, 0.0_dp, -1.0_dp]

         !> Solve the problem using the simplex method.
         call simplex(A, nleq, ngeq, neq, iposv, maxiter, info, &
                      Dantzig(), auxiliary_function())

         !> Extract the primal and slack variables.
         x = 0.0_dp; cost = A(1, 1)
         do i = 1, m
            if (iposv(i) <= n) x(iposv(i)) = A(i + 1, 1)
         end do

         !> Reference solution.
         cost_ref = 6.5_dp
         xref = [1.0_dp, 1.0_dp, 0.5_dp, 0.0_dp]

         call check(error, info == 0)              ! Optimal solution found.
         if (allocated(error)) return
         call check(error, maxval(abs(x - xref)) <= tol)  ! Matching primal.
         if (allocated(error)) return
         call check(error, abs(cost - cost_ref) <= tol)   ! Matching cost.
         if (allocated(error)) return
      end block

      ! The problem reads
      !
      !   maximize    3/4 x1 - 150 x2 + 1/50 x3 - 6 x4
      !   subject to  1/4 x1 - 60 x2 - 1/25 x3 + 9 x4 <= 0
      !               1/2 x1 - 90 x2 - 1/50 x3 + 3 x4 <= 0
      !                           x3 <= 1
      !               x1, x2, x3, x4 >= 0
      !
      ! Primal solution is x = [1/25, 0, 1, 0] with cost c.T @ x = 1/20
      ! Naive implementation of the pivot rule leads to cycling behavior.
      block
         integer(ilp), parameter :: m = 3, n = 4, maxiter = 10
         real(dp), dimension(m + 2, n + 1) :: A
         integer(ilp), parameter :: nleq = 3, ngeq = 0, neq = 0
         integer(ilp) :: iposv(m), info, i, j
         real(dp) :: x(n), cost
         real(dp) :: xref(n), cost_ref

         !> Simplex tableau.
         A(1, :) = [0.0_dp, 3.0_dp/4.0_dp, -150.0_dp, 1.0_dp/50.0_dp, -6.0_dp]
         A(2, :) = [0.0_dp, -1.0_dp/4.0_dp, 60.0_dp, 1.0_dp/25_dp, -9.0_dp]
         A(3, :) = [0.0_dp, -0.5_dp, 90.0_dp, 1.0_dp/50.0_dp, -3.0_dp]
         A(4, :) = [1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 1.0_dp]

         !> Solve the problem using the simplex method.
         call simplex(A, nleq, ngeq, neq, iposv, maxiter, info, &
                      Dantzig(), auxiliary_function())

         !> Extract the primal and slack variables.
         x = 0.0_dp; cost = A(1, 1)
         do i = 1, m
            if (iposv(i) <= n) x(iposv(i)) = A(i + 1, 1)
         end do

         !> Reference solution.
         cost_ref = 1.0/20.0_dp
         xref = [1.0_dp/25.0_dp, 0.0_dp, 1.0_dp, 0.0_dp]

         call check(error, info == 0)              ! Optimal solution found.
         if (allocated(error)) return
         call check(error, maxval(abs(x - xref)) <= tol)  ! Matching primal.
         if (allocated(error)) return
         call check(error, abs(cost - cost_ref) <= tol)   ! Matching cost.
         if (allocated(error)) return
      end block

      ! The problem reads
      !
      !   maximize      -x1 - 6 x2 + 7 x3 - x4 - 5 x5
      !   subject to    5 x1 - 4 x2 + 13 x3 - 2 x4 + x5 = 20
      !                 x1 - x2 + 5 x3 - x4 + x5 = 8
      !                 x1, x2, x3, x4, x5 >= 0
      !
      ! Primal solution is x = [0, 4/7, 12/7, 0, 0] with cost c.T @ x = 60/7
      block
         integer(ilp), parameter :: m = 2, n = 5, maxiter = 10
         real(dp), dimension(m + 2, n + 1) :: A
         integer(ilp), parameter :: nleq = 0, ngeq = 0, neq = 2
         integer(ilp) :: iposv(m), info, i, j
         real(dp) :: x(n), cost
         real(dp) :: xref(n), cost_ref

         !> Simplex tableau.
         A(1, :) = [0.0_dp, -1.0_dp, -6.0_dp, 7.0_dp, -1.0_dp, -5.0_dp]
         A(2, :) = [20.0_dp, -5.0_dp, 4.0_dp, -13.0_dp, 2.0_dp, -1.0_dp]
         A(3, :) = [8.0_dp, -1.0_dp, 1.0_dp, -5.0_dp, 1.0_dp, -1.0_dp]

         !> Solve the problem using the simplex method.
         call simplex(A, nleq, ngeq, neq, iposv, maxiter, info, &
                      Dantzig(), auxiliary_function())

         !> Extract the primal and slack variables.
         x = 0.0_dp; cost = A(1, 1)
         do i = 1, m
            if (iposv(i) <= n) x(iposv(i)) = A(i + 1, 1)
         end do

         !> Reference solution.
         cost_ref = 60.0_dp/7.0_dp
         xref = [0.0_dp, 4.0_dp/7.0_dp, 12.0_dp/7.0_dp, 0.0_dp, 0.0_dp]

         call check(error, info == 0)              ! Optimal solution found.
         if (allocated(error)) return
         call check(error, maxval(abs(x - xref)) <= tol)  ! Matching primal.
         if (allocated(error)) return
         call check(error, abs(cost - cost_ref) <= tol)   ! Matching cost.
         if (allocated(error)) return
      end block

      !
      block
         integer(ilp), parameter :: m = 3, n = 4, maxiter = 10
         real(dp), dimension(m + 2, n + 1) :: A
         integer(ilp), parameter :: nleq = 1, ngeq = 2, neq = 0
         integer(ilp) :: iposv(m), info, i, j
         real(dp) :: x(n), cost
         real(dp) :: xref(n), cost_ref

         !> Simplex tableau.
         A(1, :) = [0.0_dp, 2.0_dp, -3.0_dp, 1.0_dp, 1.0_dp]
         A(2, :) = [3.0_dp, -1.0_dp, -2.0_dp, -1.0_dp, -1.0_dp]
         A(3, :) = [2.0_dp, 1.0_dp, -2.0_dp, 2.0_dp, 1.0_dp]
         A(4, :) = [1.0_dp, 3.0_dp, -1.0_dp, 0.0_dp, -1.0_dp]

         !> Solve the problem using the simplex method.
         call simplex(A, nleq, ngeq, neq, iposv, maxiter, info, &
                      Dantzig(), auxiliary_function())

         !> Extract the primal and slack variables.
         x = 0.0_dp; cost = A(1, 1)
         do i = 1, m
            if (iposv(i) <= n) x(iposv(i)) = A(i + 1, 1)
         end do

         !> Reference solution.
         cost_ref = -3.0_dp
         xref = [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]

         call check(error, info == 0)              ! Optimal solution found.
         if (allocated(error)) return
         call check(error, maxval(abs(x - xref)) <= tol)  ! Matching primal.
         if (allocated(error)) return
         call check(error, abs(cost - cost_ref) <= tol)   ! Matching cost.
         if (allocated(error)) return
      end block
   end subroutine test_caltech_examples
end module TestLinearPrograms
