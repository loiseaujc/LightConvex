module TestLinearPrograms
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use stdlib_math, only: is_close, all_close
   use stdlib_linalg_constants, only: ilp, dp, lk
   use stdlib_math, only: all_close
   use LightConvex, only: dense_lp_type, linear_program, lp_solution, &
                          simplex, Dantzig, auxiliary_function, &
                          PrimalSimplex, solve, &
                          is_optimal, is_feasible, is_unbounded, &
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
      testsuite = [testsuite, new_unittest("Cornell example", test_cornell_example)]
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
      integer(ilp), parameter :: nleq = 2, ngeq = 1, neq = 1
      real(dp) :: xref(n), sref(m), cost_ref
      real(dp) :: c(n)
      real(dp) :: Aleq(nleq, n), bleq(nleq)
      real(dp) :: Ageq(ngeq, n), bgeq(ngeq)
      real(dp) :: Aeq(neq, n), beq(neq)
      type(dense_lp_type) :: problem
      type(lp_solution) :: solution

      !> Reference solution.
      cost_ref = 17.03_dp
      xref = [0.0_dp, 3.33_dp, 4.74_dp, 0.95_dp]
      sref = [730.55_dp, 0.0_dp, 0.0_dp, 0.0_dp]

      !> Cost vector.
      c = [1.0_dp, 1.0_dp, 3.0_dp, -0.5_dp]

      !> Linear <= inequalities.
      Aleq(1, :) = [1.0_dp, 0.0_dp, 2.0_dp, 0.0_dp]
      Aleq(2, :) = [0.0_dp, 2.0_dp, 0.0_dp, -7.0_dp]
      bleq = [740.0_dp, 0.0_dp]

      !> Linear >= inequalities.
      Ageq(1, :) = [0.0_dp, 1.0_dp, -1.0_dp, 2.0_dp]
      bgeq = 0.5_dp

      !> Equality constraints.
      Aeq(1, :) = 1.0_dp; beq = 9.0_dp

      !> Construct the LP problem for LightConvex.
      problem = linear_program(c, Aleq=Aleq, bleq=bleq, Ageq=Ageq, bgeq=bgeq, Aeq=Aeq, beq=beq)

      !> Solve the problem.
      solution = solve(problem, alg=PrimalSimplex())

      call check(error, is_optimal(problem))                                             ! Optimal solution found.
      if (allocated(error)) return
      call check(error, all_close(solution%x, xref, abs_tol=0.05_dp))                    ! Matching primal.
      if (allocated(error)) return
      call check(error, all_close(solution%s, sref, abs_tol=0.05_dp))                    ! Matching slack.
      if (allocated(error)) return
      call check(error, is_close(solution%objective_value, cost_ref, abs_tol=0.05_dp))   ! Matching cost.
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
      integer(ilp), parameter :: nleq = 2, ngeq = 0, neq = 0
      real(dp) :: xref(n), sref(m), cost_ref
      real(dp) :: c(n)
      real(dp) :: Aleq(nleq, n), bleq(nleq)
      type(dense_lp_type) :: problem
      type(lp_solution) :: solution

      !> Reference solution.
      cost_ref = 20.0_dp
      xref = [0.0_dp, 0.0_dp, 5.0_dp]
      sref = [5.0_dp, 0.0_dp]

      !> Cost vector.
      c = [2.0_dp, 3.0_dp, 4.0_dp]

      !> Linear <= inequalities.
      Aleq(1, :) = [3.0_dp, 2.0_dp, 1.0_dp]
      Aleq(2, :) = [2.0_dp, 5.0_dp, 3.0_dp]
      bleq = [10.0_dp, 15.0_dp]

      !> Construct the LP problem for LightConvex.
      problem = linear_program(c, Aleq=Aleq, bleq=bleq)

      !> Solve the problem.
      solution = solve(problem, alg=PrimalSimplex())

      call check(error, is_optimal(problem))                             ! Optimal solution found.
      if (allocated(error)) return
      call check(error, all_close(solution%x, xref))                     ! Matching primal.
      if (allocated(error)) return
      call check(error, all_close(solution%s, sref))                     ! Matching slack.
      if (allocated(error)) return
      call check(error, is_close(solution%objective_value, cost_ref))    ! Matching cost.
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
      integer(ilp), parameter :: nleq = 0, ngeq = 1, neq = 0
      real(dp) :: c(n), Ageq(ngeq, n), bgeq(ngeq)
      type(dense_lp_type) :: problem
      type(lp_solution) :: solution

      !> Cost vector.
      c = [1.0_dp, 1.0_dp]

      !> Linear >= inequalities.
      Ageq(1, :) = -1.0_dp; bgeq = 1.0_dp

      !> Construct the LP problem for LightConvex.
      problem = linear_program(c, Ageq=Ageq, bgeq=bgeq)

      !> Solve the problem.
      solution = solve(problem, alg=PrimalSimplex())

      call check(error,.not. is_feasible(problem))
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
      integer(ilp), parameter :: nleq = 0, ngeq = 1, neq = 0
      real(dp) :: c(n), Ageq(ngeq, n), bgeq(ngeq)
      type(dense_lp_type) :: problem
      type(lp_solution) :: solution

      !> Cost vector.
      c = [1.0_dp, 1.0_dp]

      !> Linear >= inequalities.
      Ageq(1, :) = 1.0_dp; bgeq = 5.0_dp

      !> Construct the LP problem for LightConvex.
      problem = linear_program(c, Ageq=Ageq, bgeq=bgeq)

      !> Solve the problem.
      solution = solve(problem, alg=PrimalSimplex())

      call check(error, is_unbounded(problem))
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
         integer(ilp), parameter :: nleq = 3, ngeq = 0, neq = 0
         real(dp) :: c(n), Aleq(nleq, n), bleq(nleq)
         real(dp) :: xref(n), cost_ref
         type(dense_lp_type) :: problem
         type(lp_solution) :: solution

         !> Reference solution.
         cost_ref = 6.5_dp
         xref = [1.0_dp, 1.0_dp, 0.5_dp, 0.0_dp]

         !> Linear cost function.
         c = [2.0_dp, 4.0_dp, 1.0_dp, 1.0_dp]

         !> Linear <= ineaquality constraints.
         Aleq(1, :) = [2.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]; bleq(1) = 3.0_dp
         Aleq(2, :) = [0.0_dp, 1.0_dp, 4.0_dp, 1.0_dp]; bleq(2) = 3.0_dp
         Aleq(3, :) = [1.0_dp, 3.0_dp, 0.0_dp, 1.0_dp]; bleq(3) = 4.0_dp

         !> Construct the problem for LightConvex.
         problem = linear_program(c, Aleq=Aleq, bleq=bleq)

         !> Solve the problem.
         solution = solve(problem, alg=PrimalSimplex())

         !> Check correcteness.
         call check(error, is_optimal(problem))                             ! Optimal solution found.
         if (allocated(error)) return
         call check(error, all_close(solution%x, xref))                     ! Matching primal.
         if (allocated(error)) return
         call check(error, is_close(solution%objective_value, cost_ref))    ! Matching cost.
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
         integer(ilp), parameter :: nleq = 3, ngeq = 0, neq = 0
         real(dp) :: c(n), Aleq(nleq, n), bleq(nleq)
         real(dp) :: xref(n), cost_ref
         type(dense_lp_type) :: problem
         type(lp_solution) :: solution

         !> Reference solution.
         cost_ref = 1.0/20.0_dp
         xref = [1.0_dp/25.0_dp, 0.0_dp, 1.0_dp, 0.0_dp]

         !> Linear cost function.
         c = [3.0_dp/4.0_dp, -150.0_dp, 1.0_dp/50.0_dp, -6.0_dp]

         !> Linear <= inequality constraints.
         Aleq(1, :) = [1.0_dp/4.0_dp, -60.0_dp, -1.0_dp/25.0_dp, 9.0_dp]; bleq(1) = 0.0_dp
         Aleq(2, :) = [0.5_dp, -90.0_dp, -1.0_dp/50.0_dp, 3.0_dp]; bleq(2) = 0.0_dp
         Aleq(3, :) = [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp]; bleq(3) = 1.0_dp

         !> Construct the problem for LightConvex.
         problem = linear_program(c, Aleq=Aleq, bleq=bleq)

         !> Solve the problem.
         solution = solve(problem, alg=PrimalSimplex())

         !> Check correcteness.
         call check(error, is_optimal(problem))                             ! Optimal solution found.
         if (allocated(error)) return
         call check(error, all_close(solution%x, xref))                     ! Matching primal.
         if (allocated(error)) return
         call check(error, is_close(solution%objective_value, cost_ref))    ! Matching cost.
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
         integer(ilp), parameter :: nleq = 0, ngeq = 0, neq = 2
         real(dp) :: xref(n), cost_ref
         real(dp) :: c(n), Aeq(neq, n), beq(neq)
         type(dense_lp_type) :: problem
         type(lp_solution) :: solution

         !> Reference solution.
         cost_ref = 60.0_dp/7.0_dp
         xref = [0.0_dp, 4.0_dp/7.0_dp, 12.0_dp/7.0_dp, 0.0_dp, 0.0_dp]

         !> Linear cost function.
         c = [-1.0_dp, -6.0_dp, 7.0_dp, -1.0_dp, -5.0_dp]

         !> Linear equality constraints.
         Aeq(1, :) = [5.0_dp, -4.0_dp, 13.0_dp, -2.0_dp, 1.0_dp]; beq(1) = 20.0_dp
         Aeq(2, :) = [1.0_dp, -1.0_dp, 5.0_dp, -1.0_dp, 1.0_dp]; beq(2) = 8.0_dp

         !> Construct the problem for LightConvex.
         problem = linear_program(c, Aeq=Aeq, beq=beq)

         !> Solve the problem.
         solution = solve(problem, alg=PrimalSimplex())

         !> Check correcteness.
         call check(error, is_optimal(problem))                             ! Optimal solution found.
         if (allocated(error)) return
         call check(error, all_close(solution%x, xref))                     ! Matching primal.
         if (allocated(error)) return
         call check(error, is_close(solution%objective_value, cost_ref))    ! Matching cost.
         if (allocated(error)) return
      end block

      ! The problem reads
      !
      !   maximize    2 x1 - 3 x2 + x3 + x4
      !   subject to  x1 + 2 x2 + x3 + x4 = 3
      !               x1 - 2 x2 + 2 x3 + x4 = -2
      !               3 x1 - x2 - x4 = -1
      !               x1, x2, x3, x4 >= 0
      !
      ! The solution is given by x = [3/16, 5/4, 0, 5/16] with
      ! optimal cost c = -49/16.
      block
         integer(ilp), parameter :: m = 3, n = 4, maxiter = 10
         integer(ilp), parameter :: nleq = 0, ngeq = 0, neq = 3
         real(dp) :: xref(n), cost_ref
         real(dp) :: c(n), Aeq(neq, n), beq(neq)
         type(dense_lp_type) :: problem
         type(lp_solution) :: solution

         !> Reference solution.
         cost_ref = -49.0_dp/16.0_dp
         xref = [3.0_dp/16.0_dp, 5.0_dp/4.0_dp, 0.0_dp, 5.0_dp/16.0_dp]

         !> Linear cost function.
         c = [2.0_dp, -3.0_dp, 1.0_dp, 1.0_dp]

         !> Linear equality constraints.
         Aeq(1, :) = [1.0_dp, 2.0_dp, 1.0_dp, 1.0_dp]; beq(1) = 3.0_dp
         Aeq(2, :) = [-1.0_dp, 2.0_dp, -2.0_dp, -1.0_dp]; beq(2) = 2.0_dp
         Aeq(3, :) = [-3.0_dp, 1.0_dp, 0.0_dp, 1.0_dp]; beq(3) = 1.0_dp

         !> Construct the problem for LightConvex.
         problem = linear_program(c, Aeq=Aeq, beq=beq)

         !> Solve the problem.
         solution = solve(problem, alg=PrimalSimplex())

         !> Check correcteness.
         call check(error, is_optimal(problem))                             ! Optimal solution found.
         if (allocated(error)) return
         call check(error, all_close(solution%x, xref))                     ! Matching primal.
         if (allocated(error)) return
         call check(error, is_close(solution%objective_value, cost_ref))    ! Matching cost.
         if (allocated(error)) return
      end block

      ! Same problem as above with the equalities replaced by inequalities <=.
      ! The solution is x = [0, 1, 0, 0] with optimal cost c = -3.
      block
         integer(ilp), parameter :: m = 3, n = 4, maxiter = 10
         integer(ilp), parameter :: nleq = 1, ngeq = 2, neq = 0
         real(dp) :: xref(n), cost_ref
         real(dp) :: c(n)
         real(dp) :: Aleq(nleq, n), bleq(nleq)
         real(dp) :: Ageq(ngeq, n), bgeq(ngeq)
         type(dense_lp_type) :: problem
         type(lp_solution) :: solution

         !> Reference solution.
         cost_ref = -3.0_dp
         xref = [0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]

         !> Linear cost function.
         c = [2.0_dp, -3.0_dp, 1.0_dp, 1.0_dp]

         !> Linear inequality constraints.
         Aleq(1, :) = [1.0_dp, 2.0_dp, 1.0_dp, 1.0_dp]; bleq(1) = 3.0_dp

         Ageq(1, :) = [-1.0_dp, 2.0_dp, -2.0_dp, -1.0_dp]; bgeq(1) = 2.0_dp
         Ageq(2, :) = [-3.0_dp, 1.0_dp, 0.0_dp, 1.0_dp]; bgeq(2) = 1.0_dp

         !> Construct the problem for LightConvex.
         problem = linear_program(c, Aleq=Aleq, bleq=bleq, Ageq=Ageq, bgeq=bgeq)

         !> Solve the problem.
         solution = solve(problem, alg=PrimalSimplex())

         !> Check correcteness.
         call check(error, is_optimal(problem))                             ! Optimal solution found.
         if (allocated(error)) return
         call check(error, all_close(solution%x, xref))                     ! Matching primal.
         if (allocated(error)) return
         call check(error, is_close(solution%objective_value, cost_ref))    ! Matching cost.
         if (allocated(error)) return
      end block
   end subroutine test_caltech_examples

   ! This example is taken lecture notes from Cornell.
   !    https://www.cs.cornell.edu/courses/cs6820/2023fa/handouts/splex.pdf
   !
   ! The problem reads
   !
   !   maximize    2 x1 + 3 x2
   !   subject to  x1 + x2 <= 8
   !               2 x1 + x2 <= 12
   !               x1 + 2 x2 <= 14
   !               x1, x2 >= 0
   !
   ! Solution is given by x = [2, 6] with cost c = 22 and
   ! slack variables s = [0, 2, 0].
   subroutine test_cornell_example(error)
      type(error_type), allocatable, intent(out) :: error
      integer(ilp), parameter :: m = 3, n = 2, maxiter = 10
      integer(ilp), parameter :: nleq = 3, ngeq = 0, neq = 0
      real(dp) :: xref(n), sref(m), cost_ref
      real(dp) :: c(n), Aleq(nleq, n), bleq(nleq)
      type(dense_lp_type) :: problem
      type(lp_solution) :: solution

      !> Reference solution.
      cost_ref = 22.0_dp; xref = [2.0_dp, 6.0_dp]
      sref = [0.0_dp, 2.0_dp, 0.0_dp]

      !> Linear cost.
      c = [2.0_dp, 3.0_dp]

      !> Linear <= inequalities.
      Aleq(1, :) = 1.0_dp; bleq(1) = 8.0_dp
      Aleq(2, :) = [2.0_dp, 1.0_dp]; bleq(2) = 12.0_dp
      Aleq(3, :) = [1.0_dp, 2.0_dp]; bleq(3) = 14.0_dp

      !> Construct LP problem for LightConvex.
      problem = linear_program(c, Aleq=Aleq, bleq=bleq)

      !> Solve the problem.
      solution = solve(problem, alg=PrimalSimplex())

      !> Check correcteness.
      call check(error, is_optimal(problem))                             ! Optimal solution found.
      if (allocated(error)) return
      call check(error, all_close(solution%x, xref))                     ! Matching primal.
      if (allocated(error)) return
      call check(error, all_close(solution%s, sref))                     ! Matching slack.
      if (allocated(error)) return
      call check(error, is_close(solution%objective_value, cost_ref))    ! Matching cost.
      if (allocated(error)) return
   end subroutine test_cornell_example
end module TestLinearPrograms
