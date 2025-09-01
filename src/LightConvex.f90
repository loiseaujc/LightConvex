module LightConvex
   use stdlib_linalg_constants, only: ilp, dp, lk
   use assert_m, only: assert
   implicit none(external)
   private

   public :: simplex

   !-----------------------------
   !-----     Constants     -----
   !-----------------------------

   real(dp), parameter :: eps = epsilon(1.0_dp)

   !----------------------------------------------
   !-----     FUNDAMENTAL ABSTRACT TYPES     -----
   !----------------------------------------------

   type, abstract :: AbstractConvexProblem
   end type AbstractConvexProblem

   type, abstract :: AbstractAlgorithm
   end type AbstractAlgorithm

   !--------------------------------------
   !-----     LINEAR PROGRAMMING     -----
   !--------------------------------------

   type, extends(AbstractConvexProblem), public :: lp_type
      private
      real(dp), allocatable :: c(:)
        !! Cost function.
      real(dp), allocatable :: A(:, :)
        !! Linear constraints.
      integer(ilp) :: nleq, ngeq, neq
        !! Number of constraints of each type.
   end type lp_type

   interface simplex
      pure module subroutine dense_simplex(A, nleq, ngeq, neq, icase, izrov, iposv)
         implicit none(external)
         real(dp), intent(inout) :: A(:, :)
        !! Simplex tableau of dimension n+2 x m
         integer(ilp), intent(in) :: nleq, ngeq, neq
        !! Number of constraints of each type.
         integer(ilp), intent(out) :: icase
        !! Return flag:
        !!  - icase = 1  : Optimal solution has been found.
        !!  - icase = 0  : Objective function is unbounded.
        !!  - icase = -1 : Problem is infeasible.
         integer(ilp), intent(out) :: iposv(:), izrov(:)
        !! Book-keeping for the primal and slack variables.
      end subroutine dense_simplex
   end interface simplex

   public :: say_hello
contains
   subroutine say_hello
      print *, "Hello, LightConvex!"
   end subroutine say_hello
end module LightConvex
