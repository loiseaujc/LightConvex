module lightconvex_abstract
   use lightconvex_constants, only: ilp, dp, lk, &
                                    optimal_status, infeasible_status, &
                                    unbounded_status, maxiter_exceeded
   implicit none(external)
   private

   !> Base type for defining convex problems.
   type, abstract, public :: abstract_cvx_problem
      private
      integer(ilp) :: status
   contains
      procedure, pass(self), public :: set_status
   end type abstract_cvx_problem

   !> Base type for defining convex solvers.
   type, abstract, public :: abstract_cvx_solver
   end type abstract_cvx_solver

   !> Base type for defining solutions to convex problems.
   type, abstract, public :: abstract_cvx_solution
   end type abstract_cvx_solution

   interface
      pure module subroutine set_status(self, status)
         implicit none(external)
         class(abstract_cvx_problem), intent(inout) :: self
         integer(ilp), intent(in) :: status
      end subroutine set_status

      pure logical(lk) module function is_optimal(problem) result(bool)
         implicit none(external)
         class(abstract_cvx_problem), intent(in) :: problem
      end function is_optimal

      pure logical(lk) module function is_feasible(problem) result(bool)
         implicit none(external)
         class(abstract_cvx_problem), intent(in) :: problem
      end function is_feasible

      pure logical(lk) module function is_unbounded(problem) result(bool)
         implicit none(external)
         class(abstract_cvx_problem), intent(in) :: problem
      end function is_unbounded
   end interface
   public :: is_optimal, is_feasible, is_unbounded

contains
   module procedure set_status
   self%status = status
   end procedure set_status

   module procedure is_optimal
   bool = problem%status == optimal_status
   end procedure is_optimal

   module procedure is_feasible
   bool = problem%status /= infeasible_status
   end procedure is_feasible

   module procedure is_unbounded
   bool = problem%status == unbounded_status
   end procedure is_unbounded
end module lightconvex_abstract
