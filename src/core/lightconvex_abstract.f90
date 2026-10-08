module lightconvex_abstract
   use LightKrylov, only: abstract_vector_rdp, abstract_linop_rdp, abstract_sym_linop_rdp
   use lightconvex_constants, only: ilp, dp, lk, &
                                    optimal_status, infeasible_status, &
                                    unsolved_status, unbounded_status, maxiter_exceeded
   implicit none(type, external)
   private

   public :: abstract_vector_rdp

   !> Base type for defining convex problems.
   type, abstract, public :: abstract_cvx_problem
      private
      integer(ilp) :: status = unsolved_status
   contains
      procedure, pass(self), public :: set_status
   end type abstract_cvx_problem

   !> Base type for defining convex solvers.
   type, abstract, public :: abstract_cvx_solver
   end type abstract_cvx_solver

   !> Base type for defining solutions to convex problems.
   type, abstract, public :: abstract_cvx_solution
   end type abstract_cvx_solution

   !------------------------------------
   !-----     ABSTRACT VECTORS     -----
   !------------------------------------

   type, abstract, public, extends(abstract_vector_rdp) :: abstract_cvx_vector
   contains
      procedure(norm_inf_iface), deferred :: norm_inf
      procedure(hadamard_iface), deferred :: hadamard
      procedure(reciprocal_iface), deferred :: reciprocal
      procedure(fill_iface), deferred :: fill
      procedure(clip_iface), deferred :: clip
      procedure(mask_le_iface), deferred :: mask_le
   end type abstract_cvx_vector

   abstract interface
      function norm_inf_iface(self) result(val)
         import abstract_cvx_vector, dp
         implicit none(type, external)
         class(abstract_cvx_vector), intent(in) :: self
         real(dp) :: val
      end function norm_inf_iface

      subroutine hadamard_iface(self, d)
         import abstract_cvx_vector
         implicit none(type, external)
         class(abstract_cvx_vector), intent(inout) :: self
         class(abstract_cvx_vector), intent(in) :: d
      end subroutine hadamard_iface

      subroutine reciprocal_iface(self)
         import abstract_cvx_vector
         implicit none(type, external)
         class(abstract_cvx_vector), intent(inout) :: self
      end subroutine reciprocal_iface

      subroutine fill_iface(self, val)
         import abstract_cvx_vector, dp
         implicit none(type, external)
         class(abstract_cvx_vector), intent(inout) :: self
         real(dp), intent(in) :: val
      end subroutine fill_iface

      subroutine clip_iface(self, l, u)
         import abstract_cvx_vector, dp
         implicit none(type, external)
         class(abstract_cvx_vector), intent(inout) :: self
         class(abstract_cvx_vector), intent(in) :: l, u
      end subroutine clip_iface

      subroutine mask_le_iface(self, val)
         import abstract_cvx_vector, dp
         implicit none(type, external)
         class(abstract_cvx_vector), intent(inout) :: self
         real(dp), intent(in) :: val
      end subroutine mask_le_iface
   end interface

   !---------------------------------------------
   !-----     ABSTRACT LINEAR OPERATORS     -----
   !---------------------------------------------

   !> Base type for general dense matrices.
   type, abstract, public :: AbstractMatrix
   contains
      procedure(matvec_iface), pass(self), deferred :: matvec
   end type AbstractMatrix

   abstract interface
      subroutine matvec_iface(self, alpha, x, beta, y, op)
         import AbstractMatrix, abstract_cvx_vector, dp
         implicit none(type, external)
         class(AbstractMatrix), intent(in) :: self
         real(dp), intent(in) :: alpha, beta
         class(abstract_cvx_vector), intent(in) :: x
         class(abstract_cvx_vector), intent(inout) :: y
         character(1), intent(in) :: op
      end subroutine matvec_iface
   end interface

   !> Base type for symmetric dense matrices.
   type, abstract, public :: AbstractSymMatrix
   contains
      procedure(symmatvec_iface), pass(self), deferred :: matvec
   end type AbstractSymMatrix

   abstract interface
      subroutine symmatvec_iface(self, alpha, x, beta, y)
         import AbstractSymMatrix, abstract_cvx_vector, dp
         implicit none(type, external)
         class(AbstractSymMatrix), intent(in) :: self
         real(dp), intent(in) :: alpha, beta
         class(abstract_cvx_vector), intent(in) :: x
         class(abstract_cvx_vector), intent(inout) :: y
      end subroutine symmatvec_iface
   end interface

   !-------------------------------------
   !-----     UTILITY FUNCTIONS     -----
   !-------------------------------------

   interface
      pure module subroutine set_status(self, status)
         implicit none(type, external)
         class(abstract_cvx_problem), intent(inout) :: self
         integer(ilp), intent(in) :: status
      end subroutine set_status

      pure logical(lk) module function is_optimal(problem) result(bool)
         implicit none(type, external)
         class(abstract_cvx_problem), intent(in) :: problem
      end function is_optimal

      pure logical(lk) module function is_feasible(problem) result(bool)
         implicit none(type, external)
         class(abstract_cvx_problem), intent(in) :: problem
      end function is_feasible

      pure logical(lk) module function is_unbounded(problem) result(bool)
         implicit none(type, external)
         class(abstract_cvx_problem), intent(in) :: problem
      end function is_unbounded

      pure logical(lk) module function is_solved(problem) result(bool)
         implicit none(type, external)
         class(abstract_cvx_problem), intent(in) :: problem
      end function is_solved
   end interface
   public :: is_optimal, is_feasible, is_unbounded, is_solved

contains
   module procedure set_status
   self%status = status
   end procedure set_status

   module procedure is_solved
   bool = problem%status /= unsolved_status
   end procedure is_solved

   module procedure is_optimal
   bool = problem%status == optimal_status
   end procedure is_optimal

   module procedure is_feasible
   bool = (problem%status == optimal_status) .and. (problem%status == unbounded_status)
   end procedure is_feasible

   module procedure is_unbounded
   bool = problem%status == unbounded_status
   end procedure is_unbounded
end module lightconvex_abstract
