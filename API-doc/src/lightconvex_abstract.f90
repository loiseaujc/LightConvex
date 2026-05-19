module lightconvex_abstract
   use stdlib_optval, only: optval
   use lightkrylov, only: abstract_vector_rdp, abstract_linop_rdp, abstract_sym_linop_rdp
   use lightconvex_constants, only: ilp, dp, lk, &
                                    optimal_status, infeasible_status, &
                                    unbounded_status, maxiter_exceeded
   implicit none(type, external)
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

   !> Base type for KKT vector.
   type, extends(abstract_vector_rdp), public :: kkt_vector
      class(abstract_vector_rdp), allocatable :: x
      !! Primal vector.
      class(abstract_vector_rdp), allocatable :: y
      !! Dual vector.
   contains
      procedure, pass(self), public :: zero => kkt_zero
      procedure, pass(self), public :: dot => kkt_dot
      procedure, pass(self), public :: scal => kkt_scal
      procedure, pass(self), public :: axpby => kkt_axpby
      procedure, pass(self), public :: rand => kkt_rand
      procedure, pass(self), public :: get_size => kkt_get_size
   end type kkt_vector

   !> Base type for Hessian matrices.
   type, abstract, extends(abstract_sym_linop_rdp), public :: abstract_hessian_rdp
      procedure(hinv_iface), pointer, nopass :: apply_hinv => null()
   contains
      procedure(hessian_matvec_iface), pass(self), deferred :: matvec
      procedure, pass(self) :: hinv_matvec
      procedure, pass(self) :: has_hinv
   end type abstract_hessian_rdp
   abstract interface
      subroutine hessian_matvec_iface(self, vec_in, vec_out)
         import :: abstract_hessian_rdp, abstract_vector_rdp
         implicit none(type, external)
         class(abstract_hessian_rdp), intent(inout) :: self
         class(abstract_vector_rdp), intent(in) :: vec_in
         class(abstract_vector_rdp), intent(out) :: vec_out
      end subroutine hessian_matvec_iface

      subroutine hinv_iface(self, vec_in, vec_out)
         import :: abstract_hessian_rdp, abstract_vector_rdp
         implicit none(type, external)
         class(abstract_hessian_rdp), intent(inout) :: self
         class(abstract_vector_rdp), intent(in) :: vec_in
         class(abstract_vector_rdp), intent(out) :: vec_out
      end subroutine hinv_iface
   end interface

   !> Base type for the KKT operator.
   type, extends(abstract_sym_linop_rdp), public :: kkt_linop
      class(abstract_hessian_rdp), allocatable :: H    ! Hessian operator.
      class(abstract_linop_rdp), allocatable :: A      ! Equality constraint operator.
      class(abstract_vector_rdp), allocatable :: wrk_x ! Working array.
   contains
      procedure, pass(self), public :: matvec => kkt_matvec
   end type kkt_linop

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

   !-----------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURE FOR THE KKT VECTOR     -----
   !-----------------------------------------------------------

   subroutine kkt_zero(self)
      implicit none(type, external)
      class(kkt_vector), intent(inout) :: self
      call self%x%zero() ! Primal vector.
      call self%y%zero() ! Dual vector.
   end subroutine kkt_zero

   real(dp) function kkt_dot(self, vec) result(alpha)
      implicit none(type, external)
      class(kkt_vector), intent(in) :: self
      class(abstract_vector_rdp), intent(in) :: vec
      select type (vec)
      type is (kkt_vector)
         alpha = self%x%dot(vec%x) + self%y%dot(vec%y)
      class default
         error stop "Type error in KKT_DOT."
      end select
   end function kkt_dot

   subroutine kkt_scal(self, alpha)
      implicit none(type, external)
      class(kkt_vector), intent(inout) :: self
      real(dp), intent(in) :: alpha
      call self%x%scal(alpha) ! Primal
      call self%y%scal(alpha) ! Dual
   end subroutine kkt_scal

   subroutine kkt_axpby(alpha, vec, beta, self)
      implicit none(type, external)
      class(kkt_vector), intent(inout) :: self
      class(abstract_vector_rdp), intent(in) :: vec
      real(dp), intent(in) :: alpha, beta
      select type (vec)
      type is (kkt_vector)
         call self%x%axpby(alpha, vec%x, beta)
         call self%y%axpby(alpha, vec%y, beta)
      class default
         error stop "Type error in KKT_AXPBY."
      end select
   end subroutine kkt_axpby

   integer(ilp) function kkt_get_size(self) result(n)
      class(kkt_vector), intent(in) :: self
      n = self%x%get_size() + self%y%get_size()
   end function kkt_get_size

   subroutine kkt_rand(self, ifnorm)
      implicit none(type, external)
      class(kkt_vector), intent(inout) :: self
      logical(lk), optional, intent(in) :: ifnorm
      logical(lk) :: normalize
      real(dp), parameter :: alpha = 1.0_dp/sqrt(2.0_dp)
      normalize = optval(ifnorm, .true.)
      ! Generate random vectors.
      call self%x%rand(ifnorm)
      call self%y%rand(ifnorm)
      ! Normalize vector.
      if (normalize) then
         call self%x%scal(alpha)
         call self%y%scal(alpha)
      end if
   end subroutine kkt_rand

   !----------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR THE HESSIAN MATRIX     -----
   !----------------------------------------------------------------

   pure logical(lk) function has_hinv(self)
      class(abstract_hessian_rdp), intent(in) :: self
      has_hinv = associated(self%apply_hinv)
   end function has_hinv

   subroutine hinv_matvec(self, vec_in, vec_out)
      class(abstract_hessian_rdp), intent(inout) :: self
      class(abstract_vector_rdp), intent(in) :: vec_in
      class(abstract_vector_rdp), intent(out) :: vec_out

      if (self%has_hinv()) then
         call self%apply_hinv(self, vec_in, vec_out)
      else
         error stop "Hessian inverse is not available."
      end if
   end subroutine hinv_matvec

   !-----------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURE FOR THE KKT MATRIX     -----
   !-----------------------------------------------------------

   subroutine kkt_matvec(self, vec_in, vec_out)
      class(kkt_linop), intent(inout) :: self
      class(abstract_vector_rdp), intent(in) :: vec_in
      class(abstract_vector_rdp), intent(out) :: vec_out
      select type (vec_in)
      type is (kkt_vector)
         select type (vec_out)
         type is (kkt_vector)
            ! Working vector.
            if (.not. allocated(self%wrk_x)) allocate (self%wrk_x, mold=vec_in%x)
            call self%wrk_x%zero()

            ! Upper blocks : x' = H @ x + A.T @ y.
            call vec_out%zero()
            call self%H%matvec(vec_in%x, vec_out%x) ! x' = H @ x
            call self%A%rmatvec(vec_in%y, self%wrk_x)      ! z = A.T @ y
            call vec_out%x%add(self%wrk_x)

            ! Lower blocks : y' = A @ x.
            call self%A%matvec(vec_in%x, vec_out%y)
         class default
            error stop "Type error in KKT_MATVEC."
         end select
      class default
         error stop "Type error in KKT_MATVEC."
      end select
   end subroutine kkt_matvec

end module lightconvex_abstract
