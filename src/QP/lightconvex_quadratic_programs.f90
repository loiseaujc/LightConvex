module lightconvex_qp
   use stdlib_optval, only: optval
   use lightconvex_constants, only: dp, ilp, lk
   use lightconvex_abstract, only: abstract_cvx_problem, abstract_cvx_solution, abstract_cvx_solver, &
                                   abstract_cvx_vector, AbstractMatrix, AbstractSymMatrix
   use lightconvex_kkt, only: abstract_kkt_solver, kkt_info, is_successful
   implicit none(type, external)
   private

   !> Base type for QP problem.
   type, extends(abstract_cvx_problem), public :: qp_problem
      class(AbstractSymMatrix), allocatable :: P
      class(AbstractMatrix), allocatable :: A
      class(abstract_cvx_vector), allocatable :: q, l, u
   end type qp_problem

   !> Base type for solution of QP.
   type, extends(abstract_cvx_solution), public :: qp_solution
      class(abstract_cvx_vector), allocatable :: x
      class(abstract_cvx_vector), allocatable :: y
      class(abstract_cvx_vector), allocatable :: z
      integer(ilp) :: iterations = 0_ilp
      real(dp) :: prim_res = huge(1.0_dp), dual_res = huge(1.0_dp)
   end type qp_solution

   !-----------------------------------------------------------
   !-----   ADMM SOLVER FOR CONVEX QUADRATIC PROGRAMS     -----
   !-----------------------------------------------------------

   type, extends(abstract_cvx_solver), public :: ADMM
      real(dp) :: rho = 0.1_dp
      real(dp) :: rho_eq_scale = 1.0e3_dp
      real(dp) :: eq_tol = 1.0e-9_dp
      real(dp) :: sigma = 1.0e-6_dp
      real(dp) :: alpha = 1.6_dp
      real(dp) :: atol = 1.0e-3_dp
      real(dp) :: rtol = 1.0e-3_dp
      integer(ilp) :: maxiter = 4000_ilp
      integer(ilp) :: check_every = 25_ilp
   contains
   end type ADMM

   interface
      module subroutine admm_solver(settings, P, q, A, l, u, kkt, x, y, z, prim_res, dual_res)
         implicit none(type, external)
         type(ADMM), intent(in) :: settings
         class(AbstractSymMatrix), intent(in) :: P
         class(AbstractMatrix), intent(in) :: A
         class(abstract_cvx_vector), intent(in) :: q, l, u
         class(abstract_kkt_solver), intent(inout) :: kkt
         class(abstract_cvx_vector), intent(inout) :: x, y, z
         real(dp), allocatable, intent(out) :: prim_res(:), dual_res(:)
      end subroutine admm_solver
   end interface
   public :: admm_solver
contains
end module lightconvex_qp
