module lightconvex_constants
   use stdlib_linalg_constants, only: ilp, dp, lk
   implicit none(external)
   public

   real(dp), parameter :: eps = epsilon(1.0_dp)
    !! Machine precision
   real(dp), parameter :: tol = sqrt(eps)
    !! Tolerance used inside the different solvers.
   integer(ilp), parameter :: unbounded_status = 1
    !! Return flag for an unbounded problem.
   integer(ilp), parameter :: optimal_status = 0
    !! Return flag for when an optimal solution has been computed.
   integer(ilp), parameter :: infeasible_status = -1
    !! Return flag for an infeasible problem.
   integer(ilp), parameter :: maxiter_exceeded = -2
    !! Return flag for excessive number of iterations.
end module lightconvex_constants
