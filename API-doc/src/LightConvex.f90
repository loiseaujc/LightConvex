module LightConvex
   use lightconvex_abstract, only: is_optimal, is_feasible, is_unbounded
   use lightconvex_lp, only: linear_program, dense_lp_type, lp_solution, &
                             Dantzig, auxiliary_function, &
                             solve, PrimalSimplex
   implicit none(external)
   public
end module LightConvex
