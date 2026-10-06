module LightConvex
   use lightconvex_abstract, only: is_optimal, is_feasible, is_unbounded, &
                                   abstract_cvx_problem, abstract_vector_rdp
   use lightconvex_dense_vectors, only: dense_vector
   use lightconvex_lp, only: linear_program, dense_lp_type, lp_solution, &
                             Dantzig, auxiliary_function, &
                             solve, PrimalSimplex
   implicit none(type, external)
   public
end module LightConvex
