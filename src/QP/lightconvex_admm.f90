submodule(lightconvex_qp) lightconvex_admm
   implicit none(type, external)
contains

   ! See https://web.stanford.edu/~boyd/papers/pdf/osqp.pdf#page=4.73
   ! for reference.
   module procedure admm_solver
   !> Primal space.
   class(abstract_cvx_vector), allocatable :: xt, rhs_x
   class(abstract_cvx_vector), allocatable :: d1, Px, Aty, tmp_n
   !> Dual space.
   class(abstract_cvx_vector), allocatable :: zt, rho
   class(abstract_cvx_vector), allocatable :: rhs_y, nu
   class(abstract_cvx_vector), allocatable :: d2, Ax, tmp_m
   !> KKT solver.
   type(kkt_info) :: kkt_status
   !> Residual.
   real(dp) :: r_prim, r_dual, eps_prim, eps_dual
   real(dp) :: px_norm_inf, ax_norm_inf, aty_norm_inf, q_norm_inf, z_norm_inf
   logical(lk) :: converged
   !> Miscellaneous.
   integer(ilp) :: i, check_counter

   !----------------------------------
   !-----     INITIALIZATION     -----
   !----------------------------------

   !> Primal-related variables.
   allocate (xt, rhs_x, d1, Px, Aty, tmp_n, source=q)
   !> Dual-related variables.
   allocate (nu, zt, rhs_y, d2, Ax, rho, tmp_m, source=l)
   !> Residual history.
   check_counter = 1
   block
      integer(ilp) :: npts
      npts = settings%maxiter/settings%check_every
      allocate (prim_res(npts), dual_res(npts), source=0.0_dp)
   end block

   !-------------------------
   !-----     SETUP     -----
   !-------------------------

   q_norm_inf = q%norm_inf()

   !> Penalty vectors.
   tmp_m = u; call tmp_m%sub(l); call tmp_m%mask_le(settings%eq_tol)
   call rho%fill(settings%rho); call rho%axpby((settings%rho_eq_scale - 1.0_dp)*settings%rho, tmp_m, 1.0_dp)
   d2 = rho; call d2%reciprocal()
   call d1%fill(settings%sigma)

   !> Initialize KKT solver.
   call kkt%update(d1, d2, kkt_status)
   if (.not. is_successful(kkt_status)) error stop "admm_solver: KKT failed to initialize."

   call A%matvec(1.0_dp, x, 0.0_dp, z, "n")
   call z%clip(l, u)

   !-----------------------------
   !-----     MAIN LOOP     -----
   !-----------------------------

   converged = .false.
   admm_loop: do i = 1, settings%maxiter

      !> Right-hand side vectors.
      call rhs_x%axpby(settings%sigma, x, 0.0_dp)
      call rhs_x%sub(q)

      tmp_m = y; call tmp_m%hadamard(d2)
      rhs_y = z; call rhs_y%sub(tmp_m)

      !> KKT solver.
      call kkt%solve(rhs_x, rhs_y, xt, nu, kkt_status)
      if (.not. is_successful(kkt_status)) error stop "admm_solver: kkt%solve failed."

      !> Link variable.
      tmp_m = nu; call tmp_m%sub(y); call tmp_m%hadamard(d2)
      zt = z; call zt%add(tmp_m)

      !> Relaxation.
      call x%axpby(settings%alpha, xt, 1 - settings%alpha)
      call zt%axpby(1 - settings%alpha, z, settings%alpha)

      !> Projection (prox of the indicator of [l, u]).
      call z%axpby(1.0_dp, d2, 0.0_dp); call z%hadamard(y)
      call z%add(zt); call z%clip(l, u)

      !> Dual update.
      tmp_m = zt; call tmp_m%sub(z)
      call tmp_m%hadamard(rho); call y%add(tmp_m)

      !> Termination test.
      if (mod(i, settings%check_every) == 0 .or. i == settings%maxiter) then
         call A%matvec(1.0_dp, x, 0.0_dp, Ax, op="n"); ax_norm_inf = Ax%norm_inf()
         call Ax%sub(z)
         r_prim = Ax%norm_inf()

         call P%matvec(1.0_dp, x, 0.0_dp, Px); px_norm_inf = Px%norm_inf()
         call A%matvec(1.0_dp, y, 0.0_dp, Aty, op="t"); aty_norm_inf = Aty%norm_inf()
         call Px%add(q); call Px%add(Aty)
         r_dual = Px%norm_inf()

         z_norm_inf = z%norm_inf()
         eps_prim = settings%atol + settings%rtol*max(ax_norm_inf, z_norm_inf, tiny(1.0_dp))
         eps_dual = settings%atol + settings%rtol*max(px_norm_inf, aty_norm_inf, q_norm_inf, tiny(1.0_dp))

         prim_res(check_counter) = eps_prim
         dual_res(check_counter) = eps_dual
         check_counter = check_counter + 1

         if ((r_prim <= eps_prim) .and. (r_dual <= eps_dual)) then
            converged = .true.
            exit admm_loop
         end if

         ! TODO: rho adaptation hook (update d2, call kkt%update again)
         ! TODO: infeasibility certificates (primal / dual)
      end if
   end do admm_loop

   !> Truncate residual arrays if early stopping.
   prim_res = prim_res(:check_counter - 1)
   dual_res = dual_res(:check_counter - 1)
   end procedure admm_solver

end submodule lightconvex_admm
