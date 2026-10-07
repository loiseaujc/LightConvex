submodule(lightconvex_kkt) lightconvex_dense_kkt
   use stdlib_linalg_lapack, only: sytrf, sytrs, symv, gemv
   use stdlib_linalg, only: norm
   use lightconvex_dense_vectors, only: dense_vector
   implicit none(type, external)
contains
   !------------------------------------
   !-----     DENSE KKT SOLVER     -----
   !------------------------------------

   module procedure create_dense_kkt_solver
   !> Allocate solver.
   allocate (solver)
   associate (n => size(P, 1), m => size(A, 1))
      !> Sanity checks.
      if (size(P, 2) /= n) error stop "dense_kkt_solver: P needs to be a square matrix."
      if (size(A, 2) /= n) error stop "dense_kkt_solver: A needs to have the same number of columns as P."
      solver%n = n; solver%m = m

      solver%max_refine = optval(max_refine, solver%max_refine)
      if (solver%max_refine < 0) error stop "dense_kkt_solver: max_refine needs to be non-negative."

      solver%tol = optval(tol, solver%tol)
      if (solver%tol < 0.0_dp) error stop "dense_kkt_solver: tol needs to be positive."

      ! --------------------

      !> Allocate matrices.
      allocate (solver%P, source=P)
      allocate (solver%A, source=A)
      allocate (solver%K(n + m, n + m), source=0.0_dp)
      allocate (solver%z(n + m, 1), solver%r(n + m), solver%dz(n + m, 1), source=0.0_dp)
      allocate (solver%d1(n), solver%d2(m), source=0.0_dp)
      allocate (solver%e1(n), solver%e2(m), source=0.0_dp)
      allocate (solver%ipiv(n + m), source=0_ilp)

      !> Workspace query.
      block
         character(1), parameter :: uplo = "L"
         integer(ilp), parameter :: lwork = -1_ilp
         integer(ilp) :: info
         real(dp) :: dummy_work(1)
         !> Workspace query.
         call sytrf(uplo, n + m, solver%K, n + m, solver%ipiv, dummy_work, lwork, info)
         if (info /= 0) error stop "dense_kkt_solver: error in sytrf."
         allocate (solver%workspace(int(dummy_work(1), kind=ilp)), source=0.0_dp)
      end block
   end associate
   end procedure create_dense_kkt_solver

   module procedure dense_update
   self%factorized = .false.
   select type (d1)
   type is (dense_vector)
      select type (d2)
      type is (dense_vector)
         associate (n => d1%get_size(), m => d2%get_size(), uplo => "L")
            !> Sanity check.
            if (n /= self%n) error stop "kkt%update: d1 has inconsistent dimensions."
            if (m /= self%m) error stop "kkt%update: d2 has inconsistent dimensions."

            !> Regularized diagonals.
            self%e1 = d1%data + optval(reg1, 0.0_dp)
            self%e2 = d2%data + optval(reg2, 0.0_dp)

            if (any(self%e1 < 0.0_dp) .or. any(self%e2 < 0.0_dp)) then
               info%status = kkt_invalid_regularization
               self%factorized = .false.
               return
            end if

            !> Assembling the regularized KKT matrix.
            call assemble_kkt_matrix(self%K, self%P, self%A, self%e1, self%e2)

            !> LDLT factorization of K.
            block
               integer(ilp) :: lwork, lapack_info
               self%ipiv = 0_ilp
               lwork = size(self%workspace, kind=ilp)
               call sytrf(uplo, n + m, self%K, n + m, self%ipiv, self%workspace, lwork, lapack_info)
               if (lapack_info < 0) error stop "kkt%update: error in sytrf."
               if (lapack_info > 0) error stop "kkt%update: D matrix in LDLT is singular."
            end block

            !> Store unregularized diagonal for iterative refinement.
            self%d1 = d1%data
            self%d2 = d2%data

            !> Update solver status.
            info%status = kkt_success
            self%factorized = .true.
         end associate
      class default
         error stop "kkt%update: Type of d2 should be dense_vector."
      end select
   class default
      error stop "kkt%update: Type of d1 should be dense_vector."
   end select
   end procedure dense_update

   module procedure dense_solve
   integer(ilp) :: i, lapack_info
   real(dp) :: rhs_norm, res, res_prev
   logical(lk) :: converged
   if (.not. self%factorized) then
      error stop "kkt%solve: KKT matrix has not been factorized. Aborting."
   end if
   select type (rhs_x)
   type is (dense_vector)
      select type (rhs_y)
      type is (dense_vector)
         select type (sol_x)
         type is (dense_vector)
            select type (sol_y)
            type is (dense_vector)

               associate (n => self%n, m => self%m, uplo => "L")
                  if (rhs_x%get_size() /= n) error stop "kkt%solve: rhs_x has inconsistent dimensions."
                  if (sol_x%get_size() /= n) error stop "kkt%solve: sol_x has inconsistent dimensions."
                  if (rhs_y%get_size() /= m) error stop "kkt%solve: rhs_y has inconsistent dimensions."
                  if (sol_y%get_size() /= m) error stop "kkt%solve: sol_y has inconsistent dimensions."

                  converged = .false.

                  !> Working vector.
                  self%z(:n, 1) = rhs_x%data; self%z(n + 1:, 1) = rhs_y%data

                  !> Solve the linear system.
                  call sytrs(uplo, n + m, 1, self%K, n + m, self%ipiv, self%z, n + m, lapack_info)
                  if (lapack_info /= 0) error stop "kkt%solve: Error in sytrs."

                  !> Compute residual.
                  call compute_residual(self%P, self%A, self%K, self%ipiv, self%d1, self%d2, &
                                        rhs_x%data, rhs_y%data, self%z, self%r)

                  rhs_norm = max(norm(rhs_x%data, "inf"), norm(rhs_y%data, "inf"), tiny(1.0_dp))
                  info%residual = norm(self%r, "inf")/rhs_norm

                  if (info%residual <= self%tol) then
                     converged = .true.
                  else
                     !> Store previous residual for comparison.
                     res_prev = info%residual

                     iterative_refinement: do i = 1, self%max_refine
                        !> Correction: K_reg dz = r, then z = z + dz.
                        self%dz(:, 1) = self%r
                        call sytrs(uplo, n + m, 1, self%K, n + m, self%ipiv, self%dz, n + m, lapack_info)
                        if (lapack_info /= 0) error stop "kkt%solve: Error sytrs (refinement)."
                        self%z = self%z + self%dz

                        !> r = rhs - K @ z (with K the unregularized matrix).
                        call compute_residual(self%P, self%A, self%K, self%ipiv, self%d1, self%d2, &
                                              rhs_x%data, rhs_y%data, self%z, self%r)

                        !> Residual norm.
                        info%residual = norm(self%r, "inf")/rhs_norm

                        if (.not. (info%residual <= huge(1.0_dp))) then
                           ! NaN or Inf: unstable factorization.
                           info%status = kkt_numerical_error
                           return
                        else if (info%residual <= self%tol) then
                           ! Solver converged.
                           converged = .true.
                           exit iterative_refinement
                        else if (info%residual > 0.5_dp*res_prev) then
                           ! Stagnation.
                           exit iterative_refinement
                        end if

                        !> Book-keeping.
                        res_prev = info%residual
                     end do iterative_refinement
                  end if

                  !> Approximated solutions.
                  info%status = merge(kkt_success, kkt_not_converged, converged)
                  sol_x%data = self%z(:n, 1)
                  sol_y%data = self%z(n + 1:, 1)
               end associate

            class default
               error stop "kkt%solve: Type of sol_y should be dense vector."
            end select
         class default
            error stop "kkt%solve: Type of sol_x should be dense vector."
         end select
      class default
         error stop "kkt%solve: Type of rhs_y should be dense_vector."
      end select
   class default
      error stop "kkt%solve: Type of rhs_x should be dense_vector."
   end select
   end procedure dense_solve

   pure subroutine assemble_kkt_matrix(K, P, A, e1, e2)
      real(dp), intent(out) :: K(:, :)
      real(dp), intent(in) :: P(:, :), A(:, :)
      real(dp), intent(in) :: e1(:), e2(:)
      integer(ilp) :: i, j
      associate (n => size(P, 1), m => size(A, 1))
         K = 0.0_dp
         do concurrent(j=1:n)
            !> Fill the (1, 1) block of K.
            do concurrent(i=j:n)
               K(i, j) = P(i, j) ! Store only the lower triangular part.
            end do
            K(j, j) = K(j, j) + e1(j)   ! Regularized diagonal.

            !> Fill the (2, 1) block of K.
            do concurrent(i=1:m)
               K(n + i, j) = A(i, j)
            end do
         end do

         !> Fill the (2, 2) diagonal block of K.
         do concurrent(j=1:m)
            K(n + j, n + j) = -e2(j)
         end do
      end associate
   end subroutine assemble_kkt_matrix

   pure subroutine compute_residual(P, A, K, ipiv, d1, d2, rhs_x, rhs_y, z, r)
      real(dp), intent(in) :: P(:, :), A(:, :), K(:, :), d1(:), d2(:)
      real(dp), intent(in) :: rhs_x(:), rhs_y(:), z(:, :)
      real(dp), intent(out) :: r(:)
      integer(ilp), intent(in) :: ipiv(:)
      integer(ilp) :: i, lda
      associate (n => size(P, 1), m => size(A, 1), uplo => "L")
         lda = max(1_ilp, m)
         !> r = rhs - K @ z (with K the unregularized matrix).
         !  -------------------------------------------------
         r(:n) = rhs_x; r(n + 1:) = rhs_y
         ! r_x = rhs_x - P @ x (P symmetric, lower triangle storage).
         call symv(uplo, n, -1.0_dp, P, n, z(:n, 1), 1, 1.0_dp, r(:n), 1)
         ! r_x = r_x - d1 .* x - A.T @ y.
         do concurrent(i=1:n)
            r(i) = r(i) - d1(i)*z(i, 1)
         end do
         call gemv("T", m, n, -1.0_dp, A, lda, z(n + 1:, 1), 1, 1.0_dp, r(:n), 1)
         ! r_y = rhs_y - A @ x + d2 .* y
         call gemv("N", m, n, -1.0_dp, A, lda, z(:n, 1), 1, 1.0_dp, r(n + 1:), 1)
         do concurrent(i=1:m)
            r(n + i) = r(n + i) - d2(i)*z(n + i, 1)
         end do
      end associate
   end subroutine compute_residual

end submodule lightconvex_dense_kkt
