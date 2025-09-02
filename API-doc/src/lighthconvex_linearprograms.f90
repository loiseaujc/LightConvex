submodule(lightconvex) lightconvex_lp
   use stdlib_math, only: arange, swap
   use stdlib_linalg, only: outer_product
   use stdlib_intrinsics, only: sum => stdlib_sum
   implicit none(external)

contains

   module procedure dense_standard_simplex
   integer(ilp) :: m, n
    !! Dimensions of the problem.
    !!  - m : Number of constraints (excluding the non-negativity).
    !!  - n : Number of variables (excluding the slack ones).
   integer(ilp) :: niter
    !! Current iteration number
   integer(ilp), allocatable           :: izrov(:)
    !! Book-keeping for variables being zeroed-out.
   integer(ilp), dimension(size(A, 2)) :: admissible_columns
    !! List of columns admissible for exchange.
   logical(lk), dimension(ngeq)        :: in_basis
    !! List of >= constraints currently in the basis.

   ! Miscellaneous
   integer(ilp) :: nl1, ip, kp, i, k, kh
   logical(lk)  :: init
   real(dp)     :: bmax

   !> Sanity checks.
   m = size(A, 1) - 2; n = size(A, 2) - 1
   call assert(assertion=m == nleq + ngeq + neq, &
               description="Tableau size is inconsistent with the number of constraints.")
   call assert(assertion=size(iposv) == m, &
               description="Dimension of iposv is inconsistant with the number of constraints.")
   call assert(assertion=all(A(2:, 1) >= 0.0_dp), &
               description="Constants b_i need to be non-negative.")

   !----- Initialization -----

   ! Index list of columns admissible for exchange.
   nl1 = n; admissible_columns(:n) = arange(n)
   izrov = admissible_columns(:n)   ! All variables are initially right-hand.
   iposv = n + arange(m)            ! Initial left-hand variables. <= constraints
   ! are represented by having their slacks left-hand with no artificial variable.
   ! >= constraints have their slack initially left-hand with a minus sign and their
   ! artificial variable handled implicitly during their first exchange.
   ! == constraints have their artificial variable initially left-hand.

   !----------------------------------------------
   !----- Phase 1 : Find a feasible solution -----
   !----------------------------------------------

   init = .true.
   phase1: do
      if (init) then
         if (ngeq + neq == 0) exit phase1
         init = .false.
         ! List of >= constraints whose slack has never been exchanged.
         in_basis = .true.
         ! Auxiliary objective function.
         A(m + 2, :) = -sum(A(nleq + 2:m + 1, :), dim=1)
      end if

      if (nl1 > 0) then
         ! Find the maximum coefficient of the auxiliary objective function.
         kp = maxloc(A(m + 2, admissible_columns(:nl1) + 1), dim=1)
         kp = admissible_columns(kp)
         bmax = A(m + 2, kp + 1)
      else
         bmax = 0.0_dp
      end if

      phase1a: do
         if ((bmax <= eps) .and. (A(m + 2, 1) < -eps)) then
            ! Auxiliary objective is still negative and can't be improved.
            ! No feasible solution exits.
            info = -1
            return
         else if ((bmax <= eps) .and. (A(m + 2, 1) <= eps)) then
            ! Auxiliary objective is zero and can't be improved. Feasible
            ! starting vector has been computed. Clean out the artificial
            ! variables corresponding to remaining equality constraints
            ! and exit phase one.
            do ip = nleq + ngeq + 1, m
               if (iposv(ip) == ip + n) then
                  if (nl1 > 0) then
                     kp = maxloc(abs(A(ip + 1, admissible_columns(:nl1) + 1)), dim=1)
                     kp = admissible_columns(kp)
                     bmax = A(ip + 1, kp + 1)
                  else
                     bmax = 0.0_dp
                  end if
                  if (bmax > eps) exit phase1a
               end if
            end do
            ! Change sign of row for any >= constraints still present from
            ! the initial basis.
            where (spread(in_basis, 2, n + 1) .eqv. .true.) &
               A(nleq + 2:nleq + ngeq + 1, :) = -A(nleq + 2:nleq + ngeq + 1, :)
            exit phase1
         end if

         ip = find_pivot(A, kp)
         if (ip == 0) then
            ! No pivot has been found. Maximum of the auxiliary function
            ! is unbounded.
            info = 1
            return
         end if

         exit phase1a
      end do phase1a

      ! Exchange a basic and non-basic variable.
      call exchange_variables(A, ip, kp, m + 1, n)

      if (iposv(ip) >= n + nleq + ngeq + 1) then
         ! Exchanged out an artifical variable for an equality constraint.
         ! Make sure it stays out by removing it from the l1 list.
         k = findloc(admissible_columns(:nl1), kp, dim=1)
         nl1 = nl1 - 1
         admissible_columns(k:nl1) = admissible_columns(k + 1:nl1 + 1)
      else
         kh = iposv(ip) - nleq - n
         ! Exchanged an >= constraints.
         if (kh >= 1) then
            ! If it's the first time, correct the pivot column for the minus
            ! sign and the implicit artifical variable.
            if (in_basis(kh)) then
               in_basis(kh) = .false.
               A(m + 2, kp + 1) = A(m + 2, kp + 1) + 1.0_dp
               A(:, kp + 1) = -A(:, kp + 1)
            end if
         end if
      end if

      ! Book-keeping.
      call swap(izrov(kp), iposv(ip)); exit phase1
   end do phase1

   !------------------------------------------------------
   !-----     Phase 2 : Compute optimal solution     -----
   !------------------------------------------------------

   niter = 0
   phase2: do while (niter <= maxiter)
      niter = niter + 1
      if (nl1 > 0) then
         kp = maxloc(A(1, admissible_columns(:nl1) + 1), dim=1)
         kp = admissible_columns(kp)
         bmax = A(1, kp + 1)
      else
         bmax = 0.0_dp
      end if

      if (bmax <= eps) then
         ! No more positive coefficient in the modified cost function.
         ! Solution is optimal.
         info = 0
         return
      end if

      ip = find_pivot(A, kp)
      if (ip == 0) then
         ! No pivot has been found. Objective function is unbounded.
         info = 1
         return
      end if

      ! Exchange a basic and non-basic variable.
      call exchange_variables(A, ip, kp, m, n)
      ! Book-keeping.
      call swap(izrov(kp), iposv(ip))
   end do phase2

   ! Return information flag info = 2 if the maximum number of
   ! iterations has been reached.
   if (niter > maxiter) then
      info = 2; return
   end if

   end procedure dense_standard_simplex

   pure integer(ilp) function find_pivot(A, k) result(ip)
      real(dp), intent(in) :: A(:, :)
        !! Simplex tableau
      integer(ilp), intent(in) :: k
        !! Current column

      ! Internal variables
      real(dp) :: q, q0, q1, qp
      integer(ilp) :: i, j

      associate (m => size(A, 1) - 2, n => size(A, 2) - 1)
         ! Determine if a pivot exist.
         i = findloc(A(2:m + 1, k + 1) < -eps, .true., dim=1)
         if (i > m) then
            ! No possible pivot. Problem is infeasible.
            ip = 0; return
         end if

         q1 = -A(i + 1, 1)/A(i + 1, k + 1); ip = i

         do i = ip + 1, m
            if (A(i + 1, k + 1) < -eps) then
               q = -A(i + 1, 1)/A(i + 1, k + 1)
               if (q < q1) then
                  ip = i; q1 = q
               else if (q == q1) then   ! Degeneracy situation.
                  do j = 1, n
                     qp = -A(i + 1, j + 1)/A(ip + 1, k + 1)
                     q0 = -A(i + 1, j + 1)/A(i + 1, k + 1)
                     if (q0 /= qp) exit
                  end do
                  if (q0 < qp) ip = i
               end if
            end if
         end do
      end associate

   end function find_pivot

   pure subroutine exchange_variables(A, ip, kp, i1, k1)
      real(dp), intent(inout) :: A(:, :)
      integer(ilp), intent(in) :: ip, kp, i1, k1
      integer(ilp) :: ip1, kp1, i, j
      real(dp) :: piv
      integer(ilp), dimension(k1) :: icol
      integer(ilp), dimension(i1) :: irow
      integer(ilp), dimension(max(i1, k1) + 1) :: itmp

      ip1 = ip + 1; kp1 = kp + 1
      piv = 1.0_dp/A(ip1, kp1)

      itmp(1:k1 + 1) = arange(k1 + 1)
      icol = pack(itmp(1:k1 + 1), itmp(1:k1 + 1) /= kp1)

      itmp(1:i1 + 1) = arange(i1 + 1)
      irow = pack(itmp(1:i1 + 1), itmp(1:i1 + 1) /= ip1)

      A(irow, kp1) = A(irow, kp1)*piv
      A(irow, icol) = A(irow, icol) - outer_product(A(irow, kp1), A(ip1, icol))

      A(ip1, icol) = -A(ip1, icol)*piv
      A(ip1, kp1) = piv
   end subroutine exchange_variables

end submodule lightconvex_lp
