submodule(lightconvex_lp) lightconvex_simplex
   use stdlib_math, only: arange, swap
   use stdlib_linalg, only: outer_product
   use stdlib_intrinsics, only: sum => stdlib_sum
   implicit none(external)
contains
   !============================================
   !============================================
   !=====                                  =====
   !=====     PRIMAL SIMPLEX ALGORITHM     =====
   !=====                                  =====
   !============================================
   !============================================

   !----------------------------------------
   !-----     HIGH-LEVEL INTERFACE     -----
   !----------------------------------------

   module procedure initialize_primal_simplex_alg
   alg%maxiter = 1000; if (present(maxiter)) alg%maxiter = maxiter
   alg%pivot = Dantzig(); if (present(pivot)) alg%pivot = pivot
   alg%initialization = auxiliary_function(); if (present(initialization)) alg%initialization = initialization
   end procedure initialize_primal_simplex_alg

   module procedure solve_with_dense_simplex
   integer(ilp), allocatable :: iposv(:)
   integer(ilp) :: i, m, n, info

   !> Problem dimension.
   m = size(problem%A, 1) - 2
   n = size(problem%A, 2) - 1

   !> Allocate variables.
   allocate (iposv(m), source=0)

   !> Solve the problem.
   call simplex(problem%A, problem%nleq, problem%ngeq, problem%neq, iposv, &
                alg%maxiter, info, alg%pivot, alg%initialization)

   !> Problem's status.
   call problem%set_status(info)

   if (is_optimal(problem)) then
      !> Extract primal variables and slacks from the tableau.
      allocate (solution%x(n), source=0.0_dp); allocate (solution%s(m), source=0.0_dp)

      do i = 1, m
         if (iposv(i) <= n) solution%x(iposv(i)) = problem%A(i + 1, 1)
         if (iposv(i) > n) solution%s(iposv(i) - n) = problem%A(i + 1, 1)
      end do

      solution%objective_value = problem%A(1, 1)
   end if
   end procedure solve_with_dense_simplex

   !-----------------------------------------------------------
   !-----     STANDARD SIMPLEX METHOD (DENSE TABLEAU)     -----
   !-----------------------------------------------------------

   module procedure dense_standard_simplex
   integer(ilp) :: m, n
    !! Dimensions of the problem.
    !!  - m : Number of constraints (excluding the non-negativity).
    !!  - n : Number of variables (excluding the slack ones).
   integer(ilp) :: niter
    !! Current iteration number
   integer(ilp), allocatable           :: izrov(:)
    !! Book-keeping for variables being zeroed-out.
   integer(ilp), dimension(size(A, 2)) :: nonbasics
    !! Indices of the non-basic variables.
   logical(lk), dimension(ngeq)        :: is_basic
    !! List of >= constraints currently in the basis.
    !! Method used to find an initial feasible point.
   integer(ilp) :: rowpiv, colpiv
    !! Index of the pivoting row and pivoting column.

   ! Miscellaneous
   integer(ilp) :: nl1
   real(dp)     :: reduced_cost

   !> Sanity checks.
   m = size(A, 1) - 2; n = size(A, 2) - 1
   call assert(assertion=m == nleq + ngeq + neq, &
               description="Tableau size is inconsistent with the number of constraints.")
   call assert(assertion=size(iposv) == m, &
               description="Dimension of iposv is inconsistant with the number of constraints.")
   ! Remove assertion as it seems to strict.
   ! call assert(assertion=all(A(2:, 1) >= -eps), &
   !             description="Constants b_i need to be non-negative.")

   !----- Initialization -----

   ! Index list of columns admissible for exchange.
   nl1 = n; nonbasics(:n) = arange(n)
   izrov = nonbasics(:n)    ! All variables are initially non-basic.
   iposv = n + arange(m)    ! Initial basic variables. <= constraints
   ! are represented by having their slacks basic with no artificial variable.
   ! >= constraints have their slack initially basic with a minus sign and their
   ! artificial variable handled implicitly during their first exchange.
   ! == constraints have their artificial variable initially basic.

   !--> Phase 1 : Find a feasible solution
   !---------------------------------------
   call find_feasible_point(initialization, A, nleq, ngeq, neq, &
                            nonbasics, iposv, izrov, pivot, nl1, info)
   if ((info == infeasible_status) .or. (info == unbounded_status)) return

   !--> Phase 2 : Compute optimal solution
   !--------------------------------------
   niter = 0
   do while (niter <= maxiter)

      ! Update iteration counter.
      niter = niter + 1

      ! Find the pivoting column according to the selected rule.
      colpiv = pivot_selection(pivot, A(1, nonbasics(:nl1) + 1), nonbasics, nl1)
      reduced_cost = merge(0.0_dp, A(1, colpiv + 1), colpiv == 0)

      if (reduced_cost <= tol) then
         ! No more positive coefficient in the reduced cost function.
         info = optimal_status; return
      end if

      ! Find the limiting row and make it the pivoting one.
      rowpiv = pivoting_row(A, colpiv)
      if (rowpiv == 0) then
         ! No pivot has been found. Objective function is unbounded.
         info = unbounded_status; return
      end if

      ! Exchange the basic and non-basic variables.
      call rank1_update(A(:m + 1, :n + 1), rowpiv + 1, colpiv + 1)

      ! Book-keeping.
      call swap(izrov(colpiv), iposv(rowpiv))
   end do

   ! Maximum number of iterations has been exceeded before
   ! an optimal solution has been found.
   if (niter > maxiter) info = maxiter_exceeded

   end procedure dense_standard_simplex

   !-----------------------------------------
   !-----     PIVOT SELECTION RULES     -----
   !-----------------------------------------

   !----- Main driver -----

   pure integer(ilp) function pivot_selection(rule, cost, nonbasics, nl1) result(colpiv)
      class(abstract_pivot_rule), intent(in) :: rule
        !! Pivot selection rule.
      real(dp), intent(in) :: cost(:)
        !! Reduced cost at the current iteration.
      integer(ilp), intent(in) :: nonbasics(:)
        !! List of current non-basic variables.
      integer(ilp), intent(in) :: nl1
        !! Yet to be understood.

      select type (rule)
      type is (Dantzig)
         colpiv = dantzig_pivot(cost, nonbasics, nl1)
      class default
         error stop "Selected rule is not implemented."
      end select
   end function pivot_selection

   !----- Dantzig's pivot -----
   pure integer(ilp) function dantzig_pivot(cost, nonbasics, nl1) result(colpiv)
      real(dp), intent(in) :: cost(:)
      integer(ilp), intent(in) :: nonbasics(:)
      integer(ilp), intent(in) :: nl1

      if (nl1 < 0) then
         colpiv = 0
      else
         colpiv = maxloc(cost, dim=1)
         colpiv = nonbasics(colpiv)
      end if
   end function dantzig_pivot

   !------------------------------------------
   !-----     INITIALIZATION METHODS     -----
   !------------------------------------------

   !----- Main driver -----

   pure subroutine find_feasible_point(method, A, nleq, ngeq, neq, &
                                       nonbasics, iposv, izrov, pivot, nl1, info)
      class(abstract_feasible_initialization), intent(in) :: method
      real(dp), intent(inout) :: A(:, :)
      integer(ilp), intent(in) :: nleq, ngeq, neq
      integer(ilp), intent(out) :: nonbasics(:)
      integer(ilp), intent(out) :: iposv(:), izrov(:)
      class(abstract_pivot_rule), intent(in) :: pivot
      integer(ilp), intent(out) :: nl1
      integer(ilp), intent(out) :: info

      select type (method)
      type is (auxiliary_function)
         call auxiliary_function_initialization(A, nleq, ngeq, neq, &
                                                nonbasics, iposv, izrov, pivot, nl1, info)
      class default
         error stop "Selected rule is not implemented."
      end select
   end subroutine find_feasible_point

   !----- Find initial feasible point using the Auxiliary function method -----

   pure subroutine auxiliary_function_initialization(A, nleq, ngeq, neq, &
                                                     nonbasics, iposv, izrov, pivot, nl1, info)
      real(dp), intent(inout) :: A(:, :)
        !! Simplex tableau with one extra row at the end for book-keeping of the
        !! auxiliary function.
      integer(ilp), intent(in) :: nleq, ngeq, neq
        !! Number of constraints of each type.
      integer(ilp), intent(out) :: nonbasics(:)
        !! Indices of the non-basic variables.
      integer(ilp), intent(out) :: iposv(:), izrov(:)
        !! Book-keeping for the variables being zeroed-out or not.
      class(abstract_pivot_rule), intent(in) :: pivot
        !! Pivot-rule used in the simplex algorithm.
      integer(ilp), intent(out) :: nl1
        !! Related to equality constraints. Yet to be fully understood.
      integer(ilp), intent(out) :: info
        !! Return flag:
        !!      - info = -1 : Problem is infeasible.
        !!      - info = 0  : Feasible initial point has been found.
      integer(ilp) :: m, n
        !! Dimensions of the problem.
        !!  - m : Number of constraints (excluding the non-negativity).
        !!  - n : Number of variables (excluding the slack ones).
      logical(lk), dimension(ngeq)        :: is_nonbasic
        !! Booking for the >= inequalities being basic or not.
      integer(ilp) :: rowpiv, colpiv
        !! Index of the pivoting row and pivoting column.

      ! Miscellaneous
      integer(ilp) :: i
      real(dp)     :: reduced_cost
      !----- Initialization -----

      m = size(A, 1) - 2; n = size(A, 2) - 1
      ! Index list of columns admissible for exchange.
      nl1 = n; nonbasics(:n) = arange(n)
      izrov = nonbasics(:n)    ! All variables are initially non-basic.
      iposv = n + arange(m)    ! Initial basic variables. <= constraints
      ! are represented by having their slacks basic with no artificial variable.
      ! >= constraints have their slack initially basic with a minus sign and their
      ! artificial variable handled implicitly during their first exchange.
      ! == constraints have their artificial variable initially basic.

      is_nonbasic = .true.   ! All >= are initially set as non-basic variables.

      ! Auxiliary objective function.
      A(m + 2, :) = -sum(A(nleq + 2:m + 1, :), dim=1)

      !--> Phase 1 : Find a feasible solution
      !---------------------------------------

      if (ngeq + neq == 0) then
         ! Origin is a feasible point.
         info = optimal_status; return
      end if

      phase_1: do

         colpiv = pivot_selection(pivot, A(m + 2, nonbasics(:nl1) + 1), nonbasics, nl1)
         reduced_cost = merge(0.0_dp, A(m + 2, colpiv + 1), colpiv == 0)

         if ((reduced_cost <= tol) .and. (A(m + 2, 1) < -tol)) then
            ! Auxiliary objective is still negative and can't be improved.
            info = infeasible_status; return
         end if

         phase_1a: block
            if ((reduced_cost <= tol) .and. (A(m + 2, 1) <= tol)) then
               ! Auxiliary objective is zero and can't be improved. Feasible
               ! starting vector has been computed. Clean out the artificial
               ! variables corresponding to remaining equality constraints
               ! and exit phase one.
               do i = nleq + ngeq + 1, m
                  if (iposv(i) == i + n) then
                     colpiv = pivot_selection(pivot, abs(A(i + 1, nonbasics(:nl1) + 1)), &
                                              nonbasics, nl1)
                     reduced_cost = merge(0.0_dp, A(i + 1, colpiv + 1), colpiv == 0)
                     if (reduced_cost > tol) exit phase_1a
                  end if
               end do
               exit phase_1
            end if
         end block phase_1a

         ! Find the limiting row and make it the pivoting one.
         rowpiv = pivoting_row(A, colpiv)
         if (rowpiv == 0) then
            ! No pivot has been found. Auxiliary function is unbounded.
            info = unbounded_status; return
         end if

         ! Exchange the basic and non-basic variables.
         call rank1_update(A, rowpiv + 1, colpiv + 1)

         if (iposv(rowpiv) >= n + nleq + ngeq + 1) then
            ! Exchanged out an artifical variable for an equality constraint.
            ! Make sure it stays out by removing it from the l1 list.
            i = findloc(nonbasics(:nl1), colpiv, dim=1); nl1 = nl1 - 1
            nonbasics(i:nl1) = nonbasics(i + 1:nl1 + 1)
         else
            i = iposv(rowpiv) - nleq - n
            ! Exchanged an >= constraints.
            if ((i >= 1) .and. (is_nonbasic(i))) then
               ! If it's the first time, correct the pivot column for the minus
               ! sign and the implicit artifical variable.
               is_nonbasic(i) = .false.
               A(m + 2, colpiv + 1) = A(m + 2, colpiv + 1) + 1.0_dp
               A(:, colpiv + 1) = -A(:, colpiv + 1)
            end if
         end if

         ! Book-keeping.
         call swap(izrov(colpiv), iposv(rowpiv))
      end do phase_1

      ! Succesfully finished the initialization.
      info = optimal_status

      ! Change sign of row for any >= constraints still present from
      ! the initial basis.
      where (spread(is_nonbasic, 2, n + 1) .eqv. .true.) &
         A(nleq + 2:nleq + ngeq + 1, :) = -A(nleq + 2:nleq + ngeq + 1, :)

   end subroutine auxiliary_function_initialization

   !-------------------------------------
   !-----     UTILITY FUNCTIONS     -----
   !-------------------------------------

   ! References:
   ! Numerical recipes
   pure integer(ilp) function pivoting_row(A, colpiv) result(rowpiv)
      real(dp), intent(in) :: A(:, :)
    !! Simplex tableau
      integer(ilp), intent(in) :: colpiv
    !! Current column

      ! Internal variables
      real(dp) :: q, q0, q1, qp
      integer(ilp) :: i, j

      associate (m => size(A, 1) - 2, n => size(A, 2) - 1)
         ! Determine if a pivot exist.
         i = findloc(A(2:m + 1, colpiv + 1) < -tol, .true., dim=1)
         if (i > m) then
            ! No possible pivot. Problem is infeasible.
            rowpiv = 0; return
         end if

         q1 = -A(i + 1, 1)/A(i + 1, colpiv + 1); rowpiv = i

         do i = rowpiv + 1, m
            if (A(i + 1, colpiv + 1) < -tol) then
               q = -A(i + 1, 1)/A(i + 1, colpiv + 1)
               if (q < q1) then
                  rowpiv = i; q1 = q
               else if (q == q1) then   ! Degeneracy situation.
                  do j = 1, n
                     qp = -A(i + 1, j + 1)/A(rowpiv + 1, colpiv + 1)
                     q0 = -A(i + 1, j + 1)/A(i + 1, colpiv + 1)
                     if (q0 /= qp) exit
                  end do
                  if (q0 < qp) rowpiv = i
               end if
            end if
         end do
      end associate

   end function pivoting_row

   pure subroutine rank1_update(A, rowpiv, colpiv)
      real(dp), intent(inout) :: A(:, :)
        !! Simplex tableau to be updated.
      integer(ilp), intent(in) :: rowpiv, colpiv
        !! Indices of the pivoting row and column.

      integer(ilp) :: m, n
        !! Dimension of the tableau.
      integer(ilp), dimension(size(A, 2) - 1) :: icol
        !! Column indices excluding the pivoting one.
      integer(ilp), dimension(size(A, 1) - 1) :: irow
        !! Row indices excluding the pivoting one.
      integer(ilp), dimension(max(size(A, 1) - 1, size(A, 2) - 1) + 1) :: itmp
        !! Temporary array.
      real(dp) :: piv
        !! Pivot value.

      !> Dimension of the matrix.
      m = size(A, 1); n = size(A, 2)

      !> Fetch the pivot value.
      piv = 1.0_dp/A(rowpiv, colpiv)

      !> Column indices excluding the pivoting one.
      itmp(:n) = arange(n); icol = pack(itmp(:n), itmp(:n) /= colpiv)

      !> Row indices excluding the pivoting one.
      itmp(:m) = arange(m); irow = pack(itmp(:m), itmp(:m) /= rowpiv)

      !> Rank-1 downdate.
      A(irow, colpiv) = piv*A(irow, colpiv)
      A(irow, icol) = A(irow, icol) - outer_product(A(irow, colpiv), A(rowpiv, icol))
      A(rowpiv, icol) = -piv*A(rowpiv, icol); A(rowpiv, colpiv) = piv
   end subroutine rank1_update

end submodule lightconvex_simplex
