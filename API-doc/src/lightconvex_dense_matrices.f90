module lightconvex_dense_matrices
   use lightconvex_constants, only: ilp, dp, lk
   use lightconvex_abstract, only: abstract_cvx_vector, abstract_vector_rdp, &
                                   AbstractMatrix, AbstractSymMatrix
   use lightconvex_dense_vectors, only: dense_vector
   use stdlib_linalg_blas, only: gemv, symv
   use stdlib_optval, only: optval
   implicit none(type, external)
   private

   !> Base type for dense matrices.
   type, public, extends(AbstractMatrix) :: dense_matrix
      real(dp), allocatable :: data(:, :)
   contains
      procedure, pass(self) :: matvec => dense_matvec
   end type dense_matrix

   interface dense_matrix
      type(dense_matrix) module function create_constant_matrix(m, n, val) result(matrix)
         implicit none(type, external)
         integer(ilp), intent(in) :: m, n
         real(dp), optional, intent(in) :: val
      end function create_constant_matrix

      type(dense_matrix) module function create_from_array(array) result(matrix)
         implicit none(type, external)
         real(dp), intent(in) :: array(:, :)
      end function create_from_array
   end interface dense_matrix

   !> Base type for symmetric dense matrices.
   type, public, extends(AbstractSymMatrix) :: dense_sym_matrix
      real(dp), allocatable :: data(:, :)
   contains
      procedure, pass(self) :: matvec => dense_sym_matvec
   end type dense_sym_matrix

   interface dense_sym_matrix
      type(dense_sym_matrix) module function create_constant_sym_matrix(n, val) result(matrix)
         implicit none(type, external)
         integer(ilp), intent(in) :: n
         real(dp), optional, intent(in) :: val
      end function create_constant_sym_matrix

      type(dense_sym_matrix) module function create_from_sym_array(array) result(matrix)
         implicit none(type, external)
         real(dp), intent(in) :: array(:, :)
      end function create_from_sym_array
   end interface dense_sym_matrix

contains

   !-----------------------------
   !-----     FACTORIES     -----
   !-----------------------------

   module procedure create_constant_matrix
   if (any([m, n] <= 0)) error stop "dense_matrix: m and n need to be positive."
   allocate (matrix%data(m, n), source=optval(val, 0.0_dp))
   end procedure create_constant_matrix

   module procedure create_from_array
   allocate (matrix%data, source=array)
   end procedure create_from_array

   module procedure create_constant_sym_matrix
   if (n <= 0) error stop "dense_sym_matrix: n needs to be positive."
   allocate (matrix%data(n, n), source=optval(val, 0.0_dp))
   end procedure create_constant_sym_matrix

   module procedure create_from_sym_array
   if (size(array, 1) /= size(array, 2)) error stop "dense_sym_matrix: array needs to be square."
   allocate (matrix%data, source=array)
   end procedure create_from_sym_array

   !------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR DENSE MATRICES     -----
   !------------------------------------------------------------

   subroutine dense_matvec(self, alpha, x, beta, y, op)
      implicit none(type, external)
      class(dense_matrix), intent(in) :: self
      real(dp), intent(in) :: alpha, beta
      class(abstract_cvx_vector), intent(in) :: x
      class(abstract_cvx_vector), intent(inout) :: y
      character(1), intent(in) :: op

      select type (x)
      type is (dense_vector)
         select type (y)
         type is (dense_vector)

            associate (m => size(self%data, 1), n => size(self%data, 2))
               block
                  !> Sanity checks.
                  if (all(op /= ["n", "t"])) then
                     error stop "matvec: op needs be 'n' or 't'."
                  end if

                  if (((op == "n") .and. (size(y%data) /= m)) &
                      .or. ((op == "t") .and. (size(y%data) /= n))) then
                     error stop "matvec: y and A have inconsistent dimensions."
                  end if

                  if (((op == "n") .and. (size(x%data) /= n)) &
                      .or. ((op == "t") .and. (size(x%data) /= m))) then
                     error stop "matvec: x and A have inconsistent dimensions."
                  end if

                  !> General matrix-vector product.
                  call gemv(op, m, n, alpha, self%data, m, x%data, 1, beta, y%data, 1)
               end block
            end associate

         class default
            error stop "matvec: y needs to be a dense_vector."
         end select
      class default
         error stop "matvec: x needs to be a dense_vector."
      end select
   end subroutine dense_matvec

   !----------------------------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR DENSE SYMMETRIC MATRICES     -----
   !----------------------------------------------------------------------

   subroutine dense_sym_matvec(self, alpha, x, beta, y)
      implicit none(type, external)
      class(dense_sym_matrix), intent(in) :: self
      real(dp), intent(in) :: alpha, beta
      class(abstract_cvx_vector), intent(in) :: x
      class(abstract_cvx_vector), intent(inout) :: y

      select type (x)
      type is (dense_vector)
         select type (y)
         type is (dense_vector)

            associate (n => size(self%data, 1), uplo => "L")
               block

                  !> Sanity checks.
                  if (size(y%data) /= n) then
                     error stop "matvec: y and A have inconsistent dimensions."
                  end if

                  if (size(x%data) /= n) then
                     error stop "matvec: x and A have inconsistent dimensions."
                  end if

                  !> Symmetric matrix-vector product.
                  call symv(uplo, n, alpha, self%data, n, x%data, 1, beta, y%data, 1)
               end block
            end associate

         class default
            error stop "matvec: y needs to be a dense_vector."
         end select
      class default
         error stop "matvec: x needs to be a dense_vector."
      end select
   end subroutine dense_sym_matvec
end module lightconvex_dense_matrices
