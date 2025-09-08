module lightconvex_linalg
   use assert_m, only: assert => assert_always
   use lightconvex_constants, only: ilp, dp, lk
   use stdlib_optval, only: optval
   use stdlib_linalg_lapack, only: geqr, gemqr, trmv
   implicit none(external)
   private

   !----- Derived-types -----!

   !> Matrix Factorization.
   type, abstract :: AbstractMatrixFactorization
   end type AbstractMatrixFactorization

   !> QR factorization.
   type, extends(AbstractMatrixFactorization), public :: qr_type
      private
      real(dp), allocatable :: data(:, :)
        !! Lapack storage of the QR factorization.
      real(dp), allocatable :: t(:)
        !! Data structuted used to represent Q in geqr.
      integer(ilp) :: tsize
        !! Dimension of the array T.
   end type

   interface QR
      type(qr_type) module function qr_fact(A) result(F)
         implicit none(external)
         real(dp), intent(in) :: A(:, :)
        !! Matrix to be factorized.
      end function qr_fact
   end interface
   public :: QR

   interface factmv
      module subroutine qrmv(A, x, y, trans)
         implicit none(external)
         type(qr_type), intent(inout) :: A
         real(dp), intent(in) :: x(:)
         character(len=*), intent(in), optional :: trans
         real(dp), intent(out), target, contiguous :: y(:)
      end subroutine qrmv
   end interface factmv
   public :: factmv

contains

   module procedure qr_fact
   integer(ilp) :: m, n, lda, lwork, info, tsize
   real(dp), allocatable :: t(:), work(:)

   !> Matrix dimension.
   m = size(A, 1); n = size(A, 2); lda = m

   !> Initialize derived-type.
   F%data = A

   !> Workspace query.
   lwork = -1; tsize = -1
   allocate (work(1), source=0.0_dp); allocate (t(1), source=0.0_dp)
   call geqr(m, n, F%data, m, t, tsize, work, lwork, info)
   call assert(assertion=info == 0, &
               description="Error during workspace query for geqr.")

   !> QR factorization.
   lwork = work(1); deallocate (work); allocate (work(lwork), source=0.0_dp)
   tsize = t(1); deallocate (t); allocate (t(tsize), source=0.0_dp)
   call geqr(m, n, F%data, m, t, tsize, work, lwork, info)
   call assert(assertion=info == 0, &
               description="Error while computing QR factorization in geqr.")

   !> Store additional data in derived type.
   F%t = t; F%tsize = tsize

   end procedure qr_fact

   module procedure qrmv
   character(len=1), parameter :: side = "L"
   integer(ilp), parameter :: nc = 1
   character(len=1) :: trans_
   integer(ilp) :: ma, na, mc
   integer(ilp) :: lda, ldc, lwork, info
   real(dp), allocatable :: work(:)
   real(dp), pointer :: ymat(:, :)

   ! By default A is not transposed.
   trans_ = optval(trans, "N")

   ma = size(A%data, 1); na = size(A%data, 2); mc = size(y)

   select case (trans_)
   case ("N")
      y(:na) = x(:na)
      !> y = R @ x
      call trmv("u", trans_, "n", na, A%data, ma, y, 1)
      !----- y = Q @ (R @ x) -----
      !> Pointer trick.
      ymat(1:na, 1:1) => y(:na)
      !> Workspace query.
      lwork = -1; allocate (work(1), source=0.0_dp)
      call gemqr(side, trans_, ma, nc, ma, A%data, ma, A%t, A%tsize, ymat, mc, work, lwork, info)
      call assert(assertion=info == 0, &
                  description="Error during workspace query for gemqr.")
      !> Actual matrix-vector product.
      lwork = work(1); deallocate (work); allocate (work(lwork), source=0.0_dp)
      call gemqr(side, trans_, ma, nc, ma, A%data, ma, A%t, A%tsize, ymat, mc, work, lwork, info)
      call assert(assertion=info == 0, &
                  description="Error while computing matrix-vector product in gemqr.")
   case ("T")
      block
         real(dp), allocatable, target :: x_tmp(:)
         real(dp), pointer :: xmat(:, :)
         !----- y = Q.T @ x -----
         !> Pointer trick.
         x_tmp = x; xmat(1:ma, 1:1) => x_tmp
         !> Workspace query.
         lwork = -1; allocate (work(1), source=0.0_dp)
         call gemqr(side, trans_, ma, nc, ma, A%data, ma, A%t, A%tsize, xmat, ma, work, lwork, info)
         call assert(assertion=info == 0, &
                     description="Error during workspace query for gemqr.")
         !> Actual matrix-vector product.
         lwork = work(1); deallocate (work); allocate (work(lwork), source=0.0_dp)
         call gemqr(side, trans_, ma, nc, ma, A%data, ma, A%t, A%tsize, xmat, ma, work, lwork, info)
         call assert(assertion=info == 0, &
                     description="Error while computing matrix-vector product in gemqr.")
         !----- y = R.T @ (Q.T @ x) -----
         y = x_tmp(:na)
         call trmv("u", trans_, "n", na, A%data, ma, y, 1)
      end block
   case default
      call assert(assertion=.false., &
                  description="Invalid value for trans.")
   end select
   end procedure
end module lightconvex_linalg
