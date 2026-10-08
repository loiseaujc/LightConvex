module lightconvex_dense_vectors
   use lightconvex_constants, only: ilp, dp, lk
   use lightconvex_abstract, only: abstract_cvx_vector, abstract_vector_rdp
   use stdlib_optval, only: optval
   use stdlib_math, only: stdlib_clip => clip
   use stdlib_linalg_blas, only: blas_scal => scal, blas_axpy => axpy
   use stdlib_linalg, only: norm
   use stdlib_intrinsics, only: stdlib_dot_product_kahan
   use stdlib_stats_distribution_normal, only: rvs_normal
   implicit none(type, external)
   private

   !> Base type for dense vectors.
   type, public, extends(abstract_cvx_vector) :: dense_vector
      real(dp), allocatable :: data(:)
   contains
      !> LightKrylov abstract vector required components.
      procedure, pass(self) :: get_size
      procedure, pass(self) :: zero
      procedure, pass(self) :: scal
      procedure, pass(self) :: axpby
      procedure, pass(self) :: dot
      procedure, pass(self) :: rand
      !> LightConvex abstract cvx vector required components.
      procedure, pass(self) :: norm_inf
      procedure, pass(self) :: hadamard
      procedure, pass(self) :: reciprocal
      procedure, pass(self) :: fill
      procedure, pass(self) :: clip
      procedure, pass(self) :: mask_le
   end type dense_vector

   interface dense_vector
      type(dense_vector) module function create_constant_vector(n, val) result(v)
         implicit none(type, external)
         integer(ilp), intent(in) :: n
         real(dp), optional, intent(in) :: val
      end function create_constant_vector

      type(dense_vector) module function create_from_array(x) result(v)
         implicit none(type, external)
         real(dp), intent(in) :: x(:)
      end function create_from_array
   end interface dense_vector

contains

   !-----------------------------
   !-----     FACTORIES     -----
   !-----------------------------

   module procedure create_constant_vector
   allocate (v%data(n), source=optval(val, 0.0_dp))
   end procedure create_constant_vector

   module procedure create_from_array
   allocate (v%data, source=x)
   end procedure create_from_array

   !-----------------------------------------------------
   !-----     TYPE-BOUND PROCEDURES FOR VECTORS     -----
   !-----------------------------------------------------

   integer(ilp) pure function get_size(self) result(n)
      implicit none(type, external)
      class(dense_vector), intent(in) :: self
      n = 0_ilp; if (allocated(self%data)) n = size(self%data, kind=ilp)
   end function get_size

   pure subroutine zero(self)
      implicit none(type, external)
      class(dense_vector), intent(inout) :: self
      self%data = 0.0_dp
   end subroutine zero

   pure subroutine scal(self, alpha)
      implicit none(type, external)
      class(dense_vector), intent(inout) :: self
      real(dp), intent(in) :: alpha
      call blas_scal(self%get_size(), alpha, self%data, 1)
   end subroutine scal

   pure subroutine axpby(alpha, vec, beta, self)
      implicit none(type, external)
      real(dp), intent(in) :: alpha, beta
      class(abstract_vector_rdp), intent(in) :: vec
      class(dense_vector), intent(inout) :: self
      select type (vec)
      type is (dense_vector)
         associate (n => self%get_size(), m => vec%get_size())
            !> Sanity check.
            if (n /= m) error stop "axpby: Vectors have inconsistent dimensions."
            if (beta == 0.0_dp) then
               call self%fill(beta)
            else if (beta /= 1.0_dp) then
               call self%scal(beta)
            end if
            call blas_axpy(n, alpha, vec%data, 1, self%data, 1)
         end associate
      class default
         error stop "axpby: self and vec have inconsistent types."
      end select
   end subroutine axpby

   real(dp) pure function dot(self, vec) result(alpha)
      implicit none(type, external)
      class(dense_vector), intent(in) :: self
      class(abstract_vector_rdp), intent(in) :: vec
      select type (vec)
      type is (dense_vector)
         associate (n => self%get_size(), m => vec%get_size())
            !> Sanity check.
            if (n /= m) error stop "dot: Vectors have inconsistent dimensions."
            alpha = stdlib_dot_product_kahan(self%data, vec%data)
         end associate
      class default
         error stop "dot: self and vec have inconsistent types."
      end select
   end function dot

   subroutine rand(self, ifnorm)
      implicit none(type, external)
      class(dense_vector), intent(inout) :: self
      logical(lk), optional, intent(in) :: ifnorm
      logical(lk) :: normalize
      self%data = rvs_normal(0.0_dp, 1.0_dp, array_size=self%get_size())
      normalize = optval(ifnorm, .false.)
      if (normalize) call self%scal(1.0_dp/norm(self%data, 2))
   end subroutine rand

   real(dp) pure function norm_inf(self) result(val)
      implicit none(type, external)
      class(dense_vector), intent(in) :: self
      val = norm(self%data, "inf")
   end function norm_inf

   pure subroutine hadamard(self, d)
      implicit none(type, external)
      class(dense_vector), intent(inout) :: self
      class(abstract_cvx_vector), intent(in) :: d
      select type (d)
      type is (dense_vector)
         associate (n => self%get_size(), m => d%get_size())
            !> Sanity check.
            if (n /= m) error stop "hadamard: Vectors have inconsistent dimensions."
            self%data = self%data*d%data
         end associate
      class default
         error stop "hadamard: self and d have inconsistent types."
      end select
   end subroutine hadamard

   pure subroutine reciprocal(self)
      implicit none(type, external)
      class(dense_vector), intent(inout) :: self
      self%data = 1.0_dp/self%data
   end subroutine reciprocal

   pure subroutine fill(self, val)
      implicit none(type, external)
      class(dense_vector), intent(inout) :: self
      real(dp), intent(in) :: val
      self%data = val
   end subroutine fill

   pure subroutine clip(self, l, u)
      implicit none(type, external)
      class(dense_vector), intent(inout) :: self
      class(abstract_cvx_vector), intent(in) :: l, u
      select type (l)
      type is (dense_vector)
         select type (u)
         type is (dense_vector)
            self%data = stdlib_clip(self%data, l%data, u%data)
         class default
            error stop "clip: u needs to be a dense_vector."
         end select
      class default
         error stop "clip: l needs to be a dense_vector."
      end select
   end subroutine clip

   pure subroutine mask_le(self, val)
      implicit none(type, external)
      class(dense_vector), intent(inout) :: self
      real(dp), intent(in) :: val
      where (self%data <= val)
         self%data = 1.0_dp
      else where
         self%data = 0.0_dp
      end where
   end subroutine mask_le

end module lightconvex_dense_vectors
