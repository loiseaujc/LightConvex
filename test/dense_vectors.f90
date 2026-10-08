module TestDenseVectors
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use stdlib_math, only: all_close, is_close
   use stdlib_linalg, only: norm
   use lightconvex_constants, only: ilp, dp
   use lightconvex, only: dense_vector
   implicit none(external)
   private

   public :: collect_dense_vectors_tests
contains
   subroutine collect_dense_vectors_tests(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [new_unittest("Dense vectors testsuite", test_dense_vectors_testsuite)]
   end subroutine collect_dense_vectors_tests

   subroutine test_dense_vectors_testsuite(error)
      type(error_type), allocatable, intent(out) :: error

      !> Check allocation from source is handled properly.
      block
         type(dense_vector), allocatable :: x, y
         x = dense_vector(10_ilp, 2.0_dp)
         allocate (y, source=x)
         call check(error, all_close(x%data, y%data))
         if (allocated(error)) return
      end block

      !> Check the elementary operations for a vector space.
      block
         integer(ilp), parameter :: n = 10
         type(dense_vector), allocatable :: x, y, z
         !> Allocate variables.
         x = dense_vector(n)
         y = dense_vector(n)

         !> Check get size.
         call check(error, x%get_size() == n)
         if (allocated(error)) return

         !> Check zero vector.
         block
            real(dp), parameter :: zero_vector(n) = 0.0_dp
            call x%zero()
            call check(error, all_close(x%data, zero_vector))
            if (allocated(error)) return
         end block

         !> Create random vectors and check unit-norm.
         call x%rand(ifnorm=.true.)
         call y%rand(ifnorm=.true.)
         call check(error, is_close(x%norm(), 1.0_dp))
         if (allocated(error)) return

         !> Check vector addition.
         allocate (z, source=x)
         call z%add(y)
         call check(error, all_close(z%data, x%data + y%data))
         if (allocated(error)) return

         !> Check vector subtraction.
         z = x; call z%sub(y)
         call check(error, all_close(z%data, x%data - y%data))
         if (allocated(error)) return

         !> Check scaling.
         z = x; call z%scal(10.0_dp)
         call check(error, all_close(z%data, 10.0_dp*x%data))
         if (allocated(error)) return

         !> Check dot.
         call check(error, is_close(z%dot(x), dot_product(z%data, x%data)))
         if (allocated(error)) return

         !> Check infinity norm.
         call check(error, is_close(x%norm_inf(), norm(x%data, "inf")))
         if (allocated(error)) return

         !> Check reciprocal.
         z = x; call z%reciprocal()
         call check(error, all_close(z%data, 1.0_dp/x%data))
         if (allocated(error)) return

         !> Check Hadamard product.
         z = x; call z%hadamard(y)
         call check(error, all_close(z%data, x%data*y%data))
         if (allocated(error)) return

         !> Check fill.
         block
            real(dp), parameter :: fill_vector(n) = 10.0_dp
            call z%fill(10.0_dp)
            call check(error, all_close(z%data, fill_vector))
            if (allocated(error)) return
         end block

      end block
   end subroutine test_dense_vectors_testsuite
end module TestDenseVectors
