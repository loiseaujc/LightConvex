module LightConvex
   implicit none(external)
   private

   public :: say_hello
contains
   subroutine say_hello
      print *, "Hello, LightConvex!"
   end subroutine say_hello
end module LightConvex
