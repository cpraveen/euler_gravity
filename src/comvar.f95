module comvar
   implicit none

   integer :: nx, ny
   real    :: xmin, xmax, ymin, ymax, dx, dy, dt
   integer :: itmax
   integer :: itsave
   real    :: cfl
   real    :: final_time

   real    :: mu, Prandtl, gas_const, kth, Cp

   integer :: nrk
   real    :: ark(3)

   real :: gamma = 1.4
   real :: M_PI = 4.0*atan(1.0)

   integer :: fileid_sol, fileid_omg

   integer :: fluxtype
   integer :: iroe=1, irusanov=2, imiczek=3

   integer :: ikepes_diss

   integer :: limtype
   integer :: ford=0, muscl3=1, mmod=2

   integer :: scheme
   integer :: fvm=1, gmd=2, kep=3, mvf=4

   integer :: vconf

   integer :: no=0, yes=1

   integer :: xperiod, yperiod

   contains

      real function logavg(a,b)
         implicit none
         real :: a, b

         real :: xi, f, u, u2, u3, ff

         xi = b/a
         f = (xi-1.0)/(xi+1.0)
         u = f*f

         if(u.lt.1.0e-2)then
            u2 = u*u
            u3 = u2*u
            ff = 1.0 + u/3.0 + u2/5.0 + u3/7.0
         else
            ff = 0.5*log(xi)/f
         endif

         logavg = 0.5*(a+b)/ff
      end function logavg

      real function hbeta(tl, tr)
         implicit none
         real :: tl, tr

         hbeta = logavg(tl, tr)
         hbeta = 0.5/(gas_const*hbeta)

      end function hbeta

      real function potential(x,y)
         implicit none
         real :: x, y

         potential = y
         !potential = x + y
      end function potential

end module comvar
