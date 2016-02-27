! NOTE: you must switch off solid wall treatment by setting
!  wall = 0
! in solveFVM.f95 file.
subroutine init_cond_gresho(rho, vex, vey, pre, phi, phix, phiy)
   use comvar
   implicit none

   real    :: rho(-1:nx+2, -1:ny+2)
   real    :: vex(-1:nx+2, -1:ny+2)
   real    :: vey(-1:nx+2, -1:ny+2)
   real    :: pre(-1:nx+2, -1:ny+2)
   real    :: phi(-1:nx+2, -1:ny+2)
   real    :: phix(-1:nx+2, -1:ny+2)
   real    :: phiy(-1:nx+2, -1:ny+2)

   real    :: x, y, r, p0
   integer :: i, j

   ! Gresho
   xmin = -1.0
   ymin = -1.0
   xmax =  1.0
   ymax =  1.0

   dx = (xmax - xmin)/nx
   dy = (ymax - ymin)/ny

   final_time = 1.0
   xperiod = yes
   yperiod = yes

   ! Gresho
   p0 = 1.0 / GAMMA / 0.01**2 - 0.5;

   do i=-1,nx+2
      do j=-1,ny+2
         x = xmin + (i-1)*dx + 0.5*dx
         y = ymin + (j-1)*dy + 0.5*dy

         phi(i,j) = potential(x,y)

         ! gresho vortex
         r = sqrt(x*x + y*y)
         rho(i,j) = 1.0
         if (r < 0.2) then
            pre(i,j) = p0 + 12.5 * r**2
            vex(i,j) = -5.0 * y
            vey(i,j) =  5.0 * x
         else if (r < 0.4) then
            pre(i,j) = p0 + 12.5 * r**2 + 4.0*(1.0 - 5.0*r - log(0.2) + log(r))
            vex(i,j) = -(2.0 - 5.0*r) * y/r
            vey(i,j) =  (2.0 - 5.0*r) * x/r
         else 
            pre(i,j) = p0 + 4.0*log(2.0) - 2.0
            vex(i,j) = 0.0
            vey(i,j) = 0.0
         end if

      enddo
   enddo
end subroutine init_cond_gresho
