subroutine init_cond(rho, vex, vey, pre, phi, phix, phiy)
   use comvar
   implicit none

   real    :: rho(-1:nx+2, -1:ny+2)
   real    :: vex(-1:nx+2, -1:ny+2)
   real    :: vey(-1:nx+2, -1:ny+2)
   real    :: pre(-1:nx+2, -1:ny+2)
   real    :: phi(-1:nx+2, -1:ny+2)
   real    :: phix(-1:nx+2, -1:ny+2)
   real    :: phiy(-1:nx+2, -1:ny+2)

   real    :: x, y, r, theta, ri
   real    :: rho0, p0, nu, T, alpha
   real    :: drho, r0, alp, ff
   integer :: i, j

   !call init_cond_hydrostatic(rho, vex, vey, pre, phi, phix, phiy)
   call init_cond_gresho(rho, vex, vey, pre, phi, phix, phiy)

end subroutine init_cond

subroutine add_perturbation(rho, vex, vey, pre, phi)
   use comvar
   implicit none

   real    :: rho(-1:nx+2, -1:ny+2)
   real    :: vex(-1:nx+2, -1:ny+2)
   real    :: vey(-1:nx+2, -1:ny+2)
   real    :: pre(-1:nx+2, -1:ny+2)
   real    :: phi(-1:nx+2, -1:ny+2)

   real    :: x, y
   real    :: rho0, p0
   integer :: i, j

   rho0 = 1.21
   p0   = 1.0

   do i=-1,nx+2
      do j=-1,ny+2
         x = xmin + (i-1)*dx + 0.5*dx
         y = ymin + (j-1)*dy + 0.5*dy

         pre(i,j) = pre(i,j) + 1.0e-1 * exp(-100.0*rho0/p0*((x-0.3)**2 + (y-0.3)**2))
      enddo
   enddo
end subroutine add_perturbation
