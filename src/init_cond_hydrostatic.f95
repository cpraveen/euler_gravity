subroutine init_cond_hydrostatic(rho, vex, vey, pre, phi, phix, phiy)
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

   xmin = 0.0
   xmax = 1.0
   ymin = 0.0
   ymax = 1.0

   ! Miczek test
   !xmin =-0.1
   !xmax = 0.1
   !ymin =-0.1
   !ymax = 0.1

   !xmin = -1.0
   !xmax = +1.0
   !ymin = -1.0
   !ymax = +1.0

   dx = (xmax - xmin)/nx
   dy = (ymax - ymin)/ny

   xperiod = yes
   yperiod = no
   final_time = 25.0

   rho0 = 1.21
   p0   = 1.0
   nu   = 1.2

   ! polytropic miczek
   alpha = 4.0

   ! Rayleigh-taylor parameters
   drho = 0.1
   r0   = 0.6
   alp  = exp(-r0)/(exp(-r0) + drho)
   ff   = exp(r0*(1.0-alp)/alp)

   do i=-1,nx+2
      do j=-1,ny+2
         x = xmin + (i-1)*dx + 0.5*dx
         y = ymin + (j-1)*dy + 0.5*dy

         phi(i,j) = potential(x,y)

         ! isothermal
         !rho(i,j) = rho0 * exp(-rho0/p0*phi(i,j))
         !pre(i,j) = p0 * exp(-rho0/p0*phi(i,j))

         ! polytropic
         T        = 1.0 - (nu-1.0)*phi(i,j)/nu
         rho(i,j) = T**(1.0/(nu-1.0))
         pre(i,j) = T**(nu/(nu-1.0))

         ! polytropic Miczek: unstable atmosphere
         !T        = 1.0 - alpha*phi(i,j)
         !pre(i,j) = T**(1.0/alpha)
         !rho(i,j) = pre(i,j)/T

         ! rayleigh taylor
         !r = sqrt(x*x + y*y)
         !phi(i,j) = r
         !phix(i,j) = x/r
         !phiy(i,j) = y/r
         !theta = atan2(y,x)
         !ri = r0*(1.0 + 0.02*cos(20.0*theta))
         !if(r.lt.r0)then
         !   pre(i,j) = exp(-r)
         !else
         !   pre(i,j) = ff*exp(-r/alp)
         !endif
         !if(r.lt.ri)then
         !   rho(i,j) = exp(-r)
         !else
         !   rho(i,j) = ff*exp(-r/alp)/alp
         !endif

         vex(i,j) = 0.0
         vey(i,j) = 0.0
      enddo
   enddo
end subroutine init_cond_hydrostatic

