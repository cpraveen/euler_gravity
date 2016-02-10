subroutine compute_error(rho, vex, vey, pre, rho0, pre0)
   use comvar
   implicit none

   real    :: rho(-1:nx+2, -1:ny+2)
   real    :: vex(-1:nx+2, -1:ny+2)
   real    :: vey(-1:nx+2, -1:ny+2)
   real    :: pre(-1:nx+2, -1:ny+2)
   real    :: rho0(-1:nx+2, -1:ny+2)
   real    :: pre0(-1:nx+2, -1:ny+2)

   integer :: i, j
   real    :: x, y, q, a, m
   real    :: err_rho, err_vex, err_vey, err_pre

   err_rho = 0.0
   err_vex = 0.0
   err_vey = 0.0
   err_pre = 0.0

   do j=1,ny
      do i=1,nx
         err_rho = err_rho + abs(rho(i,j) - rho0(i,j))
         err_vex = err_vex + abs(vex(i,j))
         err_vey = err_vey + abs(vey(i,j))
         err_pre = err_pre + abs(pre(i,j) - pre0(i,j))
      enddo
   enddo

   err_rho = err_rho / (nx*ny)
   err_vex = err_vex / (nx*ny)
   err_vey = err_vey / (nx*ny)
   err_pre = err_pre / (nx*ny)

   write(*,'(4E18.8)') err_rho, err_vex, err_vey, err_pre

   close(10)

end subroutine compute_error
