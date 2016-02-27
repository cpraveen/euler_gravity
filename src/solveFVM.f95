subroutine solveFVM(rho, vex, vey, pre, omg, co0, co1, res)

   use comvar

   implicit none

   real :: rho(-1:nx+2, -1:ny+2)
   real :: vex(-1:nx+2, -1:ny+2)
   real :: vey(-1:nx+2, -1:ny+2)
   real :: pre(-1:nx+2, -1:ny+2)
   real :: tem(-1:nx+2, -1:ny+2)
   real :: omg( 1:nx+1,  1:ny+1)
   real :: co0(4, -1:nx+2, -1:ny+2)
   real :: co1(4, -1:nx+2, -1:ny+2)
   real :: phi(-1:nx+2, -1:ny+2)
   real :: phix(-1:nx+2, -1:ny+2)
   real :: phiy(-1:nx+2, -1:ny+2)
   real :: res(4,0:nx+1,0:ny+1)
   real :: fd(4,nx+1,ny+1)
   real :: gd(4,nx+1,ny+1)

   real :: rho0(-1:nx+2, -1:ny+2)
   real :: pre0(-1:nx+2, -1:ny+2)

   integer :: it, i, j, rks
   real    :: lambda
   real    :: wim1(4), wi(4), wip1(4), wip2(4)
   real    :: wjm1(4), wj(4), wjp1(4), wjp2(4)
   real    :: betaimh, betaiph, betaip3h
   real    :: betajmh, betajph, betajp3h
   real    :: psiim1, psii, psiip1, psiip2, phiiph
   real    :: psijm1, psij, psijp1, psijp2, phijph
   real    :: priml(4), primr(4)
   real    :: s, sx, sy
   real    :: xflux(4), yflux(4)
   real    :: dxflux(4), dyflux(4)
   real    :: time
   real    :: resid(4), resid1(4)
   logical :: tostop

   integer :: wall, source

   wall   = 1  ! 0 = no wall, 1 = wall
   source = 1  ! 0 = non wb, 1 = wb

   ! set initial condition
   call init_cond(rho, vex, vey, pre, phi, phix, phiy)

   ! store stationary solution
   rho0 = rho
   pre0 = pre

   ! Perturb initial condition
   !call add_perturbation(rho, vex, vey, pre, phi)

   call prim2cons(rho, vex, vey, pre, co1)
   call periodic(co1)
   call cons2prim(co1, rho, vex, vey, pre, tem)
   call saveprim(0.0, rho, vex, vey, pre, rho0, pre0)
   call vorticity(rho, vex, vey, pre, omg)
   call savevort(0.0, omg)
   call timestep(rho, vex, vey, pre)

   time   = 0.0
   it     = 0

   do while(time < final_time .and. it < itmax)

      ! Exactly match final time
      tostop = .false.
      if(time + dt > final_time)then
         dt = final_time - time
         tostop = .true.
      endif

      lambda = dt/dx/dy

      co0(:,:,:) = co1(:,:,:)

      do rks=1,nrk

         res = 0.0

         ! x fluxes
         do i=0,nx
            do j=1,ny
               if(i.eq.0)then
                  betaiph = hbeta(tem(i+1,j), tem(i+1,j))
                  betaip3h= hbeta(tem(i+1,j), tem(i+2,j))

                  phiiph = 1.5*phi(i+1,j) - 0.5*phi(i+2,j)

                  psiip1 =-2.0*betaiph*(phi(i+1,j)-phiiph)
                  psiip2 =-2.0*betaiph*(phi(i+1,j)-phiiph) - 2.0*betaip3h*(phi(i+2,j)-phi(i+1,j))

                  wip1(1) = rho(i+1,j) * exp(-psiip1 * source)
                  wip1(2) = vex(i+1,j)
                  wip1(3) = vey(i+1,j)
                  wip1(4) = pre(i+1,j) * exp(-psiip1 * source)

                  wip2(1) = rho(i+2,j) * exp(-psiip2 * source)
                  wip2(2) = vex(i+2,j)
                  wip2(3) = vey(i+2,j)
                  wip2(4) = pre(i+2,j) * exp(-psiip2 * source)

                  primr = 1.5*wip1 - 0.5*wip2
                  priml = primr
                  if(wall.eq.1) priml(2) = -primr(2)
               else if(i.eq.1)then
                  betaiph = hbeta(tem(i  ,j), tem(i+1,j))
                  betaip3h= hbeta(tem(i+1,j), tem(i+2,j))

                  psii   = betaiph*(phi(i+1,j)-phi(i,j))
                  psiip1 =-betaiph*(phi(i+1,j)-phi(i,j))
                  psiip2 =-betaiph*(phi(i+1,j)-phi(i,j)) - 2.0*betaip3h*(phi(i+2,j)-phi(i+1,j))

                  wi(1) = rho(i,j) * exp(-psii * source)
                  wi(2) = vex(i,j)
                  wi(3) = vey(i,j)
                  wi(4) = pre(i,j) * exp(-psii * source)

                  wip1(1) = rho(i+1,j) * exp(-psiip1 * source)
                  wip1(2) = vex(i+1,j)
                  wip1(3) = vey(i+1,j)
                  wip1(4) = pre(i+1,j) * exp(-psiip1 * source)

                  wip2(1) = rho(i+2,j) * exp(-psiip2 * source)
                  wip2(2) = vex(i+2,j)
                  wip2(3) = vey(i+2,j)
                  wip2(4) = pre(i+2,j) * exp(-psiip2 * source)

                  priml = 0.5*(wi + wip1)
                  call reconstruct(wip2, wip1, wi, primr)
               else if(i.eq.nx-1)then
                  betaimh = hbeta(tem(i-1,j), tem(i,  j))
                  betaiph = hbeta(tem(i  ,j), tem(i+1,j))

                  psiim1 = 2.0*betaimh*(phi(i,j)-phi(i-1,j)) + betaiph*(phi(i+1,j)-phi(i,j))
                  psii   = betaiph*(phi(i+1,j)-phi(i,j))
                  psiip1 =-betaiph*(phi(i+1,j)-phi(i,j))

                  wim1(1) = rho(i-1,j) * exp(-psiim1 * source)
                  wim1(2) = vex(i-1,j)
                  wim1(3) = vey(i-1,j)
                  wim1(4) = pre(i-1,j) * exp(-psiim1 * source)

                  wi(1) = rho(i,j) * exp(-psii * source)
                  wi(2) = vex(i,j)
                  wi(3) = vey(i,j)
                  wi(4) = pre(i,j) * exp(-psii * source)

                  wip1(1) = rho(i+1,j) * exp(-psiip1 * source)
                  wip1(2) = vex(i+1,j)
                  wip1(3) = vey(i+1,j)
                  wip1(4) = pre(i+1,j) * exp(-psiip1 * source)

                  call reconstruct(wim1, wi, wip1, priml)
                  primr = 0.5*(wi + wip1)
               else if(i.eq.nx)then
                  betaimh = hbeta(tem(i-1,j), tem(i,  j))
                  betaiph = hbeta(tem(i  ,j), tem(i  ,j))

                  phiiph = 1.5*phi(i,j) - 0.5*phi(i-1,j)

                  psiim1 = 2.0*betaimh*(phi(i,j)-phi(i-1,j)) + 2.0*betaiph*(phiiph-phi(i,j))
                  psii   = 2.0*betaiph*(phiiph-phi(i,j))

                  wim1(1) = rho(i-1,j) * exp(-psiim1 * source)
                  wim1(2) = vex(i-1,j)
                  wim1(3) = vey(i-1,j)
                  wim1(4) = pre(i-1,j) * exp(-psiim1 * source)

                  wi(1) = rho(i,j) * exp(-psii * source)
                  wi(2) = vex(i,j)
                  wi(3) = vey(i,j)
                  wi(4) = pre(i,j) * exp(-psii * source)

                  priml = 1.5*wi - 0.5*wim1
                  primr = priml
                  if(wall.eq.1)primr(2) = -priml(2)
               else
                  betaimh = hbeta(tem(i-1,j), tem(i,  j))
                  betaiph = hbeta(tem(i  ,j), tem(i+1,j))
                  betaip3h= hbeta(tem(i+1,j), tem(i+2,j))

                  psiim1 = 2.0*betaimh*(phi(i,j)-phi(i-1,j)) + betaiph*(phi(i+1,j)-phi(i,j))
                  psii   = betaiph*(phi(i+1,j)-phi(i,j))
                  psiip1 =-betaiph*(phi(i+1,j)-phi(i,j))
                  psiip2 =-betaiph*(phi(i+1,j)-phi(i,j)) - 2.0*betaip3h*(phi(i+2,j)-phi(i+1,j))

                  wim1(1) = rho(i-1,j) * exp(-psiim1 * source)
                  wim1(2) = vex(i-1,j)
                  wim1(3) = vey(i-1,j)
                  wim1(4) = pre(i-1,j) * exp(-psiim1 * source)

                  wi(1) = rho(i,j) * exp(-psii * source)
                  wi(2) = vex(i,j)
                  wi(3) = vey(i,j)
                  wi(4) = pre(i,j) * exp(-psii * source)

                  wip1(1) = rho(i+1,j) * exp(-psiip1 * source)
                  wip1(2) = vex(i+1,j)
                  wip1(3) = vey(i+1,j)
                  wip1(4) = pre(i+1,j) * exp(-psiip1 * source)

                  wip2(1) = rho(i+2,j) * exp(-psiip2 * source)
                  wip2(2) = vex(i+2,j)
                  wip2(3) = vey(i+2,j)
                  wip2(4) = pre(i+2,j) * exp(-psiip2 * source)

                  call reconstruct(wim1, wi, wip1, priml)
                  call reconstruct(wip2, wip1, wi, primr)

               endif
               call numflux(1.0, 0.0, priml, primr, xflux, dxflux)
               res(:,i,j)   = res(:,i,j)   + dy*xflux(:)
               res(:,i+1,j) = res(:,i+1,j) - dy*xflux(:)
            enddo
         enddo

         ! y fluxes
         do j=0,ny
            do i=1,nx
               if(j.eq.0)then
                  betajph = hbeta(tem(i,j+1), tem(i,j+1))
                  betajp3h= hbeta(tem(i,j+1), tem(i,j+2))

                  phijph = 1.5*phi(i,j+1) - 0.5*phi(i,j+2)

                  psijp1 =-2.0*betajph*(phi(i,j+1)-phijph)
                  psijp2 =-2.0*betajph*(phi(i,j+1)-phijph) - 2.0*betajp3h*(phi(i,j+2)-phi(i,j+1))

                  wjp1(1) = rho(i,j+1) * exp(-psijp1 * source)
                  wjp1(2) = vex(i,j+1)
                  wjp1(3) = vey(i,j+1)
                  wjp1(4) = pre(i,j+1) * exp(-psijp1 * source)

                  wjp2(1) = rho(i,j+2) * exp(-psijp2 * source)
                  wjp2(2) = vex(i,j+2)
                  wjp2(3) = vey(i,j+2)
                  wjp2(4) = pre(i,j+2) * exp(-psijp2 * source)

                  primr = 1.5*wjp1 - 0.5*wjp2
                  priml = primr
                  if(wall.eq.1)priml(3) = -primr(3)
               else if(j.eq.1)then
                  betajph = hbeta(tem(i  ,j), tem(i,j+1))
                  betajp3h= hbeta(tem(i,j+1), tem(i,j+2))

                  psij   = betajph*(phi(i,j+1)-phi(i,j))
                  psijp1 =-betajph*(phi(i,j+1)-phi(i,j))
                  psijp2 =-betajph*(phi(i,j+1)-phi(i,j)) - 2.0*betajp3h*(phi(i,j+2)-phi(i,j+1))

                  wj(1) = rho(i,j) * exp(-psij * source)
                  wj(2) = vex(i,j)
                  wj(3) = vey(i,j)
                  wj(4) = pre(i,j) * exp(-psij * source)

                  wjp1(1) = rho(i,j+1) * exp(-psijp1 * source)
                  wjp1(2) = vex(i,j+1)
                  wjp1(3) = vey(i,j+1)
                  wjp1(4) = pre(i,j+1) * exp(-psijp1 * source)

                  wjp2(1) = rho(i,j+2) * exp(-psijp2 * source)
                  wjp2(2) = vex(i,j+2)
                  wjp2(3) = vey(i,j+2)
                  wjp2(4) = pre(i,j+2) * exp(-psijp2 * source)

                  priml = 0.5*(wj + wjp1)
                  call reconstruct(wjp2, wjp1, wj, primr)
               else if(j.eq.ny-1)then
                  betajmh = hbeta(tem(i,j-1), tem(i,  j))
                  betajph = hbeta(tem(i,j),   tem(i,j+1))

                  psijm1 = 2.0*betajmh*(phi(i,j)-phi(i,j-1)) + betaiph*(phi(i,j+1)-phi(i,j))
                  psij   = betajph*(phi(i,j+1)-phi(i,j))
                  psijp1 =-betajph*(phi(i,j+1)-phi(i,j))

                  wjm1(1) = rho(i,j-1) * exp(-psijm1 * source)
                  wjm1(2) = vex(i,j-1)
                  wjm1(3) = vey(i,j-1)
                  wjm1(4) = pre(i,j-1) * exp(-psijm1 * source)

                  wj(1) = rho(i,j) * exp(-psij * source)
                  wj(2) = vex(i,j)
                  wj(3) = vey(i,j)
                  wj(4) = pre(i,j) * exp(-psij * source)

                  wjp1(1) = rho(i,j+1) * exp(-psijp1 * source)
                  wjp1(2) = vex(i,j+1)
                  wjp1(3) = vey(i,j+1)
                  wjp1(4) = pre(i,j+1) * exp(-psijp1 * source)

                  call reconstruct(wjm1, wj, wjp1, priml)
                  primr = 0.5*(wj + wjp1)
               else if(j.eq.ny)then
                  betajmh = hbeta(tem(i,j-1), tem(i,j))
                  betajph = hbeta(tem(i,j),   tem(i,j))

                  phijph = 1.5*phi(i,j) - 0.5*phi(i,j-1)

                  psijm1 = 2.0*betajmh*(phi(i,j)-phi(i,j-1)) + 2.0*betajph*(phijph-phi(i,j))
                  psij   = 2.0*betajph*(phijph-phi(i,j))

                  wjm1(1) = rho(i,j-1) * exp(-psijm1 * source)
                  wjm1(2) = vex(i,j-1)
                  wjm1(3) = vey(i,j-1)
                  wjm1(4) = pre(i,j-1) * exp(-psijm1 * source)

                  wj(1) = rho(i,j) * exp(-psij * source)
                  wj(2) = vex(i,j)
                  wj(3) = vey(i,j)
                  wj(4) = pre(i,j) * exp(-psij * source)

                  priml = 1.5*wj - 0.5*wjm1
                  primr = priml
                  if(wall.eq.1) primr(3) = -priml(3)
               else
                  betajmh = hbeta(tem(i,j-1), tem(i,j))
                  betajph = hbeta(tem(i,j),   tem(i,j+1))
                  betajp3h= hbeta(tem(i,j+1), tem(i,j+2))

                  psijm1 = 2.0*betajmh*(phi(i,j)-phi(i,j-1)) + betajph*(phi(i,j+1)-phi(i,j))
                  psij   = betajph*(phi(i,j+1)-phi(i,j))
                  psijp1 =-betajph*(phi(i,j+1)-phi(i,j))
                  psijp2 =-betajph*(phi(i,j+1)-phi(i,j)) - 2.0*betajp3h*(phi(i,j+2)-phi(i,j+1))

                  wjm1(1) = rho(i,j-1) * exp(-psijm1 * source)
                  wjm1(2) = vex(i,j-1)
                  wjm1(3) = vey(i,j-1)
                  wjm1(4) = pre(i,j-1) * exp(-psijm1 * source)

                  wj(1) = rho(i,j) * exp(-psij * source)
                  wj(2) = vex(i,j)
                  wj(3) = vey(i,j)
                  wj(4) = pre(i,j) * exp(-psij * source)

                  wjp1(1) = rho(i,j+1) * exp(-psijp1 * source)
                  wjp1(2) = vex(i,j+1)
                  wjp1(3) = vey(i,j+1)
                  wjp1(4) = pre(i,j+1) * exp(-psijp1 * source)

                  wjp2(1) = rho(i,j+2) * exp(-psijp2 * source)
                  wjp2(2) = vex(i,j+2)
                  wjp2(3) = vey(i,j+2)
                  wjp2(4) = pre(i,j+2) * exp(-psijp2 * source)

                  call reconstruct(wjm1, wj, wjp1, priml)
                  call reconstruct(wjp2, wjp1, wj, primr)
               endif
               call numflux(0.0, 1.0, priml, primr, yflux, dyflux)
               res(:,i,j)   = res(:,i,j)   + dx*yflux(:)
               res(:,i,j+1) = res(:,i,j+1) - dx*yflux(:)
            enddo
         enddo

         ! Add gravitational source terms
         if(source.eq.0)then
            do i=1,nx
               do j=1,ny
                  sx = phix(i,j) ! x derivative of phi
                  sy = phiy(i,j) ! y derivative of phi
                  res(2,i,j) = res(2,i,j) + dx * dy * sx * rho(i,j)
                  res(3,i,j) = res(3,i,j) + dx * dy * sy * rho(i,j)
                  res(4,i,j) = res(4,i,j) + dx * dy * sx * rho(i,j) * vex(i,j) &
                                          + dx * dy * sy * rho(i,j) * vey(i,j)
               enddo
            enddo

         else
            do i=1,nx
               do j=1,ny
                  ! x derivative terms
                  if(i.eq.1)then
                     betaimh = hbeta(tem(i,j), tem(i,  j))
                     betaiph = hbeta(tem(i,j), tem(i+1,j))
                     s = exp(betaiph*(phi(i,j)-phi(i+1,j))) - exp(betaimh*(phi(i+1,j)-phi(i,j)))
                  else if(i.eq.nx)then
                     betaimh = hbeta(tem(i-1,j), tem(i,j))
                     betaiph = hbeta(tem(i,  j), tem(i,j))
                     s = exp(betaiph*(phi(i-1,j)-phi(i,j))) - exp(betaimh*(phi(i,j)-phi(i-1,j)))
                  else
                     betaimh = hbeta(tem(i-1,j), tem(i,j))
                     betaiph = hbeta(tem(i+1,j), tem(i,j))
                     s = exp(betaiph*(phi(i,j)-phi(i+1,j))) - exp(betaimh*(phi(i,j)-phi(i-1,j)))
                  endif

                  res(2,i,j) = res(2,i,j) - dy * s * pre(i,j)
                  res(4,i,j) = res(4,i,j) - dy * s * pre(i,j) * vex(i,j)

                  ! y derivative terms
                  if(j.eq.1)then
                     betajmh = hbeta(tem(i,j), tem(i,j))
                     betajph = hbeta(tem(i,j), tem(i,j+1))
                     s = exp(betajph*(phi(i,j)-phi(i,j+1))) - exp(betajmh*(phi(i,j+1)-phi(i,j)))
                  else if(j.eq.ny)then
                     betajmh = hbeta(tem(i,j-1), tem(i,j))
                     betajph = hbeta(tem(i,j),   tem(i,j))
                     s = exp(betajph*(phi(i,j-1)-phi(i,j))) - exp(betajmh*(phi(i,j)-phi(i,j-1)))
                  else
                     betajmh = hbeta(tem(i,j-1), tem(i,j))
                     betajph = hbeta(tem(i,j+1), tem(i,j))
                     s = exp(betajph*(phi(i,j)-phi(i,j+1))) - exp(betajmh*(phi(i,j)-phi(i,j-1)))
                  endif

                  res(3,i,j) = res(3,i,j) - dx * s * pre(i,j)
                  res(4,i,j) = res(4,i,j) - dx * s * pre(i,j) * vey(i,j)

                  !write(*,'(2I6,4E12.4)')i,j,res(:,i,j)
               enddo
            enddo
         endif

         ! update conserved variables
         resid = 0.0
         do i=1,nx
            do j=1,ny
               co1(:,i,j) = ark(rks)*co0(:,i,j) + &
                            (1.0-ark(rks))*(co1(:,i,j) - lambda*res(:,i,j))
               resid = resid + res(:,i,j)**2
            enddo
         enddo
         resid = sqrt(resid)

         call periodic(co1)
         call cons2prim(co1, rho, vex, vey, pre, tem)

      enddo ! Rk stage loop

      it = it + 1
      if(it==1)then
         resid1 = 1.0
         if(resid(1) > 0.0) resid1(1) = resid(1)
         if(resid(2) > 0.0) resid1(2) = resid(2)
         if(resid(3) > 0.0) resid1(3) = resid(3)
         if(resid(4) > 0.0) resid1(4) = resid(4)
      endif

      time = time + dt
      write(*,'(I6,F14.6,4E12.4)')it,time,resid(:)

      if(mod(it,itsave)==0 .or. it==itmax .or. tostop)then
         call cons2prim(co1,rho,vex,vey,pre,tem)
         call saveprim(time, rho, vex, vey, pre, rho0, pre0)
         call vorticity(rho, vex, vey, pre, omg)
         call savevort(time, omg)
      endif

   enddo ! time iteration loop

   call compute_error(rho, vex, vey, pre, rho0, pre0)


end subroutine solveFVM
