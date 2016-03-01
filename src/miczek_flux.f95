! Computes Miczek flux in the direction (normalX,normalY)
! (normalX,normalY)      = unit normal vector to cell face
!                vector points from left state towards right state
! primitiveLeft, primitiveRight = left/right primitive state (density, x velocity, y velocity, pressure)
! flux         = numerical flux (mass, x momentum, y momentum, energy)
! dflux        = ignore this
subroutine miczek_flux(normalX, normalY, primitiveLeft, primitiveRight, flux, dflux)
  use comvar
  implicit none

  real :: normalX, normalY, primitiveLeft(4), primitiveRight(4), flux(4), dflux(4)

  integer :: i,j
  real    :: rhoLeft, uLeft, vLeft, pressureLeft, enthalpyLeft, &
       rhoRight, uRight, vRight, pressureRight, enthalpyRight, &
       uAveraged, vAveraged, velSqrAveraged, soundspeedSqrAveraged, soundspeedAveraged, enthalpyAveraged, &
       sqrtRhoLeft, sqrtRhoRight, roeAverageFactor, &
       velNormalLeft, velNormalRight, velNormalAveraged, velTangentialAveraged, FluxAverage(4), upwinding(4), &
       m1, m2, a1, a2, a3, a4, l1, l2, l3, l4, &
       a1l1, a2l2, a3l3, a4l4, aact, aast, &
       diffRho, diffMomX, diffMomY, diffEnergy, &
       absVelNormalAveraged, &
       dUdV(0:3, 0:3), dVdU(0:3, 0:3), miczekMatrixPrim(0:3, 0:3), diff(0:3), res(0:3), res2(0:3), res3(0:3)
 
  rhoLeft = primitiveLeft(1)
  uLeft = primitiveLeft(2)
  vLeft = primitiveLeft(3)
  pressureLeft = primitiveLeft(4)
  enthalpyLeft = GAMMA*pressureLeft/rhoLeft/(GAMMA-1.0) + 0.5*(uLeft**2 + vLeft**2)

  rhoRight = primitiveRight(1)
  uRight = primitiveRight(2)
  vRight = primitiveRight(3)
  pressureRight = primitiveRight(4)
  enthalpyRight = GAMMA*pressureRight/rhoRight/(GAMMA-1.0) + 0.5*(uRight**2 + vRight**2)

  !     Rotated velocity
  velNormalLeft = uLeft*normalX + vLeft*normalY
  velNormalRight = uRight*normalX + vRight*normalY

  !     Average of flux functions
  FluxAverage(1) = rhoLeft*velNormalLeft                              + rhoRight*velNormalRight
  FluxAverage(2) = pressureLeft*normalX + rhoLeft*uLeft*velNormalLeft + pressureRight*normalX + rhoRight*uRight*velNormalRight
  FluxAverage(3) = pressureLeft*normalY + rhoLeft*vLeft*velNormalLeft + pressureRight*normalY + rhoRight*vRight*velNormalRight
  FluxAverage(4) = rhoLeft*enthalpyLeft*velNormalLeft                 + rhoRight*enthalpyRight*velNormalRight

  !     Roe average
  sqrtRhoLeft = 1.0 !sqrt(rhoLeft)
  sqrtRhoRight = 1.0 !sqrt(rhoRight)
  roeAverageFactor   = 0.5 !1.0/(sqrtRhoLeft + sqrtRhoRight)

  uAveraged          = (uLeft*sqrtRhoLeft + uRight*sqrtRhoRight)*roeAverageFactor
  vAveraged          = (vLeft*sqrtRhoLeft + vRight*sqrtRhoRight)*roeAverageFactor
  enthalpyAveraged   = (enthalpyLeft*sqrtRhoLeft + enthalpyRight*sqrtRhoRight)*roeAverageFactor
  velSqrAveraged       = uAveraged**2 + vAveraged**2
  soundspeedSqrAveraged  = (GAMMA-1.0)*(enthalpyAveraged - 0.5*velSqrAveraged)

  if(soundspeedSqrAveraged .le. 0.0)then
     print*,'Sonic speed is imaginary'
     print*,'Left/right conserved values'
     print*,primitiveLeft(:)
     print*,primitiveRight(:)
     print*
     print*,'Left/right primitive values'
     print*,rhoLeft,uLeft,vLeft,pressureLeft
     print*,rhoRight,uRight,vRight,pressureRight
     stop
  endif

  soundspeedAveraged  = sqrt(soundspeedSqrAveraged)
  velNormalAveraged = uAveraged*normalX + vAveraged*normalY
  velTangentialAveraged =-uAveraged*normalY + vAveraged*normalX

  absVelNormalAveraged = abs(velNormalAveraged)

  !     Miczek matrix (simplified form)
  miczekMatrixPrim = 0

  miczekMatrixPrim(0,0) = absVelNormalAveraged
  miczekMatrixPrim(1,1) = miczekMatrixPrim(0,0)
  miczekMatrixPrim(2,2) = miczekMatrixPrim(0,0)
  miczekMatrixPrim(3,3) = miczekMatrixPrim(0,0)

  miczekMatrixPrim(1,3) = normalX
  miczekMatrixPrim(3,1) = -soundspeedSqrAveraged*normalX

  miczekMatrixPrim(2,3) = normalY
  miczekMatrixPrim(3,2) = -soundspeedSqrAveraged*normalY

  !     Difference of conserved variables
  diff(0) = rhoRight                                 - rhoLeft
  diff(1) = rhoRight*uRight                          - rhoLeft*uLeft
  diff(2) = rhoRight*vRight                          - rhoLeft*vLeft
  diff(3) = (rhoRight*enthalpyRight - pressureRight) - (rhoLeft*enthalpyLeft - pressureLeft)

  !     Prim to Cons trafo for matrix
  dUdV = 0
  dVdU = 0

  dUdV(0,0) = 1.0
  dUdV(1,0) = uAveraged
  dUdV(1,1) = 1.0
  dUdV(2,0) = vAveraged
  dUdV(2,2) = 1.0
  dUdV(3,0) = velSqrAveraged*0.5
  dUdV(3,1) = uAveraged
  dUdV(3,2) = vAveraged
  dUdV(3,3) = 1.0/(GAMMA - 1)

  dVdU(0,0) = 1.0
  dVdU(1,0) = -uAveraged
  dVdU(1,1) = 1.0
  dVdU(2,0) = -vAveraged
  dVdU(2,2) = 1.0
  dVdU(3,0) = (GAMMA-1)*velSqrAveraged*0.5
  dVdU(3,1) = -(GAMMA-1)*uAveraged
  dVdU(3,2) = -(GAMMA-1)*vAveraged
  dVdU(3,3) = (GAMMA - 1)

  do i = 0, 3
     res(i) = 0
     do j = 0, 3
        res(i) = res(i) + dVdU(i, j)*diff(j)
     end do
  end do

  do i = 0, 3
     res2(i) = 0
     do j = 0, 3
        res2(i) = res2(i) + miczekMatrixPrim(i, j)*res(j)
     end do
  end do


  do i = 0, 3
     res3(i) = 0
     do j = 0, 3
        res3(i) = res3(i) + dUdV(i, j)*res2(j)
     end do
  end do

  upwinding(1) = res3(0)
  upwinding(2) = res3(1)
  upwinding(3) = res3(2)
  upwinding(4) = res3(3)

  !     Total flux
  flux    =  0.5*( FluxAverage - upwinding )
  dflux   = -0.5*upwinding

end subroutine miczek_flux
