subroutine reconstruct(conjm1, conj, conjp1, conl)
   use comvar
   implicit none

   real    :: conjm1(4), conj(4), conjp1(4), conjp2(4)
   real    :: conl(4)

   integer :: i
   real    :: kkk=1.0/3.0
   real    :: minmod
   real    :: theta = 2.0

   ! reconstructed states
   if(limtype == ford)then
   ! first order
      do i=1,4
         conl(i) = conj(i)
      enddo
   else if(limtype == muscl3)then
   !muscl scheme
      do i=1,4
         conl(i) = conj(i)   + 0.25*( (1.0-kkk)*(conj(i) - conjm1(i)) &
                                    + (1.0+kkk)*(conjp1(i) - conj(i)) )
      enddo
   else if(limtype == mmod)then
   ! minmod limiter
      do i=1,4
         conl(i) = conj(i) + 0.5*minmod( theta*(conj(i)-conjm1(i)), &
                                          0.5*(conjp1(i)-conjm1(i)), &
                                          theta*(conjp1(i)-conj(i)) )
      enddo
   endif

end subroutine reconstruct
