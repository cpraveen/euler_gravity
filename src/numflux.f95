subroutine numflux(ct, st, priml, primr, flux, dflux)
   use comvar
   implicit none

   real :: ct, st
   real :: priml(4), primr(4), flux(4), dflux(4)

   if(fluxtype == iroe)then
      call roe_flux(ct, st, priml, primr, flux, dflux)
   else if(fluxtype == irusanov)then
      call rusanov_flux(ct, st, priml, primr, flux, dflux)
   else if(fluxtype == imiczek)then
      call miczek_flux(ct, st, priml, primr, flux, dflux)
   else
      write(*,*)'Uknown flux type fluxtype =', fluxtype
      stop
   endif

end subroutine numflux
