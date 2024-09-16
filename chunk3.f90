! melting of snow at +2 C
  do k=1,nlev

     do i=1,mgncol

        if (t(i,k)+tlat(i,k)/cpp*deltat > snowmelt) then
           if (dums(i,k) > 0._r8) then

              ! make sure melting snow doesn't reduce temperature below threshold
              dum = -xlf/cpp*dums(i,k)
              if (t(i,k)+tlat(i,k)/cpp*deltat+dum.lt. snowmelt) then
                 dum = (t(i,k)+tlat(i,k)/cpp*deltat-snowmelt)*cpp/xlf
                 dum = dum/dums(i,k)
                 dum = max(0._r8,dum)
                 dum = min(1._r8,dum)
              else
                 dum = 1._r8
              end if

              qstend(i,k)=qstend(i,k)-dum*dums(i,k)/deltat
              nstend(i,k)=nstend(i,k)-dum*dumns(i,k)/deltat
              qrtend(i,k)=qrtend(i,k)+dum*dums(i,k)/deltat
              nrtend(i,k)=nrtend(i,k)+dum*dumns(i,k)/deltat

              dum1=-xlf*dum*dums(i,k)/deltat
              tlat(i,k)=tlat(i,k)+dum1
              meltsdttot(i,k)=meltsdttot(i,k) + dum1
           end if
        end if
     enddo
  enddo

  ! melting of graupel at +2 C
  do k=1,nlev

     do i=1,mgncol

        if (t(i,k)+tlat(i,k)/cpp*deltat > snowmelt) then
           if (dumg(i,k) > 0._r8) then

              ! make sure melting graupel doesn't reduce temperature below threshold
              dum = -xlf/cpp*dumg(i,k)
              if (t(i,k)+tlat(i,k)/cpp*deltat+dum.lt. snowmelt) then
                 dum = (t(i,k)+tlat(i,k)/cpp*deltat-snowmelt)*cpp/xlf
                 dum = dum/dumg(i,k)
                 dum = max(0._r8,dum)
                 dum = min(1._r8,dum)
              else
                 dum = 1._r8
              end if


              qgtend(i,k)=qgtend(i,k)-dum*dumg(i,k)/deltat
              ngtend(i,k)=ngtend(i,k)-dum*dumng(i,k)/deltat
              qrtend(i,k)=qrtend(i,k)+dum*dumg(i,k)/deltat
              nrtend(i,k)=nrtend(i,k)+dum*dumng(i,k)/deltat

              dum1=-xlf*dum*dumg(i,k)/deltat
              tlat(i,k)=tlat(i,k)+dum1
              meltsdttot(i,k)=meltsdttot(i,k) + dum1
           end if
        end if
     enddo
  enddo

   do k=1,nlev
      do i=1,mgncol

        ! freezing of rain at -5 C

        if (t(i,k)+tlat(i,k)/cpp*deltat < rainfrze) then

           if (dumr(i,k) > 0._r8) then

              ! make sure freezing rain doesn't increase temperature above threshold
              dum = xlf/cpp*dumr(i,k)
              if (t(i,k)+tlat(i,k)/cpp*deltat+dum.gt.rainfrze) then
                 dum = -(t(i,k)+tlat(i,k)/cpp*deltat-rainfrze)*cpp/xlf
                 dum = dum/dumr(i,k)
                 dum = max(0._r8,dum)
                 dum = min(1._r8,dum)
              else
                 dum = 1._r8
              end if

              qrtend(i,k)=qrtend(i,k)-dum*dumr(i,k)/deltat
              nrtend(i,k)=nrtend(i,k)-dum*dumnr(i,k)/deltat

              ! get mean size of rain = 1/lamr, add frozen rain to either snow or cloud ice
              ! depending on mean rain size
              ! add to graupel if using that option....

              call size_dist_param_basic(mg_rain_props, dumr(i,k), dumnr(i,k), &
                   lamr(i,k))

              if (lamr(i,k) < 1._r8/Dcs) then

                 if(do_hail.or.do_graupel) then
                    qgtend(i,k)=qgtend(i,k)+dum*dumr(i,k)/deltat
                    ngtend(i,k)=ngtend(i,k)+dum*dumnr(i,k)/deltat
                 else
                    qstend(i,k)=qstend(i,k)+dum*dumr(i,k)/deltat
                    nstend(i,k)=nstend(i,k)+dum*dumnr(i,k)/deltat
                 end if

              else
                 qitend(i,k)=qitend(i,k)+dum*dumr(i,k)/deltat
                 nitend(i,k)=nitend(i,k)+dum*dumnr(i,k)/deltat
              end if

              ! heating tendency
              dum1 = xlf*dum*dumr(i,k)/deltat
              frzrdttot(i,k)=frzrdttot(i,k) + dum1
              tlat(i,k)=tlat(i,k)+dum1

           end if
        end if

      enddo
   enddo
