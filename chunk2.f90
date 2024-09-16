do k=1,nlev
    do i=1,mgncol
        ! freezing of rain at -5 C

        if (t(i,k) < rainfrze) then

           if (qr(i,k) > 0._r8) then

              ! make sure freezing rain doesn't increase temperature above threshold
              dum = xlf/cpp*qr(i,k)
              if (t(i,k)+dum > rainfrze) then
                 dum = -(t(i,k)-rainfrze)*cpp/xlf
                 dum = dum/qr(i,k)
                 dum = max(0._r8,dum)
                 dum = min(1._r8,dum)
              else
                 dum = 1._r8
              end if

              minstrf(i,k) = dum*qr(i,k)
              ninstrf(i,k) = dum*nr(i,k)

              ! heating tendency
              dum1 = xlf*minstrf(i,k)/deltat
              tlat(i,k)=tlat(i,k)+dum1
              frzrdttot(i,k)=frzrdttot(i,k) + dum1


              qr(i,k) = max(qr(i,k) - minstrf(i,k), 0._r8)
              nr(i,k) = max(nr(i,k) - ninstrf(i,k), 0._r8)

              ! freeze rain to graupel not snow.
              if(do_hail.or.do_graupel) then
                 qg(i,k) = max(qg(i,k) + minstrf(i,k), 0._r8)
                 ng(i,k) = max(ng(i,k) + ninstrf(i,k), 0._r8)
              else
                 qs(i,k) = max(qs(i,k) + minstrf(i,k), 0._r8)
                 ns(i,k) = max(ns(i,k) + ninstrf(i,k), 0._r8)
              end if
           end if
        end if
     end do
  end do 

  do k=1,nlev
    do i=1,mgncol
        ! obtain in-cloud values of cloud water/ice mixing ratios and number concentrations
        !-------------------------------------------------------
        ! for microphysical process calculations
        ! units are kg/kg for mixing ratio, 1/kg for number conc

        if (qc(i,k).ge.qsmall) then
           ! limit in-cloud values to 0.005 kg/kg
           qcic(i,k)=min(qc(i,k)/lcldm(i,k),5.e-3_r8)
           ncic(i,k)=max(nc(i,k)/lcldm(i,k),0._r8)

           ! specify droplet concentration
           if (nccons) then
              ncic(i,k)=ncnst/rho(i,k)
           end if
        else
           qcic(i,k)=0._r8
           ncic(i,k)=0._r8
        end if

        if (qi(i,k).ge.qsmall) then
           ! limit in-cloud values to 0.005 kg/kg
           qiic(i,k)=min(qi(i,k)/icldm(i,k),5.e-3_r8)
           niic(i,k)=max(ni(i,k)/icldm(i,k),0._r8)

           ! switch for specification of cloud ice number
           if (nicons) then
              niic(i,k)=ninst/rho(i,k)
           end if
        else
           qiic(i,k)=0._r8
           niic(i,k)=0._r8
        end if

     end do
  end do

  !========================================================================

  ! for sub-columns cldm has already been set to 1 if cloud
  ! water or ice is present, so precip_frac will be correctly set below
  ! and nothing extra needs to be done here

  precip_frac = cldm

  micro_vert_loop: do k=1,nlev

     if (trim(micro_mg_precip_frac_method) == 'in_cloud') then

        if (k /= 1) then
           where (qc(:,k) < qsmall .and. qi(:,k) < qsmall)
              precip_frac(:,k) = precip_frac(:,k-1)
           end where
        endif

     else if (trim(micro_mg_precip_frac_method) == 'max_overlap') then

        ! calculate precip fraction based on maximum overlap assumption

        ! if rain or snow mix ratios are smaller than threshold,
        ! then leave precip_frac as cloud fraction at current level
        if (k /= 1) then
           where (qr(:,k-1) >= qsmall .or. qs(:,k-1) >= qsmall)
              precip_frac(:,k)=max(precip_frac(:,k-1),precip_frac(:,k))
           end where
        end if

     endif


     !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     ! get size distribution parameters based on in-cloud cloud water
     ! these calculations also ensure consistency between number and mixing ratio
     !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

     ! cloud liquid
     !-------------------------------------------

     call size_dist_param_liq(mg_liq_props, qcic(1:mgncol,k), ncic(1:mgncol,k),& 
          rho(1:mgncol,k), pgam(1:mgncol,k), lamc(1:mgncol,k), mgncol)


     !========================================================================
     ! autoconversion of cloud liquid water to rain
     ! formula from Khrouditnov and Kogan (2000), modified for sub-grid distribution of qc
     ! minimum qc of 1 x 10^-8 prevents floating point error

     if (.not. do_sb_physics) then
       call kk2000_liq_autoconversion(microp_uniform, qcic(1:mgncol,k), &
          ncic(:,k), rho(:,k), relvar(:,k), prc(:,k), nprc(:,k), nprc1(:,k), mgncol)
     endif

     ! assign qric based on prognostic qr, using assumed precip fraction
     ! note: this could be moved above for consistency with qcic and qiic calculations
     qric(:,k) = qr(:,k)/precip_frac(:,k)
     nric(:,k) = nr(:,k)/precip_frac(:,k)

     ! limit in-precip mixing ratios to 10 g/kg
     qric(:,k)=min(qric(:,k),0.01_r8)

     ! add autoconversion to precip from above to get provisional rain mixing ratio
     ! and number concentration (qric and nric)

     where (qric(:,k).lt.qsmall)
        qric(:,k)=0._r8
        nric(:,k)=0._r8
     end where

     ! make sure number concentration is a positive number to avoid
     ! taking root of negative later

     nric(:,k)=max(nric(:,k),0._r8)

     ! Get size distribution parameters for cloud ice

     call size_dist_param_basic(mg_ice_props, qiic(:,k), niic(:,k), &
          lami(:,k), mgncol, n0=n0i(:,k))
	  
     ! Alternative autoconversion 
     if (do_sb_physics) then
       call sb2001v2_liq_autoconversion(pgam(:,k),qcic(:,k),ncic(:,k), &
            qric(:,k),rho(:,k),relvar(:,k),prc(:,k),nprc(:,k),nprc1(:,k), mgncol)     
     endif	  

     !.......................................................................
     ! Autoconversion of cloud ice to snow
     ! similar to Ferrier (1994)

     if (do_cldice) then
        call ice_autoconversion(t(:,k), qiic(:,k), lami(:,k), n0i(:,k), &
             dcs, prci(:,k), nprci(:,k), mgncol)
     else
        ! Add in the particles that we have already converted to snow, and
        ! don't do any further autoconversion of ice.
        prci(:,k)  = tnd_qsnow(:,k) / cldm(:,k)
        nprci(:,k) = tnd_nsnow(:,k) / cldm(:,k)
     end if

     ! note, currently we don't have this
     ! inside the do_cldice block, should be changed later
     ! assign qsic based on prognostic qs, using assumed precip fraction
     qsic(:,k) = qs(:,k)/precip_frac(:,k)
     nsic(:,k) = ns(:,k)/precip_frac(:,k)

     ! limit in-precip mixing ratios to 10 g/kg
     qsic(:,k)=min(qsic(:,k),0.01_r8)

     ! if precip mix ratio is zero so should number concentration

     where (qsic(:,k) < qsmall)
        qsic(:,k)=0._r8
        nsic(:,k)=0._r8
     end where

     ! make sure number concentration is a positive number to avoid
     ! taking root of negative later

     nsic(:,k)=max(nsic(:,k),0._r8)


     ! also do this for graupel, which is assumed to be 'precip_frac'
     qgic(:,k) = qg(:,k)/precip_frac(:,k)
     ngic(:,k) = ng(:,k)/precip_frac(:,k)

     ! limit in-precip mixing ratios to 10 g/kg
     qgic(:,k)=min(qgic(:,k),0.01_r8)

     ! if precip mix ratio is zero so should number concentration
     where (qgic(:,k) < qsmall)
        qgic(:,k)=0._r8
        ngic(:,k)=0._r8
     end where

     ! make sure number concentration is a positive number to avoid
     ! taking root of negative later

     ngic(:,k)=max(ngic(:,k),0._r8)    

     !.......................................................................
     ! get size distribution parameters for precip
     !......................................................................
     ! rain

     call size_dist_param_basic(mg_rain_props, qric(:,k), nric(:,k), &
          lamr(:,k), mgncol, n0=n0r(:,k))

     where (lamr(:,k) >= qsmall)

        ! provisional rain number and mass weighted mean fallspeed (m/s)

        unr(:,k) = min(arn(:,k)*gamma_br_plus1/lamr(:,k)**br,9.1_r8*rhof(:,k))
        umr(:,k) = min(arn(:,k)*gamma_br_plus4/(6._r8*lamr(:,k)**br),9.1_r8*rhof(:,k))

     elsewhere
        umr(:,k) = 0._r8
        unr(:,k) = 0._r8
     end where

     !......................................................................
     ! snow

     call size_dist_param_basic(mg_snow_props, qsic(:,k), nsic(:,k), &
          lams(:,k), mgncol, n0=n0s(:,k))

     where (lams(:,k) > 0._r8)

        ! provisional snow number and mass weighted mean fallspeed (m/s)

        ums(:,k) = min(asn(:,k)*gamma_bs_plus4/(6._r8*lams(:,k)**bs),1.2_r8*rhof(:,k))
        uns(:,k) = min(asn(:,k)*gamma_bs_plus1/lams(:,k)**bs,1.2_r8*rhof(:,k))

     elsewhere
        ums(:,k) = 0._r8
        uns(:,k) = 0._r8
     end where

     !......................................................................
     !       graupel/hail density set (Hail = 400, Graupel = 500 from M2005)
        
     if (do_hail) then 
        bgtmp = bh 
        rhogtmp = rhoh
     end if
     if (do_graupel) then 
        bgtmp = bg
        rhogtmp = rhog
     end if

     !  graupel/hail size distributions and properties

     if (do_hail) then
        call size_dist_param_basic(mg_hail_props, qgic(:,k), ngic(:,k), &
          lamg(:,k), mgncol, n0=n0g(:,k))
     end if
     if (do_graupel) then
        call size_dist_param_basic(mg_graupel_props, qgic(:,k), ngic(:,k), &
          lamg(:,k), mgncol, n0=n0g(:,k))
     end if
        
     where (lamg(:,k) > 0._r8)

        ! provisional graupel/hail number and mass weighted mean fallspeed (m/s)
        umg(:,k) = min(agn(:,k)*gamma(4._r8+bgtmp)/(6._r8*lamg(:,k)**bgtmp),20._r8*rhof(:,k))
        ung(:,k) = min(agn(:,k)*gamma(1._r8+bgtmp)/lamg(:,k)**bgtmp,20._r8*rhof(:,k))

     elsewhere
        umg(:,k) = 0._r8
        ung(:,k) = 0._r8
     end where

     if (do_cldice) then
        if (.not. use_hetfrz_classnuc) then

           ! heterogeneous freezing of cloud water
           !----------------------------------------------

           call immersion_freezing(microp_uniform, t(:,k), pgam(:,k), lamc(:,k), &
                qcic(1:mgncol,k), ncic(:,k), relvar(:,k), mnuccc(:,k), nnuccc(:,k), mgncol)

           ! make sure number of droplets frozen does not exceed available ice nuclei concentration
           ! this prevents 'runaway' droplet freezing

           where (qcic(1:mgncol,k).ge.qsmall .and. t(:,k).lt.269.15_r8)
              where (nnuccc(:,k)*lcldm(:,k).gt.nnuccd(:,k))
                 ! scale mixing ratio of droplet freezing with limit
                 mnuccc(:,k)=mnuccc(:,k)*(nnuccd(:,k)/(nnuccc(:,k)*lcldm(:,k)))
                 nnuccc(:,k)=nnuccd(:,k)/lcldm(:,k)
              end where
           end where

           mdust = size(rndst,3)
           call contact_freezing(microp_uniform, t(:,k), p(:,k), rndst(:,k,:), &
                nacon(:,k,:), pgam(:,k), lamc(:,k), qcic(1:mgncol,k), ncic(:,k), &
                relvar(:,k), mnucct(:,k), nnucct(:,k), mgncol, mdust)

           mnudep(:,k)=0._r8
           nnudep(:,k)=0._r8

        else

           ! Mass of droplets frozen is the average droplet mass, except
           ! with two limiters: concentration must be at least 1/cm^3, and
           ! mass must be at least the minimum defined above.
           mi0l = qcic(:,k)/max(ncic(:,k), 1.0e6_r8/rho(:,k))
           mi0l = max(mi0l_min, mi0l)

           where (qcic(:,k) >= qsmall)
              nnuccc(:,k) = frzimm(:,k)*1.0e6_r8/rho(:,k)
              mnuccc(:,k) = nnuccc(:,k)*mi0l

              nnucct(:,k) = frzcnt(:,k)*1.0e6_r8/rho(:,k)
              mnucct(:,k) = nnucct(:,k)*mi0l

              nnudep(:,k) = frzdep(:,k)*1.0e6_r8/rho(:,k)
              mnudep(:,k) = nnudep(:,k)*mi0
           elsewhere
              nnuccc(:,k) = 0._r8
              mnuccc(:,k) = 0._r8

              nnucct(:,k) = 0._r8
              mnucct(:,k) = 0._r8

              nnudep(:,k) = 0._r8
              mnudep(:,k) = 0._r8
           end where

        end if

     else
        mnuccc(:,k)=0._r8
        nnuccc(:,k)=0._r8
        mnucct(:,k)=0._r8
        nnucct(:,k)=0._r8
        mnudep(:,k)=0._r8
        nnudep(:,k)=0._r8
     end if

     call snow_self_aggregation(t(:,k), rho(:,k), asn(:,k), rhosn, qsic(:,k), nsic(:,k), &
          nsagg(:,k), mgncol)

     call accrete_cloud_water_snow(t(:,k), rho(:,k), asn(:,k), uns(:,k), mu(:,k), &
          qcic(1:mgncol,k), ncic(:,k), qsic(:,k), pgam(:,k), lamc(:,k), lams(:,k), n0s(:,k), &
          psacws(:,k), npsacws(:,k), mgncol)

     if (do_cldice) then
        call secondary_ice_production(t(:,k), psacws(:,k), msacwi(:,k), nsacwi(:,k), mgncol)
     else
        nsacwi(:,k) = 0.0_r8
        msacwi(:,k) = 0.0_r8
     end if

     call accrete_rain_snow(t(:,k), rho(:,k), umr(:,k), ums(:,k), unr(:,k), uns(:,k), &
          qric(:,k), qsic(:,k), lamr(:,k), n0r(:,k), lams(:,k), n0s(:,k), &
          pracs(:,k), npracs(:,k), mgncol)

     call heterogeneous_rain_freezing(t(:,k), qric(:,k), nric(:,k), lamr(:,k), &
          mnuccr(:,k), nnuccr(:,k), mgncol)

     if (do_sb_physics) then
       call sb2001v2_accre_cld_water_rain(qcic(1:mgncol,k), ncic(:,k), qric(:,k), &
            rho(:,k), relvar(:,k), pra(:,k), npra(:,k), mgncol)     
     else
       call accrete_cloud_water_rain(microp_uniform, qric(:,k), qcic(1:mgncol,k), &
            ncic(:,k), relvar(:,k), accre_enhan(:,k), pra(:,k), npra(:,k), mgncol)
     endif

     call self_collection_rain(rho(:,k), qric(:,k), nric(:,k), nragg(:,k), mgncol)

     if (do_cldice) then
        call accrete_cloud_ice_snow(t(:,k), rho(:,k), asn(:,k), qiic(:,k), niic(:,k), &
             qsic(:,k), lams(:,k), n0s(:,k), prai(:,k), nprai(:,k), mgncol)
     else
        prai(:,k) = 0._r8
        nprai(:,k) = 0._r8
     end if

     call bergeron_process_snow(t(:,k), rho(:,k), dv(:,k), mu(:,k), sc(:,k), &
          qvl(:,k), qvi(:,k), asn(:,k), qcic(1:mgncol,k), qsic(:,k), lams(:,k), n0s(:,k), &
          bergs(:,k), mgncol)

     bergs(:,k)=bergs(:,k)*micro_mg_berg_eff_factor

     if (do_cldice) then

        call ice_deposition_sublimation(t(:,k), q(:,k), qi(:,k), ni(:,k), &
             icldm(:,k), rho(:,k), dv(:,k), qvl(:,k), qvi(:,k), &
             berg(:,k), vap_dep(:,k), ice_sublim(:,k), mgncol)

        berg(:,k)=berg(:,k)*micro_mg_berg_eff_factor

        where (ice_sublim(:,k) < 0._r8 .and. qi(:,k) > qsmall .and. icldm(:,k) > mincld)
           nsubi(:,k) = sublim_factor*ice_sublim(:,k) / qi(:,k) * ni(:,k) / icldm(:,k)
        elsewhere
           nsubi(:,k) = 0._r8
        end where

        ! bergeron process should not reduce nc unless
        ! all ql is removed (which is handled elsewhere)
        !in fact, nothing in this entire file makes nsubc nonzero.
        nsubc(:,k) = 0._r8

     end if !do_cldice

! Process rate calls for graupel   
!===================================================================

     if(do_hail.or.do_graupel) then
        call graupel_collecting_snow(qsic(:,k),qric(:,k),umr(:,k),ums(:,k), &
             rho(:,k),lamr(:,k),n0r(:,k),lams(:,k),n0s(:,k), &
             psacr(:,k), mgncol)

       call graupel_collecting_cld_water(qgic(:,k),qcic(:,k),ncic(:,k),rho(:,k), &
             n0g(:,k),lamg(:,k),bgtmp,agn(:,k), psacwg(:,k), npsacwg(:,k), mgncol)
        
        call graupel_riming_liquid_snow(psacws(:,k),qsic(:,k),qcic(:,k),nsic(:,k), &
             rho(:,k),rhosn,rhogtmp,asn(:,k),lams(:,k),n0s(:,k),deltat, &
             pgsacw(:,k),nscng(:,k),mgncol)

        call graupel_collecting_rain(qric(:,k),qgic(:,k),umg(:,k), &
             umr(:,k),ung(:,k),unr(:,k),rho(:,k),n0r(:,k),lamr(:,k),n0g(:,k), &
             lamg(:,k), pracg(:,k),npracg(:,k),mgncol)
        
!AG note: Graupel rain riming snow changes  
!    pracs, npracs, (accretion of rain by snow)  psacr (collection of snow by rain)

       call graupel_rain_riming_snow(pracs(:,k),npracs(:,k),psacr(:,k),qsic(:,k), &
             qric(:,k),nric(:,k),nsic(:,k),n0s(:,k),lams(:,k),n0r(:,k),lamr(:,k), &
             deltat,pgracs(:,k),ngracs(:,k),mgncol)
        
        call graupel_rime_splintering(t(:,k),qcic(:,k),qric(:,k),qgic(:,k), &
             psacwg(:,k),pracg(:,k),qmultg(:,k),nmultg(:,k),qmultrg(:,k), &
             nmultrg(:,k),mgncol)


        call evaporate_sublimate_precip_graupel(t(:,k), rho(:,k), &
             dv(:,k), mu(:,k), sc(:,k), q(:,k), qvl(:,k), qvi(:,k), &
             lcldm(:,k), precip_frac(:,k), arn(:,k), asn(:,k), agn(:,k), bgtmp, &
             qcic(1:mgncol,k), qiic(:,k), qric(:,k), qsic(:,k), qgic(:,k), &
             lamr(:,k), n0r(:,k), lams(:,k), n0s(:,k), lamg(:,k), n0g(:,k), &
             pre(:,k), prds(:,k), prdg(:,k), am_evp_st(:,k), mgncol)   
     else
             
        ! Routine without Graupel (original)        
        call evaporate_sublimate_precip(t(:,k), rho(:,k), &
          dv(:,k), mu(:,k), sc(:,k), q(:,k), qvl(:,k), qvi(:,k), &
          lcldm(:,k), precip_frac(:,k), arn(:,k), asn(:,k), qcic(1:mgncol,k), qiic(:,k), &
          qric(:,k), qsic(:,k), lamr(:,k), n0r(:,k), lams(:,k), n0s(:,k), &
          pre(:,k), prds(:,k), am_evp_st(:,k), mgncol)


     end if ! end do_graupel/hail loop

     do i=1,mgncol

        ! conservation to ensure no negative values of cloud water/precipitation
        ! in case microphysical process rates are large
        !===================================================================

        ! note: for check on conservation, processes are multiplied by omsm
        ! to prevent problems due to round off error

        ! conservation of qc
        !-------------------------------------------------------------------

        dum = ((prc(i,k)+pra(i,k)+mnuccc(i,k)+mnucct(i,k)+msacwi(i,k)+ &
             psacws(i,k)+bergs(i,k)+qmultg(i,k)+psacwg(i,k)+pgsacw(i,k))*lcldm(i,k)+ &
             berg(i,k))*deltat 

        if (dum.gt.qc(i,k)) then
                
           ratio = qc(i,k)/deltat/((prc(i,k)+pra(i,k)+mnuccc(i,k)+mnucct(i,k)+ &
                msacwi(i,k)+psacws(i,k)+bergs(i,k)+qmultg(i,k)+psacwg(i,k)+pgsacw(i,k))*lcldm(i,k)+&
                berg(i,k))*omsm

           qmultg(i,k)=qmultg(i,k)*ratio
           psacwg(i,k)=psacwg(i,k)*ratio
           pgsacw(i,k)=pgsacw(i,k)*ratio

           prc(i,k) = prc(i,k)*ratio
           pra(i,k) = pra(i,k)*ratio
           mnuccc(i,k) = mnuccc(i,k)*ratio
           mnucct(i,k) = mnucct(i,k)*ratio
           msacwi(i,k) = msacwi(i,k)*ratio
           psacws(i,k) = psacws(i,k)*ratio
           bergs(i,k) = bergs(i,k)*ratio
           berg(i,k) = berg(i,k)*ratio
           qcrat(i,k) = ratio
        else
           qcrat(i,k) = 1._r8
        end if

        !PMC 12/3/12: ratio is also frac of step w/ liquid.
        !thus we apply berg for "ratio" of timestep and vapor
        !deposition for the remaining frac of the timestep.
        if (qc(i,k) >= qsmall) then
           vap_dep(i,k) = vap_dep(i,k)*(1._r8-qcrat(i,k))
        end if

     end do

     do i=1,mgncol

        !=================================================================
        ! apply limiter to ensure that ice/snow sublimation and rain evap
        ! don't push conditions into supersaturation, and ice deposition/nucleation don't
        ! push conditions into sub-saturation
        ! note this is done after qc conservation since we don't know how large
        ! vap_dep is before then
        ! estimates are only approximate since other process terms haven't been limited
        ! for conservation yet

        ! first limit ice deposition/nucleation vap_dep + mnuccd
        dum1 = vap_dep(i,k) + mnuccd(i,k)
        if (dum1 > 1.e-20_r8) then
           dum = (q(i,k)-qvi(i,k))/(1._r8 + xxls_squared*qvi(i,k)/(cpp*rv*t(i,k)**2))/deltat
           dum = max(dum,0._r8)
           if (dum1 > dum) then
              ! Allocate the limited "dum" tendency to mnuccd and vap_dep
              ! processes. Don't divide by cloud fraction; these are grid-
              ! mean rates.
              dum1 = mnuccd(i,k) / (vap_dep(i,k)+mnuccd(i,k))
              mnuccd(i,k) = dum*dum1
              vap_dep(i,k) = dum - mnuccd(i,k)
           end if
        end if

     end do

     do i=1,mgncol

        !===================================================================
        ! conservation of nc
        !-------------------------------------------------------------------

        dum = (nprc1(i,k)+npra(i,k)+nnuccc(i,k)+nnucct(i,k)+ &
             npsacws(i,k)-nsubc(i,k)+npsacwg(i,k))*lcldm(i,k)*deltat

        if (dum.gt.nc(i,k)) then
           ratio = nc(i,k)/deltat/((nprc1(i,k)+npra(i,k)+nnuccc(i,k)+nnucct(i,k)+&
                npsacws(i,k)-nsubc(i,k)+npsacwg(i,k))*lcldm(i,k))*omsm
           npsacwg(i,k)=npsacwg(i,k)*ratio

           nprc1(i,k) = nprc1(i,k)*ratio
           npra(i,k) = npra(i,k)*ratio
           nnuccc(i,k) = nnuccc(i,k)*ratio
           nnucct(i,k) = nnucct(i,k)*ratio
           npsacws(i,k) = npsacws(i,k)*ratio
           nsubc(i,k)=nsubc(i,k)*ratio
        end if

        mnuccri(i,k)=0._r8
        nnuccri(i,k)=0._r8

        if (do_cldice) then

           ! freezing of rain to produce ice if mean rain size is smaller than Dcs
           if (lamr(i,k) > qsmall .and. 1._r8/lamr(i,k) < Dcs) then
              mnuccri(i,k)=mnuccr(i,k)
              nnuccri(i,k)=nnuccr(i,k)
              mnuccr(i,k)=0._r8
              nnuccr(i,k)=0._r8
           end if
        end if

     end do

     do i=1,mgncol

        ! conservation of rain mixing ratio
        !-------------------------------------------------------------------
        dum = ((-pre(i,k)+pracs(i,k)+mnuccr(i,k)+mnuccri(i,k) &
             +qmultrg(i,k)+pracg(i,k)+pgracs(i,k))*precip_frac(i,k)- &
             (pra(i,k)+prc(i,k))*lcldm(i,k))*deltat

        ! note that qrtend is included below because of instantaneous freezing/melt
        if (dum.gt.qr(i,k).and. &
             (-pre(i,k)+pracs(i,k)+mnuccr(i,k)+mnuccri(i,k)+qmultrg(i,k)+pracg(i,k)+pgracs(i,k)).ge.qsmall) then

           ratio = (qr(i,k)/deltat+(pra(i,k)+prc(i,k))*lcldm(i,k))/   &
                precip_frac(i,k)/(-pre(i,k)+pracs(i,k)+mnuccr(i,k)+mnuccri(i,k) &
                +qmultrg(i,k)+pracg(i,k)+pgracs(i,k))*omsm

           qmultrg(i,k)= qmultrg(i,k)*ratio
           pracg(i,k)=pracg(i,k)*ratio
           pgracs(i,k)=pgracs(i,k)*ratio

           pre(i,k)=pre(i,k)*ratio
           pracs(i,k)=pracs(i,k)*ratio
           mnuccr(i,k)=mnuccr(i,k)*ratio
           mnuccri(i,k)=mnuccri(i,k)*ratio
        end if

     end do

     do i=1,mgncol

        ! conservation of rain number
        !-------------------------------------------------------------------

        ! Add evaporation of rain number.
        if (pre(i,k) < 0._r8) then
           nsubr(i,k) = pre(i,k)*nr(i,k)/qr(i,k)
        else
           nsubr(i,k) = 0._r8
        end if

     end do

     do i=1,mgncol

        dum = ((-nsubr(i,k)+npracs(i,k)+nnuccr(i,k)+nnuccri(i,k)-nragg(i,k)+npracg(i,k)+ngracs(i,k)) &
             *precip_frac(i,k)- nprc(i,k)*lcldm(i,k))*deltat

        if (dum.gt.nr(i,k)) then

           ratio = (nr(i,k)/deltat+nprc(i,k)*lcldm(i,k))/precip_frac(i,k)/ &
                (-nsubr(i,k)+npracs(i,k)+nnuccr(i,k)+nnuccri(i,k)-nragg(i,k)+npracg(i,k)+ngracs(i,k))*omsm

           npracg(i,k)=npracg(i,k)*ratio
           ngracs(i,k)=ngracs(i,k)*ratio
           nragg(i,k)=nragg(i,k)*ratio
           npracs(i,k)=npracs(i,k)*ratio
           nnuccr(i,k)=nnuccr(i,k)*ratio
           nsubr(i,k)=nsubr(i,k)*ratio
           nnuccri(i,k)=nnuccri(i,k)*ratio
        end if

     end do

     if (do_cldice) then

        do i=1,mgncol

           ! conservation of qi
           !-------------------------------------------------------------------

           dum = ((-mnuccc(i,k)-mnucct(i,k)-mnudep(i,k)-msacwi(i,k)-qmultg(i,k))*lcldm(i,k)+(prci(i,k)+ &
                prai(i,k))*icldm(i,k)+(-qmultrg(i,k)-mnuccri(i,k))*precip_frac(i,k) &
                -ice_sublim(i,k)-vap_dep(i,k)-berg(i,k)-mnuccd(i,k))*deltat

           if (dum.gt.qi(i,k)) then

              ratio = (qi(i,k)/deltat+vap_dep(i,k)+berg(i,k)+mnuccd(i,k)+ &
                   (mnuccc(i,k)+mnucct(i,k)+mnudep(i,k)+msacwi(i,k)+qmultg(i,k))*lcldm(i,k)+ &
                   (qmultrg(i,k)+mnuccri(i,k))*precip_frac(i,k))/ &
                   ((prci(i,k)+prai(i,k))*icldm(i,k)-ice_sublim(i,k))*omsm

              prci(i,k) = prci(i,k)*ratio
              prai(i,k) = prai(i,k)*ratio
              ice_sublim(i,k) = ice_sublim(i,k)*ratio
           end if

        end do

     end if

     if (do_cldice) then

        do i=1,mgncol

           ! conservation of ni
           !-------------------------------------------------------------------
           if (use_hetfrz_classnuc) then
              tmpfrz = nnuccc(i,k)
           else
              tmpfrz = 0._r8
           end if

           dum = ((-nnucct(i,k)-tmpfrz-nnudep(i,k)-nsacwi(i,k)-nmultg(i,k))*lcldm(i,k)+(nprci(i,k)+ &
                nprai(i,k)-nsubi(i,k))*icldm(i,k)+(-nmultrg(i,k)-nnuccri(i,k))*precip_frac(i,k)- &
                nnuccd(i,k))*deltat

           if (dum.gt.ni(i,k)) then

              ratio = (ni(i,k)/deltat+nnuccd(i,k)+ &
                 (nnucct(i,k)+tmpfrz+nnudep(i,k)+nsacwi(i,k)+nmultg(i,k))*lcldm(i,k)+ &
                 (nnuccri(i,k)+nmultrg(i,k))*precip_frac(i,k))/ &
                 ((nprci(i,k)+nprai(i,k)-nsubi(i,k))*icldm(i,k))*omsm

              nprci(i,k) = nprci(i,k)*ratio
              nprai(i,k) = nprai(i,k)*ratio
              nsubi(i,k) = nsubi(i,k)*ratio
           end if

        end do

     end if

     do i=1,mgncol

        ! conservation of snow mixing ratio
        !-------------------------------------------------------------------

        if (do_hail .or. do_graupel) then
        ! NOTE: mnuccr is moved to graupel when active
        ! psacr is a positive value, but a loss for snow
        !HM: psacr is positive in dum (two negatives)

           dum = (-(prds(i,k)+pracs(i,k)-psacr(i,k))*precip_frac(i,k)-(prai(i,k)+prci(i,k))*icldm(i,k) &
             -(bergs(i,k)+psacws(i,k))*lcldm(i,k))*deltat
        else
           dum = (-(prds(i,k)+pracs(i,k)+mnuccr(i,k))*precip_frac(i,k)-(prai(i,k)+prci(i,k))*icldm(i,k) &
             -(bergs(i,k)+psacws(i,k))*lcldm(i,k))*deltat
        end if 

        if (dum.gt.qs(i,k).and.(psacr(i,k)-prds(i,k)).ge.qsmall) then

           if (do_hail .or. do_graupel) then        
          
              ratio = (qs(i,k)/deltat+(prai(i,k)+prci(i,k))*icldm(i,k)+ &
                   (bergs(i,k)+psacws(i,k))*lcldm(i,k)+pracs(i,k)*precip_frac(i,k))/ &
                   precip_frac(i,k)/(psacr(i,k)-prds(i,k))*omsm
            
              psacr(i,k)=psacr(i,k)*ratio
           else 
            
              ratio = (qs(i,k)/deltat+(prai(i,k)+prci(i,k))*icldm(i,k)+ &
                   (bergs(i,k)+psacws(i,k))*lcldm(i,k)+(pracs(i,k)+mnuccr(i,k))*precip_frac(i,k))/ &
                   precip_frac(i,k)/(-prds(i,k))*omsm
           end if

           prds(i,k)=prds(i,k)*ratio
        end if

     end do

     do i=1,mgncol

        ! conservation of snow number
        !-------------------------------------------------------------------
        ! calculate loss of number due to sublimation
        ! for now neglect sublimation of ns
        nsubs(i,k)=0._r8

        if (do_hail .or. do_graupel) then        
           dum = ((-nsagg(i,k)-nsubs(i,k)+ngracs(i,k))*precip_frac(i,k)-nprci(i,k)*icldm(i,k)-nscng(i,k)*lcldm(i,k))*deltat
        else
           dum = ((-nsagg(i,k)-nsubs(i,k)-nnuccr(i,k))*precip_frac(i,k)-nprci(i,k)*icldm(i,k))*deltat
        end if

        if (dum.gt.ns(i,k)) then
           
           if (do_hail .or. do_graupel) then        
              ratio = (ns(i,k)/deltat+nprci(i,k)*icldm(i,k))/precip_frac(i,k)/ &
                   (-nsubs(i,k)-nsagg(i,k)+ngracs(i,k)+lcldm(i,k)/precip_frac(i,k)*nscng(i,k))*omsm
              nscng(i,k)=nscng(i,k)*ratio
              ngracs(i,k)=ngracs(i,k)*ratio
           else
              ratio = (ns(i,k)/deltat+nnuccr(i,k)* &
                   precip_frac(i,k)+nprci(i,k)*icldm(i,k))/precip_frac(i,k)/ &
                   (-nsubs(i,k)-nsagg(i,k))*omsm
           endif

           nsubs(i,k)=nsubs(i,k)*ratio
           nsagg(i,k)=nsagg(i,k)*ratio

        end if

     end do

! Graupel Conservation Checks
!-------------------------------------------------------------------

     if(do_hail.or.do_graupel) then

        ! conservation of graupel mass
        !-------------------------------------------------------------------
        do i=1,mgncol

           dum= ((-pracg(i,k)-pgracs(i,k)-prdg(i,k)-psacr(i,k)-mnuccr(i,k))*precip_frac(i,k) &
                + (-psacwg(i,k)-pgsacw(i,k))*lcldm(i,k))*deltat
           
           if (dum.gt.qg(i,k)) then
                   
              ! note: prdg is always negative (like prds), so it needs to be subtracted in ratio
              ratio = (qg(i,k)/deltat + (pracg(i,k)+pgracs(i,k)+psacr(i,k)+mnuccr(i,k))*precip_frac(i,k) &
                   + (psacwg(i,k)+pgsacw(i,k))*lcldm(i,k)) / ((-prdg(i,k))*precip_frac(i,k))  *omsm

              prdg(i,k)= prdg(i,k)*ratio

           end if

        end do

        ! conservation of graupel number: not needed, no sinks
        !-------------------------------------------------------------------
     end if

     do i=1,mgncol

        ! next limit ice and snow sublimation and rain evaporation
        ! get estimate of q and t at end of time step
        ! don't include other microphysical processes since they haven't
        ! been limited via conservation checks yet

        if ((pre(i,k)+prds(i,k)+prdg(i,k))*precip_frac(i,k)+ice_sublim(i,k) < -1.e-20_r8) then
           qtmp=q(i,k)-(ice_sublim(i,k)+vap_dep(i,k)+mnuccd(i,k)+ &
                (pre(i,k)+prds(i,k)+prdg(i,k))*precip_frac(i,k))*deltat
           ttmp=t(i,k)+((pre(i,k)*precip_frac(i,k))*xxlv+ &
                ((prds(i,k)+prdg(i,k))*precip_frac(i,k)+vap_dep(i,k)+ice_sublim(i,k)+mnuccd(i,k))*xxls)*deltat/cpp

           ! use rhw to allow ice supersaturation
           call qsat_water(ttmp, p(i,k), esn, qvn)

           ! modify ice/precip evaporation rate if q > qsat
           if (qtmp > qvn) then

              dum1=pre(i,k)*precip_frac(i,k)/((pre(i,k)+prds(i,k)+prdg(i,k))*precip_frac(i,k)+ice_sublim(i,k))
              dum2=prds(i,k)*precip_frac(i,k)/((pre(i,k)+prds(i,k)+prdg(i,k))*precip_frac(i,k)+ice_sublim(i,k))
              dum3=prdg(i,k)*precip_frac(i,k)/((pre(i,k)+prds(i,k)+prdg(i,k))*precip_frac(i,k)+ice_sublim(i,k))

              ! recalculate q and t after vap_dep and mnuccd but without evap or sublim
              qtmp=q(i,k)-(vap_dep(i,k)+mnuccd(i,k))*deltat
              ttmp=t(i,k)+((vap_dep(i,k)+mnuccd(i,k))*xxls)*deltat/cpp

              ! use rhw to allow ice supersaturation
              call qsat_water(ttmp, p(i,k), esn, qvn)

              dum=(qtmp-qvn)/(1._r8 + xxlv_squared*qvn/(cpp*rv*ttmp**2))
              dum=min(dum,0._r8)

              ! modify rates if needed, divide by precip_frac to get local (in-precip) value
              pre(i,k)=dum*dum1/deltat/precip_frac(i,k)

              ! do separately using RHI for prds and ice_sublim
              call qsat_ice(ttmp, p(i,k), esn, qvn)

              dum=(qtmp-qvn)/(1._r8 + xxls_squared*qvn/(cpp*rv*ttmp**2))
              dum=min(dum,0._r8)

              ! modify rates if needed, divide by precip_frac to get local (in-precip) value
              prds(i,k) = dum*dum2/deltat/precip_frac(i,k)
              prdg(i,k) = dum*dum3/deltat/precip_frac(i,k)

              ! don't divide ice_sublim by cloud fraction since it is grid-averaged
              dum1 = (1._r8-dum1-dum2-dum3)
              ice_sublim(i,k) = dum*dum1/deltat
           end if
        end if

     end do

     ! Big "administration" loop enforces conservation, updates variables
     ! that accumulate over substeps, and sets output variables.

     do i=1,mgncol

        ! get tendencies due to microphysical conversion processes
        !==========================================================
        ! note: tendencies are multiplied by appropriate cloud/precip
        ! fraction to get grid-scale values
        ! note: vap_dep is already grid-average values

        ! The net tendencies need to be added to rather than overwritten,
        ! because they may have a value already set for instantaneous
        ! melting/freezing.

        qvlat(i,k) = qvlat(i,k)-(pre(i,k)+prds(i,k))*precip_frac(i,k)-&
             vap_dep(i,k)-ice_sublim(i,k)-mnuccd(i,k)-mnudep(i,k)*lcldm(i,k) &
             -prdg(i,k)*precip_frac(i,k) 

        tlat(i,k) = tlat(i,k)+((pre(i,k)*precip_frac(i,k))*xxlv+ &
             ((prds(i,k)+prdg(i,k))*precip_frac(i,k)+vap_dep(i,k)+ice_sublim(i,k)+ &
                 mnuccd(i,k)+mnudep(i,k)*lcldm(i,k))*xxls+ &
             ((bergs(i,k)+psacws(i,k)+mnuccc(i,k)+mnucct(i,k)+msacwi(i,k)+psacwg(i,k)+ &
                  qmultg(i,k)+pgsacw(i,k))*lcldm(i,k)+ &
             (mnuccr(i,k)+pracs(i,k)+mnuccri(i,k)+pracg(i,k)+pgracs(i,k)+qmultrg(i,k))*precip_frac(i,k)+ &
                  berg(i,k))*xlf)  

        qctend(i,k) = qctend(i,k)+ &
             (-pra(i,k)-prc(i,k)-mnuccc(i,k)-mnucct(i,k)-msacwi(i,k)- &
             psacws(i,k)-bergs(i,k)-qmultg(i,k)-psacwg(i,k)-pgsacw(i,k))*lcldm(i,k)-berg(i,k)


        if (do_cldice) then

           qitend(i,k) = qitend(i,k)+ &
              (mnuccc(i,k)+mnucct(i,k)+mnudep(i,k)+msacwi(i,k)+qmultg(i,k))*lcldm(i,k)+(-prci(i,k)- &
              prai(i,k))*icldm(i,k)+vap_dep(i,k)+berg(i,k)+ice_sublim(i,k)+ &
              mnuccd(i,k)+(mnuccri(i,k)+qmultrg(i,k))*precip_frac(i,k)
 
        end if

        qrtend(i,k) = qrtend(i,k)+ &
             (pra(i,k)+prc(i,k))*lcldm(i,k)+(pre(i,k)-pracs(i,k)- &
             mnuccr(i,k)-mnuccri(i,k)-qmultrg(i,k)-pracg(i,k)-pgracs(i,k))*precip_frac(i,k)

        if (do_hail.or.do_graupel) then
           qgtend(i,k) = qgtend(i,k) + (pracg(i,k)+pgracs(i,k)+prdg(i,k)+psacr(i,k)+mnuccr(i,k))*precip_frac(i,k) &
                + (psacwg(i,k)+pgsacw(i,k))*lcldm(i,k)

           qstend(i,k) = qstend(i,k)+ &
                (prai(i,k)+prci(i,k))*icldm(i,k)+(psacws(i,k)+bergs(i,k))*lcldm(i,k)+(prds(i,k)+ &
                pracs(i,k)-psacr(i,k))*precip_frac(i,k)

        else
           !necessary since mnuccr moved to graupel
           qstend(i,k) = qstend(i,k)+ &
                (prai(i,k)+prci(i,k))*icldm(i,k)+(psacws(i,k)+bergs(i,k))*lcldm(i,k)+(prds(i,k)+ &
                pracs(i,k)+mnuccr(i,k))*precip_frac(i,k)

        end if

        cmeout(i,k) = vap_dep(i,k) + ice_sublim(i,k) + mnuccd(i,k)

        ! add output for cmei (accumulate)
        cmeitot(i,k) = vap_dep(i,k) + ice_sublim(i,k) + mnuccd(i,k)

        !-------------------------------------------------------------------
        ! evaporation/sublimation is stored here as positive term
        ! Add to evapsnow via prdg
        evapsnow(i,k) = (-prds(i,k)-prdg(i,k))*precip_frac(i,k)
        nevapr(i,k) = -pre(i,k)*precip_frac(i,k)
        prer_evap(i,k) = -pre(i,k)*precip_frac(i,k)

        ! change to make sure prain is positive: do not remove snow from
        ! prain used for wet deposition

        prain(i,k) = (pra(i,k)+prc(i,k))*lcldm(i,k)+(-pracs(i,k)- &
             mnuccr(i,k)-mnuccri(i,k))*precip_frac(i,k)
        if (do_hail .or. do_graupel) then
           prodsnow(i,k) = (prai(i,k)+prci(i,k))*icldm(i,k)+(psacws(i,k)+bergs(i,k))*lcldm(i,k)+(&
                pracs(i,k))*precip_frac(i,k)
        else
           prodsnow(i,k) = (prai(i,k)+prci(i,k))*icldm(i,k)+(psacws(i,k)+bergs(i,k))*lcldm(i,k)+(&
                pracs(i,k)+mnuccr(i,k))*precip_frac(i,k)
        end if

        ! following are used to calculate 1st order conversion rate of cloud water
        !    to rain and snow (1/s), for later use in aerosol wet removal routine
        ! previously, wetdepa used (prain/qc) for this, and the qc in wetdepa may be smaller than the qc
        !    used to calculate pra, prc, ... in this routine
        ! qcsinksum_rate1ord = { rate of direct transfer of cloud water to rain & snow }
        !                      (no cloud ice or bergeron terms)

        qcsinksum_rate1ord(i,k) = (pra(i,k)+prc(i,k)+psacws(i,k)+psacwg(i,k)+pgsacw(i,k))*lcldm(i,k) 
        ! Avoid zero/near-zero division.
        qcsinksum_rate1ord(i,k) = qcsinksum_rate1ord(i,k) / &
             max(qc(i,k),1.0e-30_r8)


        ! microphysics output, note this is grid-averaged
        pratot(i,k) = pra(i,k)*lcldm(i,k)
        prctot(i,k) = prc(i,k)*lcldm(i,k)
        mnuccctot(i,k) = mnuccc(i,k)*lcldm(i,k)
        mnuccttot(i,k) = mnucct(i,k)*lcldm(i,k)
        msacwitot(i,k) = msacwi(i,k)*lcldm(i,k)
        psacwstot(i,k) = psacws(i,k)*lcldm(i,k)
        bergstot(i,k) = bergs(i,k)*lcldm(i,k)
        bergtot(i,k) = berg(i,k)
        prcitot(i,k) = prci(i,k)*icldm(i,k)
        praitot(i,k) = prai(i,k)*icldm(i,k)
        mnuccdtot(i,k) = mnuccd(i,k)*icldm(i,k)

        pracstot(i,k) = pracs(i,k)*precip_frac(i,k)
        mnuccrtot(i,k) = mnuccr(i,k)*precip_frac(i,k)
        mnuccritot(i,k) = mnuccri(i,k)*precip_frac(i,k)

        psacrtot(i,k) = psacr(i,k)*precip_frac(i,k)
        pracgtot(i,k) = pracg(i,k)*precip_frac(i,k)
        psacwgtot(i,k) = psacwg(i,k)*lcldm(i,k)
        pgsacwtot(i,k) = pgsacw(i,k)*lcldm(i,k)
        pgracstot(i,k) = pgracs(i,k)*precip_frac(i,k)
        prdgtot(i,k) = prdg(i,k)*precip_frac(i,k)
        qmultgtot(i,k) = qmultg(i,k)*lcldm(i,k)
        qmultrgtot(i,k) = qmultrg(i,k)*precip_frac(i,k)
        npracgtot(i,k) = npracg(i,k)*precip_frac(i,k) 
        nscngtot(i,k) = nscng(i,k)*lcldm(i,k) 
        ngracstot(i,k) = ngracs(i,k)*precip_frac(i,k)
        nmultgtot(i,k) = nmultg(i,k)*lcldm(i,k)
        nmultrgtot(i,k) = nmultrg(i,k)*precip_frac(i,k)
        npsacwgtot(i,k) = npsacwg(i,k)*lcldm(i,k)

        nctend(i,k) = nctend(i,k)+&
             (-nnuccc(i,k)-nnucct(i,k)-npsacws(i,k)+nsubc(i,k) &
             -npra(i,k)-nprc1(i,k)-npsacwg(i,k))*lcldm(i,k)

        if (do_cldice) then
           if (use_hetfrz_classnuc) then
              tmpfrz = nnuccc(i,k)
           else
              tmpfrz = 0._r8
           end if

           nitend(i,k) = nitend(i,k)+ nnuccd(i,k)+ &
                (nnucct(i,k)+tmpfrz+nnudep(i,k)+nsacwi(i,k)+nmultg(i,k))*lcldm(i,k)+(nsubi(i,k)-nprci(i,k)- &
                nprai(i,k))*icldm(i,k)+(nnuccri(i,k)+nmultrg(i,k))*precip_frac(i,k)

        end if

        if(do_graupel.or.do_hail) then

           nstend(i,k) = nstend(i,k)+(nsubs(i,k)+ &
                nsagg(i,k)-ngracs(i,k))*precip_frac(i,k)+nprci(i,k)*icldm(i,k)-nscng(i,k)*lcldm(i,k)
           ngtend(i,k) = ngtend(i,k)+nscng(i,k)*lcldm(i,k)+(ngracs(i,k)+nnuccr(i,k))*precip_frac(i,k)

        else
           !necessary since mnuccr moved to graupel
           nstend(i,k) = nstend(i,k)+(nsubs(i,k)+ &
                nsagg(i,k)+nnuccr(i,k))*precip_frac(i,k)+nprci(i,k)*icldm(i,k)
        end if

        nrtend(i,k) = nrtend(i,k)+ &
             nprc(i,k)*lcldm(i,k)+(nsubr(i,k)-npracs(i,k)-nnuccr(i,k) &
             -nnuccri(i,k)+nragg(i,k)-npracg(i,k)-ngracs(i,k))*precip_frac(i,k)

        ! make sure that ni at advanced time step does not exceed
        ! maximum (existing N + source terms*dt), which is possible if mtime < deltat
        ! note that currently mtime = deltat
        !================================================================

        if (do_cldice .and. nitend(i,k).gt.0._r8.and.ni(i,k)+nitend(i,k)*deltat.gt.nimax(i,k)) then
           nitend(i,k)=max(0._r8,(nimax(i,k)-ni(i,k))/deltat)
        end if

     end do

     ! End of "administration" loop

  end do micro_vert_loop ! end k loop