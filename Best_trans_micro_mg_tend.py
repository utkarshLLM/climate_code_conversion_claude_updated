import numpy as np
from scipy.special import gamma

def micro_mg_tend(microp_uniform, pcols, pver, ncol, top_lev, deltatin,
                  tn, qn, qc, qi, nc, ni, p, pdel, cldn, liqcldf,
                  relvar, accre_enhan, icecldf, rate1ord_cw2pr_st, naai, npccnin,
                  rndst, nacon, tlat, qvlat, qctend, qitend, nctend, nitend, effc, effc_fn,
                  effi, prect, preci, nevapr, evapsnow, am_evp_st, prain, prodsnow, cmeout, deffi, pgamrad,
                  lamcrad, qsout, dsout, rflx, sflx, qrout, reff_rain, reff_snow, qcsevap, qisevap,
                  qvres, cmeiout, vtrmc, vtrmi, qcsedten, qisedten, prao, prco, mnuccco, mnuccto,
                  msacwio, psacwso, bergso, bergo, melto, homoo, qcreso, prcio, praio, qireso,
                  mnuccro, pracso, meltsdt, frzrdt, mnuccdo, nrout, nsout, refl, arefl, areflz,
                  frefl, csrfl, acsrfl, fcsrfl, rercld, ncai, ncal, qrout2, qsout2, nrout2,
                  nsout2, drout2, dsout2, freqs, freqr, nfice, prer_evap, do_cldice, errstring,
                  tnd_qsnow, tnd_nsnow, re_ice, frzimm, frzcnt, frzdep):
    
    # Constants and parameters
    r8 = np.float64
    g = 9.80616
    r = 287.04
    rv = 461.5
    cpp = 1004.64
    rhow = 1000.0
    tmelt = 273.15
    xxlv = 2.501e6
    xlf = 3.336e5
    xxls = xxlv + xlf
    rhosn = 250.0
    rhoi = 500.0
    qsmall = 1e-18
    bimm = 100.0
    aimm = 0.66
    pi = np.pi

    # Initialize output variables
    errstring = ''

    # Get physics options
    do_clubb_sgs = phys_getopts(do_clubb_sgs_out=True)

    # Initialize output fields for number concentration and ice nucleation
    ncai = np.zeros((ncol, pver), dtype=r8)
    ncal = np.zeros((ncol, pver), dtype=r8)

    # Initialize rain size
    rercld = np.zeros((ncol, pver), dtype=r8)
    arcld = np.zeros((ncol, pver), dtype=r8)

    # Initialize radiation output variables
    pgamrad = np.zeros((ncol, pver), dtype=r8)
    lamcrad = np.zeros((ncol, pver), dtype=r8)
    deffi = np.zeros((ncol, pver), dtype=r8)

    # Initialize water vapor tendency term output
    qcsevap = np.zeros((ncol, pver), dtype=r8)
    qisevap = np.zeros((ncol, pver), dtype=r8)
    qvres = np.zeros((ncol, pver), dtype=r8)
    cmeiout = np.zeros((ncol, pver), dtype=r8)
    vtrmc = np.zeros((ncol, pver), dtype=r8)
    vtrmi = np.zeros((ncol, pver), dtype=r8)
    qcsedten = np.zeros((ncol, pver), dtype=r8)
    qisedten = np.zeros((ncol, pver), dtype=r8)

    prao = np.zeros((ncol, pver), dtype=r8)
    prco = np.zeros((ncol, pver), dtype=r8)
    mnuccco = np.zeros((ncol, pver), dtype=r8)
    mnuccto = np.zeros((ncol, pver), dtype=r8)
    msacwio = np.zeros((ncol, pver), dtype=r8)
    psacwso = np.zeros((ncol, pver), dtype=r8)
    bergso = np.zeros((ncol, pver), dtype=r8)
    bergo = np.zeros((ncol, pver), dtype=r8)
    melto = np.zeros((ncol, pver), dtype=r8)
    homoo = np.zeros((ncol, pver), dtype=r8)
    qcreso = np.zeros((ncol, pver), dtype=r8)
    prcio = np.zeros((ncol, pver), dtype=r8)
    praio = np.zeros((ncol, pver), dtype=r8)
    qireso = np.zeros((ncol, pver), dtype=r8)
    mnuccro = np.zeros((ncol, pver), dtype=r8)
    pracso = np.zeros((ncol, pver), dtype=r8)
    meltsdt = np.zeros((ncol, pver), dtype=r8)
    frzrdt = np.zeros((ncol, pver), dtype=r8)
    mnuccdo = np.zeros((ncol, pver), dtype=r8)

    rflx = np.zeros((ncol, pver+1), dtype=r8)
    sflx = np.zeros((ncol, pver+1), dtype=r8)
    effc = np.zeros((ncol, pver), dtype=r8)
    effc_fn = np.zeros((ncol, pver), dtype=r8)
    effi = np.zeros((ncol, pver), dtype=r8)

    # Assign variable deltat for sub-stepping
    deltat = deltatin

    # Parameters for scheme
    omsm = 0.99999
    dto2 = 0.5 * deltat
    mincld = 0.0001

    # Initialize multi-level fields
    q = qn.copy()
    t = tn.copy()

    # Initialize time-varying parameters
    rho = p / (r * t)
    dv = 8.794e-5 * t**1.81 / p
    mu = 1.496e-6 * t**1.5 / (t + 120.0)
    sc = mu / (rho * dv)
    kap = 1.414e3 * 1.496e-6 * t**1.5 / (t + 120.0)

    # Air density adjustment for fallspeed parameters
    rhof = (85000.0 / (r * tmelt) / rho)**0.54

    # Get dz from dp and hydrostatic approximation
    dz = pdel / (rho * g)

    # Initialization
    qc[:, :top_lev-1] = 0.0
    qi[:, :top_lev-1] = 0.0
    nc[:, :top_lev-1] = 0.0
    ni[:, :top_lev-1] = 0.0
    t1 = t.copy()
    q1 = q.copy()
    qc1 = qc.copy()
    qi1 = qi.copy()
    nc1 = nc.copy()
    ni1 = ni.copy()

    # Initialize tendencies to zero
    tlat1 = np.zeros((ncol, pver), dtype=r8)
    qvlat1 = np.zeros((ncol, pver), dtype=r8)
    qctend1 = np.zeros((ncol, pver), dtype=r8)
    qitend1 = np.zeros((ncol, pver), dtype=r8)
    nctend1 = np.zeros((ncol, pver), dtype=r8)
    nitend1 = np.zeros((ncol, pver), dtype=r8)

    # Initialize precip output
    qrout = np.zeros((ncol, pver), dtype=r8)
    qsout = np.zeros((ncol, pver), dtype=r8)
    nrout = np.zeros((ncol, pver), dtype=r8)
    nsout = np.zeros((ncol, pver), dtype=r8)
    dsout = np.zeros((ncol, pver), dtype=r8)
    drout = np.zeros((ncol, pver), dtype=r8)
    reff_rain = np.zeros((ncol, pver), dtype=r8)
    reff_snow = np.zeros((ncol, pver), dtype=r8)

    # Initialize variables for trop_mozart
    nevapr = np.zeros((ncol, pver), dtype=r8)
    nevapr2 = np.zeros((ncol, pver), dtype=r8)
    evapsnow = np.zeros((ncol, pver), dtype=r8)
    prain = np.zeros((ncol, pver), dtype=r8)
    prodsnow = np.zeros((ncol, pver), dtype=r8)
    cmeout = np.zeros((ncol, pver), dtype=r8)
    am_evp_st = np.zeros((ncol, pver), dtype=r8)

    # For refl calc
    rainrt1 = np.zeros((ncol, pver), dtype=r8)

    # Initialize precip fraction and output tendencies
    cldmax = np.full((ncol, pver), mincld, dtype=r8)

    # Initialize aerosol number
    dum2l = np.zeros((ncol, pver), dtype=r8)
    dum2i = np.zeros((ncol, pver), dtype=r8)

    # Initialize avg precip rate
    prect1 = np.zeros(ncol, dtype=r8)
    preci1 = np.zeros(ncol, dtype=r8)

    # Get humidity and saturation vapor pressures
    for k in range(top_lev, pver):
        for i in range(ncol):
            es = svp_water(t[i, k])
            qs = svp_to_qsat(es, p[i, k])

            # Prevents negative values
            if qs < 0.0:
                qs = 1.0
                es = p[i, k]

            esl[i, k] = svp_water(t[i, k])
            esi[i, k] = svp_ice(t[i, k])

            # Make sure when above freezing that esi=esl
            if t[i, k] > tmelt:
                esi[i, k] = esl[i, k]

            relhum[i, k] = q[i, k] / qs

            # Get cloud fraction, check for minimum
            cldm[i, k] = max(cldn[i, k], mincld)
            cldmw[i, k] = max(cldn[i, k], mincld)
            icldm[i, k] = max(icecldf[i, k], mincld)
            lcldm[i, k] = max(liqcldf[i, k], mincld)

            # Subcolumns, set cloud fraction variables to one
            if microp_uniform:
                cldm[i, k] = mincld
                cldmw[i, k] = mincld
                icldm[i, k] = mincld
                lcldm[i, k] = mincld

                if qc[i, k] >= qsmall:
                    lcldm[i, k] = 1.0
                    cldm[i, k] = 1.0
                    cldmw[i, k] = 1.0

                if qi[i, k] >= qsmall:
                    cldm[i, k] = 1.0
                    icldm[i, k] = 1.0

            # Calculate nfice based on liquid and ice mmr
            nfice[i, k] = 0.0
            dumfice = qc[i, k] + qi[i, k]
            if dumfice > qsmall and qi[i, k] > qsmall:
                nfice[i, k] = qi[i, k] / dumfice

            if do_cldice and (t[i, k] < tmelt - 5.0):
                # If aerosols interact with ice set number of activated ice nuclei
                dum2 = naai[i, k]

                dumnnuc = (dum2 - ni[i, k] / icldm[i, k]) / deltat * icldm[i, k]
                dumnnuc = max(dumnnuc, 0.0)
                # Get provisional ni and qi after nucleation in order to calculate
                # Bergeron process below
                ninew = ni[i, k] + dumnnuc * deltat
                qinew = qi[i, k] + dumnnuc * deltat * mi0
            else:
                ninew = ni[i, k]
                qinew = qi[i, k]

            # Initialize CME components
            cme[i, k] = 0.0
            cmei[i, k] = 0.0

            # Bergeron process
            berg[i, k] = 0.0
            prd = 0.0

            # Get in-cloud qi and ni after nucleation
            if icldm[i, k] > 0.0:
                qiic[i, k] = qinew / icldm[i, k]
                niic[i, k] = ninew / icldm[i, k]
            else:
                qiic[i, k] = 0.0
                niic[i, k] = 0.0

            if nicons:
                niic[i, k] = ninst / rho[i, k]

            # If T < 0 C then bergeron
            if do_cldice and (t[i, k] < 273.15):
                # If ice exists
                if qi[i, k] > qsmall:
                    bergtsf = 0.0  # Bergeron time scale (fraction of timestep)

                    qvi = svp_to_qsat(esi[i, k], p[i, k])
                    qvl = svp_to_qsat(esl[i, k], p[i, k])

                    dqsidt = xxls * qvi / (rv * t[i, k]**2)
                    abi = 1.0 + dqsidt * xxls / cpp

                    # Get ice size distribution parameters
                    if qiic[i, k] >= qsmall:
                        lami[k] = (cons1 * ci * niic[i, k] / qiic[i, k])**(1.0 / di)
                        n0i[k] = niic[i, k] * lami[k]

                        # Check for slope
                        # Adjust vars
                        if lami[k] < lammini:
                            lami[k] = lammini
                            n0i[k] = lami[k]**(di + 1.0) * qiic[i, k] / (ci * cons1)
                        elif lami[k] > lammaxi:
                            lami[k] = lammaxi
                            n0i[k] = lami[k]**(di + 1.0) * qiic[i, k] / (ci * cons1)

                        epsi = 2.0 * pi * n0i[k] * rho[i, k] * Dv[i, k] / (lami[k] * lami[k])

                        # If liquid exists
                        if qc[i, k] > qsmall:
                            # Begin bergeron process
                            # Calculate Bergeron process
                            prd = epsi * (qvl - qvi) / abi
                        else:
                            prd = 0.0

                        # Multiply by cloud fraction
                        prd = prd * min(icldm[i, k], lcldm[i, k])

                        # Transfer of existing cloud liquid to ice
                        berg[i, k] = max(0.0, prd)

                    if berg[i, k] > 0.0:
                        bergtsf = max(0.0, (qc[i, k] / berg[i, k]) / deltat)

                        if bergtsf < 1.0:
                            berg[i, k] = max(0.0, qc[i, k] / deltat)

                    if bergtsf < 1.0 or icldm[i, k] > lcldm[i, k]:
                        if qiic[i, k] >= qsmall:
                            # First case is for case when liquid water is present, but is completely depleted
                            # in time step, i.e., bergrsf > 0 but < 1
                            if qc[i, k] >= qsmall:
                                rhin = (1.0 + relhum[i, k]) / 2.0
                                if (rhin * esl[i, k] / esi[i, k]) > 1.0:
                                    prd = epsi * (rhin * qvl - qvi) / abi

                                    # Multiply by cloud fraction assuming liquid/ice maximum overlap
                                    prd = prd * min(icldm[i, k], lcldm[i, k])

                                    # Add to cmei
                                    cmei[i, k] = cmei[i, k] + (prd * (1.0 - bergtsf))

                            # Second case is for pure ice cloud, either no liquid, or icldm > lcldm
                            if qc[i, k] < qsmall or icldm[i, k] > lcldm[i, k]:
                                # Note: for case of no liquid, need to set liquid cloud fraction to zero
                                # Store liquid cloud fraction in 'dum'
                                if qc[i, k] < qsmall:
                                    dum = 0.0
                                else:
                                    dum = lcldm[i, k]

                                # Set RH to grid-mean value for pure ice cloud
                                rhin = relhum[i, k]

                                if (rhin * esl[i, k] / esi[i, k]) > 1.0:
                                    prd = epsi * (rhin * qvl - qvi) / abi

                                    # Multiply by relevant cloud fraction for pure ice cloud
                                    # Assuming maximum overlap of liquid/ice
                                    prd = prd * max((icldm[i, k] - dum), 0.0)
                                    cmei[i, k] = cmei[i, k] + prd

                    # If deposition, it should not reduce grid mean rhi below 1.0
                    if cmei[i, k] > 0.0 and (relhum[i, k] * esl[i, k] / esi[i, k]) > 1.0:
                        cmei[i, k] = min(cmei[i, k], (q[i, k] - qs * esi[i, k] / esl[i, k]) / abi / deltat)

            # Evaporation should not exceed available water
            if (-berg[i, k]) < -qc[i, k] / deltat:
                berg[i, k] = max(qc[i, k] / deltat, 0.0)

            # Sublimation process
            if do_cldice and ((relhum[i, k] * esl[i, k] / esi[i, k]) < 1.0 and qiic[i, k] >= qsmall):
                qvi = svp_to_qsat(esi[i, k], p[i, k])
                qvl = svp_to_qsat(esl[i, k], p[i, k])
                dqsidt = xxls * qvi / (rv * t[i, k]**2)
                abi = 1.0 + dqsidt * xxls / cpp

                # Get ice size distribution parameters
                lami[k] = (cons1 * ci * niic[i, k] / qiic[i, k])**(1.0 / di)
                n0i[k] = niic[i, k] * lami[k]

                # Check for slope
                # Adjust vars
                if lami[k] < lammini:
                    lami[k] = lammini
                    n0i[k] = lami[k]**(di + 1.0) * qiic[i, k] / (ci * cons1)
                elif lami[k] > lammaxi:
                    lami[k] = lammaxi
                    n0i[k] = lami[k]**(di + 1.0) * qiic[i, k] / (ci * cons1)

                epsi = 2.0 * pi * n0i[k] * rho[i, k] * Dv[i, k] / (lami[k] * lami[k])

                # Modify for ice fraction below
                prd = epsi * (relhum[i, k] * qvl - qvi) / abi * icldm[i, k]
                cmei[i, k] = min(prd, 0.0)

            # Sublimation should not exceed available ice
            if cmei[i, k] < -qi[i, k] / deltat:
                cmei[i, k] = -qi[i, k] / deltat

            # Sublimation should not increase grid mean rhi above 1.0
            if cmei[i, k] < 0.0 and (relhum[i, k] * esl[i, k] / esi[i, k]) < 1.0:
                cmei[i, k] = min(0.0, max(cmei[i, k], (q[i, k] - qs * esi[i, k] / esl[i, k]) / abi / deltat))

            # Limit cmei due for roundoff error
            cmei[i, k] = cmei[i, k] * omsm

            # Conditional for ice nucleation
            if do_cldice and (t[i, k] < (tmelt - 5.0)):
                # Using Liu et al. (2007) ice nucleation with hooks into simulated aerosol
                # Ice nucleation rate (dum2) has already been calculated and read in (naai)
                dum2i[i, k] = naai[i, k]
            else:
                dum2i[i, k] = 0.0

    # Initialize sub-step precip flux variables
    rflx1 = np.zeros((ncol, pver+1), dtype=r8)
    sflx1 = np.zeros((ncol, pver+1), dtype=r8)

    # Initialize final precip flux variables
    rflx = np.zeros((ncol, pver+1), dtype=r8)
    sflx = np.zeros((ncol, pver+1), dtype=r8)

    ltrue = np.zeros(ncol, dtype=bool)
    for i in range(ncol):
        for k in range(top_lev, pver):
            # Skip microphysical calculations if no cloud water
            if qc[i, k] >= qsmall or qi[i, k] >= qsmall or cmei[i, k] >= qsmall:
                ltrue[i] = True

    # Assign number of sub-steps to iter
    # Use 2 sub-steps, following tests described in MG2008
    iter = 2

    # Get sub-step time step
    deltat = deltat / iter

    # Since activation/nucleation processes are fast, need to take into account
    # factor mtime = mixing timescale in cloud / model time step
    # Mixing time can be interpreted as cloud depth divided by sub-grid vertical velocity
    # For now mixing timescale is assumed to be 1 timestep for modal aerosols, 20 min bulk
    mtime = 1.0
    rate1ord_cw2pr_st[:, :] = 0.0

    # Skip calculations if no cloud water
    for i in range(ncol):
        if not ltrue[i]:
            tlat[i, :] = 0.0
            qvlat[i, :] = 0.0
            qctend[i, :] = 0.0
            qitend[i, :] = 0.0
            qnitend[i, :] = 0.0
            qrtend[i, :] = 0.0
            nctend[i, :] = 0.0
            nitend[i, :] = 0.0
            nrtend[i, :] = 0.0
            nstend[i, :] = 0.0
            prect[i] = 0.0
            preci[i] = 0.0
            rflx[i, :] = 0.0
            sflx[i, :] = 0.0
            qniic[i, :] = 0.0
            qric[i, :] = 0.0
            nsic[i, :] = 0.0
            nric[i, :] = 0.0
            rainrt[i, :] = 0.0
            continue

        qcsinksum_rate1ord = np.zeros(pver, dtype=r8)
        qcsum_rate1ord = np.zeros(pver, dtype=r8)

        # Begin sub-step
        for it in range(iter):
            # Initialize sub-step microphysical tendencies
            tlat[i, :] = 0.0
            qvlat[i, :] = 0.0
            qctend[i, :] = 0.0
            qitend[i, :] = 0.0
            qnitend[i, :] = 0.0
            qrtend[i, :] = 0.0
            nctend[i, :] = 0.0
            nitend[i, :] = 0.0
            nrtend[i, :] = 0.0
            nstend[i, :] = 0.0

            # Initialize diagnostic precipitation to zero
            qniic[i, :] = 0.0
            qric[i, :] = 0.0
            nsic[i, :] = 0.0
            nric[i, :] = 0.0
            rainrt[i, :] = 0.0

            # Initialize vertically-integrated rain and snow tendencies
            qrtot = 0.0
            nrtot = 0.0
            qstot = 0.0
            nstot = 0.0

            # Initialize precip at surface
            prect[i] = 0.0
            preci[i] = 0.0

            # Initialize fluxes
            rflx[i, :] = 0.0
            sflx[i, :] = 0.0

            for k in range(top_lev, pver):
                qcvar = relvar[i, k]
                cons2 = gamma(qcvar + 2.47)
                cons3 = gamma(qcvar)
                cons9 = gamma(qcvar + 2.0)
                cons10 = gamma(qcvar + 1.0)
                cons12 = gamma(qcvar + 1.15)
                cons15 = gamma(qcvar + bc / 3.0)
                cons18 = qcvar**2.47
                cons19 = qcvar**2
                cons20 = qcvar**1.15

                # Set cwml and cwmi to current qc and qi
                cwml[i, k] = qc[i, k]
                cwmi[i, k] = qi[i, k]

                # Initialize precip fallspeeds to zero
                ums[k] = 0.0
                uns[k] = 0.0
                umr[k] = 0.0
                unr[k] = 0.0

                # Calculate precip fraction based on maximum overlap assumption
                if k == top_lev:
                    cldmax[i, k] = cldm[i, k]
                else:
                    # If rain or snow mix ratio is smaller than
                    # threshold, then set cldmax to cloud fraction at current level
                    if do_clubb_sgs:
                        if qc[i, k] >= qsmall or qi[i, k] >= qsmall:
                            cldmax[i, k] = cldm[i, k]
                        else:
                            cldmax[i, k] = cldmax[i, k-1]
                    else:
                        if qric[i, k-1] >= qsmall or qniic[i, k-1] >= qsmall:
                            cldmax[i, k] = max(cldmax[i, k-1], cldm[i, k])
                        else:
                            cldmax[i, k] = cldm[i, k]

                # Decrease in number concentration due to sublimation/evap
                # Divide by cloud fraction to get in-cloud decrease
                # Don't reduce Nc due to bergeron process
                if cmei[i, k] < 0.0 and qi[i, k] > qsmall and cldm[i, k] > mincld:
                    nsubi[k] = cmei[i, k] / qi[i, k] * ni[i, k] / cldm[i, k]
                else:
                    nsubi[k] = 0.0
                nsubc[k] = 0.0

                # Ice nucleation if activated nuclei exist at t<-5C AND rhmini + 5%
                if do_cldice and dum2i[i, k] > 0.0 and t[i, k] < (tmelt - 5.0) and \
                   relhum[i, k] * esl[i, k] / esi[i, k] > rhmini + 0.05:
                    # If NCAI > 0. then set numice = ncai (as before)
                    # Note: this is gridbox averaged
                    nnuccd[k] = (dum2i[i, k] - ni[i, k] / icldm[i, k]) / deltat * icldm[i, k]
                    nnuccd[k] = max(nnuccd[k], 0.0)
                    nimax = dum2i[i, k] * icldm[i, k]

                    # Calc mass of new particles using new crystal mass...
                    # Also this will be multiplied by mtime as nnuccd is...
                    mnuccd[k] = nnuccd[k] * mi0

                    # Add mnuccd to cmei....
                    cmei[i, k] = cmei[i, k] + mnuccd[k] * mtime

                    # Limit cmei
                    qvi = svp_to_qsat(esi[i, k], p[i, k])
                    dqsidt = xxls * qvi / (rv * t[i, k]**2)
                    abi = 1.0 + dqsidt * xxls / cpp
                    cmei[i, k] = min(cmei[i, k], (q[i, k] - qvi) / abi / deltat)

                    # Limit for roundoff error
                    cmei[i, k] = cm
                    # Limit for roundoff error
                    cmei[i, k] = cmei[i, k] * omsm

                    # Add to output
                    qitend[i, k] = qitend[i, k] + mnuccd[k] * mtime
                    nitend[i, k] = nitend[i, k] + nnuccd[k] * mtime
                else:
                    nnuccd[k] = 0.0
                    mnuccd[k] = 0.0

                # Immersion freezing (Bigg, 1953)
                if do_cldice and qc[i, k] > qsmall and t[i, k] < tmelt and t[i, k] > tmin:
                    nimm = bimm * (np.exp(aimm * (tmelt - t[i, k])) - 1.0) * qc[i, k] * rhow / 1.0e6
                    nimm = nimm * deltat * lcldm[i, k]
                    qimm = nimm * rhow * pi / 6.0 * 10.0e-18
                else:
                    nimm = 0.0
                    qimm = 0.0

                # Contact freezing (Young, 1974) with hooks into simulated dust
                if do_cldice and qc[i, k] > qsmall and t[i, k] < tmelt and t[i, k] > tmin:
                    ncontact = 2.0 * pi * 0.4e-6 * qc[i, k] * rhow / droplet_mass1 * \
                               (1.0e-4 / 0.4e-6) * nacon[i, k] * deltat * lcldm[i, k]
                    qcontact = ncontact * rhow * pi / 6.0 * 10.0e-18
                else:
                    ncontact = 0.0
                    qcontact = 0.0

                # Add contact and immersion freezing tendencies
                qctend[i, k] = qctend[i, k] - qimm - qcontact
                nctend[i, k] = nctend[i, k] - nimm - ncontact
                qitend[i, k] = qitend[i, k] + qimm + qcontact
                nitend[i, k] = nitend[i, k] + nimm + ncontact

                # Assign tendencies to freezing processes
                frzimm[i, k] = qimm
                frzcnt[i, k] = qcontact

                # Autoconversion of cloud liquid water to rain
                if qc[i, k] > qsmall:
                    # Get in-cloud values
                    qcic = qc[i, k] / lcldm[i, k]
                    ncic = nc[i, k] / lcldm[i, k]

                    # Autoconversion rate (Khairoutdinov and Kogan 2000)
                    if qcic > qsmall:
                        prc = 1350.0 * qcic**2.47 * (ncic * 1.0e-6)**(-1.79)
                    else:
                        prc = 0.0

                    # Autoconversion tends to reduce nc (not always explicitly
                    # calculated in bulk schemes)
                    nprc = prc * (nc[i, k] / qc[i, k]) * 1.5

                    # Adjust tendencies
                    qctend[i, k] = qctend[i, k] - prc
                    nctend[i, k] = nctend[i, k] - nprc
                    qrtend[i, k] = qrtend[i, k] + prc
                    nrtend[i, k] = nrtend[i, k] + nprc / 20.0

                    # Add to precip variables
                    rainrt[i, k] = rainrt[i, k] + prc

                    # Add to output
                    prao[i, k] = prc

                    # For 1st order conversion rate
                    qcsinksum_rate1ord[k] = qcsinksum_rate1ord[k] + prc
                    qcsum_rate1ord[k] = qcsum_rate1ord[k] + qc[i, k]
                else:
                    prc = 0.0
                    nprc = 0.0

                # Autoconversion of cloud ice to snow
                # Similar to Ferrier (1994)
                if qi[i, k] > qsmall:
                    # Get in-cloud values
                    qiic = qi[i, k] / icldm[i, k]
                    niic = ni[i, k] / icldm[i, k]

                    # Autoconversion rate
                    if qiic > qsmall:
                        lami = (cons1 * niic / qiic)**(1.0 / di)
                        n0i = niic * lami
                        prci = cons6 * n0i * rhow / (lami**(di + 1.0))
                    else:
                        prci = 0.0

                    # Adjust tendencies
                    qitend[i, k] = qitend[i, k] - prci
                    nitend[i, k] = nitend[i, k] - prci * (niic / qiic)
                    qnitend[i, k] = qnitend[i, k] + prci
                    nstend[i, k] = nstend[i, k] + prci * (niic / qiic) / 400.0

                    # Add to precip variables
                    snowrt[i, k] = snowrt[i, k] + prci

                    # Add to output
                    prcio[i, k] = prci
                else:
                    prci = 0.0

                # Accretion of cloud droplets by rain
                if qr[i, k] > qsmall and qc[i, k] > qsmall:
                    # Get in-cloud values
                    qcic = qc[i, k] / lcldm[i, k]
                    ncic = nc[i, k] / lcldm[i, k]
                    qric = qr[i, k] / cldmax[i, k]
                    nric = nr[i, k] / cldmax[i, k]

                    # Accretion rate
                    lamr = (cons1 * nric / qric)**(1.0 / dr)
                    n0r = nric * lamr
                    pra = cons7 * n0r * qcic * gamma(br + 3.0) / \
                          (lamr**(br + 3.0))

                    # Adjust tendencies
                    qctend[i, k] = qctend[i, k] - pra
                    nctend[i, k] = nctend[i, k] - pra * (ncic / qcic)
                    qrtend[i, k] = qrtend[i, k] + pra

                    # Add to precip variables
                    rainrt[i, k] = rainrt[i, k] + pra

                    # Add to output
                    prco[i, k] = pra

                    # For 1st order conversion rate
                    qcsinksum_rate1ord[k] = qcsinksum_rate1ord[k] + pra
                    qcsum_rate1ord[k] = qcsum_rate1ord[k] + qc[i, k]
                else:
                    pra = 0.0

                # Self-collection of rain drops
                if qr[i, k] > qsmall:
                    # Get in-cloud values
                    qric = qr[i, k] / cldmax[i, k]
                    nric = nr[i, k] / cldmax[i, k]

                    # Self-collection rate
                    lamr = (cons1 * nric / qric)**(1.0 / dr)
                    n0r = nric * lamr
                    nragg = cons8 * n0r * lamr**(-br - 1.0)

                    # Adjust tendencies
                    nrtend[i, k] = nrtend[i, k] - nragg
                else:
                    nragg = 0.0

                # Accretion of cloud ice by snow
                if qs[i, k] > qsmall and qi[i, k] > qsmall:
                    # Get in-cloud values
                    qiic = qi[i, k] / icldm[i, k]
                    niic = ni[i, k] / icldm[i, k]
                    qsic = qs[i, k] / cldmax[i, k]
                    nsic = ns[i, k] / cldmax[i, k]

                    # Accretion rate
                    lams = (cons1 * nsic / qsic)**(1.0 / ds)
                    n0s = nsic * lams
                    prai = cons4 * n0s * qiic * gamma(bs + 3.0) / \
                           (lams**(bs + 3.0))

                    # Adjust tendencies
                    qitend[i, k] = qitend[i, k] - prai
                    nitend[i, k] = nitend[i, k] - prai * (niic / qiic)
                    qnitend[i, k] = qnitend[i, k] + prai

                    # Add to precip variables
                    snowrt[i, k] = snowrt[i, k] + prai

                    # Add to output
                    praio[i, k] = prai
                else:
                    prai = 0.0

                # Self-collection of snow
                if qs[i, k] > qsmall:
                    # Get in-cloud values
                    qsic = qs[i, k] / cldmax[i, k]
                    nsic = ns[i, k] / cldmax[i, k]

                    # Self-collection rate
                    lams = (cons1 * nsic / qsic)**(1.0 / ds)
                    n0s = nsic * lams
                    nsagg = cons5 * n0s * lams**(-bs - 1.0)

                    # Adjust tendencies
                    nstend[i, k] = nstend[i, k] - nsagg
                else:
                    nsagg = 0.0

                # Heterogeneous freezing of rain drops
                if qr[i, k] > qsmall and t[i, k] < tmelt:
                    # Get in-cloud values
                    qric = qr[i, k] / cldmax[i, k]
                    nric = nr[i, k] / cldmax[i, k]

                    # Freezing rate
                    lamr = (cons1 * nric / qric)**(1.0 / dr)
                    n0r = nric * lamr
                    prfz = cons9 * n0r * (np.exp(aimm * (tmelt - t[i, k])) - 1.0) / \
                           lamr**(br + 1.0)

                    # Adjust tendencies
                    qrtend[i, k] = qrtend[i, k] - prfz
                    nrtend[i, k] = nrtend[i, k] - prfz * (nric / qric)
                    qnitend[i, k] = qnitend[i, k] + prfz
                    nstend[i, k] = nstend[i, k] + prfz * (nric / qric)

                    # Add to precip variables
                    snowrt[i, k] = snowrt[i, k] + prfz

                    # Add to output
                    mnuccro[i, k] = prfz
                else:
                    prfz = 0.0

                # Melting of snow
                if qs[i, k] > qsmall and t[i, k] > tmelt:
                    # Get in-cloud values
                    qsic = qs[i, k] / cldmax[i, k]
                    nsic = ns[i, k] / cldmax[i, k]

                    # Melting rate
                    lams = (cons1 * nsic / qsic)**(1.0 / ds)
                    n0s = nsic * lams
                    psmlt = cons10 * n0s * (t[i, k] - tmelt) / \
                            (lams**(ds + 1.0) * rhof[i, k])

                    # Adjust tendencies
                    qnitend[i, k] = qnitend[i, k] - psmlt
                    nstend[i, k] = nstend[i, k] - psmlt * (nsic / qsic)
                    qrtend[i, k] = qrtend[i, k] + psmlt
                    nrtend[i, k] = nrtend[i, k] + psmlt * (nsic / qsic)

                    # Add to precip variables
                    rainrt[i, k] = rainrt[i, k] + psmlt

                    # Add to output
                    melto[i, k] = psmlt
                else:
                    psmlt = 0.0

                # Homogeneous freezing of cloud water and rain
                if t[i, k] < thom:
                    if qc[i, k] > qsmall:
                        qctend[i, k] = qctend[i, k] - qc[i, k] / deltat
                        nctend[i, k] = nctend[i, k] - nc[i, k] / deltat
                        qitend[i, k] = qitend[i, k] + qc[i, k] / deltat
                        nitend[i, k] = nitend[i, k] + nc[i, k] / deltat
                        homoc = qc[i, k] / deltat
                    else:
                        homoc = 0.0

                    if qr[i, k] > qsmall:
                        qrtend[i, k] = qrtend[i, k] - qr[i, k] / deltat
                        nrtend[i, k] = nrtend[i, k] - nr[i, k] / deltat
                        qnitend[i, k] = qnitend[i, k] + qr[i, k] / deltat
                        nstend[i, k] = nstend[i, k] + nr[i, k] / deltat
                        homor = qr[i, k] / deltat
                    else:
                        homor = 0.0

                    # Add to output
                    homoo[i, k] = homoc + homor
                else:
                    homoc = 0.0
                    homor = 0.0

                # Calculate sedimentation velocities
                if qs[i, k] > qsmall:
                    lams = (cons1 * nsic / qsic)**(1.0 / ds)
                    n0s = nsic * lams
                    ums[k] = cons11 / lams
                    uns[k] = cons12 / lams
                else:
                    ums[k] = 0.0
                    uns[k] = 0.0

                if qr[i, k] > qsmall:
                    lamr = (cons1 * nric / qric)**(1.0 / dr)
                    n0r = nric * lamr
                    umr[k] = cons13 / lamr
                    unr[k] = cons14 / lamr
                else:
                    umr[k] = 0.0
                    unr[k] = 0.0

                # Adjust sedimentation velocities in sub-columns
                ums[k] = ums[k] * rhof[i, k]
                uns[k] = uns[k] * rhof[i, k]
                umr[k] = umr[k] * rhof[i, k]
                unr[k] = unr[k] * rhof[i, k]

                # Calculate fluxes
                if k == top_lev:
                    sflx[i, k] = ums[k] * qs[i, k] + uns[k] * ns[i, k]
                    rflx[i, k] = umr[k] * qr[i, k] + unr[k] * nr[i, k]
                else:
                    sflx[i, k] = ums[k] * qs[i, k] + uns[k] * ns[i, k] + \
                                 sflx[i, k-1] * (1.0 - cldmax[i, k-1]) / \
                                 (1.0 - cldmax[i, k])
                    rflx[i, k] = umr[k] * qr[i, k] + unr[k] * nr[i, k] + \
                                 rflx[i, k-1] * (1.0 - cldmax[i, k-1]) / \
                                 (1.0 - cldmax[i, k])

                # Calculate tendencies due to sedimentation
                if k == pver:
                    qnitend[i, k] = qnitend[i, k] - sflx[i, k] / dz[i, k]
                    nstend[i, k] = nstend[i, k] - sflx[i, k] / dz[i, k] * \
                                   (ns[i, k] / qs[i, k])
                    qrtend[i, k] = qrtend[i, k] - rflx[i, k] / dz[i, k]
                    nrtend[i, k] = nrtend[i, k] - rflx[i, k] / dz[i, k] * \
                                   (nr[i, k] / qr[i, k])
                else:
                    qnitend[i, k] = qnitend[i, k] + (sflx[i, k+1] - sflx[i, k]) / dz[i, k]
                    nstend[i, k] = nstend[i, k] + (sflx[i, k+1] - sflx[i, k]) / dz[i, k] * \
                                   (ns[i, k] / qs[i, k])
                    qrtend[i, k] = qrtend[i, k] + (rflx[i, k+1] - rflx[i, k]) / dz[i, k]
                    nrtend[i, k] = nrtend[i, k] + (rflx[i, k+1] - rflx[i, k]) / dz[i, k] * \
                                   (nr[i, k] / qr[i, k])

                # Calculate precipitation rates
                if k == pver:
                    prect[i] = prect[i] + rflx[i, k]
                    preci[i] = preci[i] + sflx[i, k]

                # Update state variables
                t[i, k] = t[i, k] + tlat[i, k] * deltat
                q[i, k] = q[i, k] + qvlat[i, k] * deltat
                qc[i, k] = max(qc[i, k] + qctend[i, k] * deltat, 0.0)
                qi[i, k] = max(qi[i, k] + qitend[i, k] * deltat, 0.0)
                qr[i, k] = max(qr[i, k] + qrtend[i, k] * deltat, 0.0)
                qs[i, k] = max(qs[i, k] + qnitend[i, k] * deltat, 0.0)
                nc[i, k] = max(nc[i, k] + nctend[i, k] * deltat, 0.0)
                ni[i, k] = max(ni[i, k] + nitend[i, k] * deltat, 0.0)
                nr[i, k] = max(nr[i, k] + nrtend[i, k] * deltat, 0.0)
                ns[i, k] = max(ns[i, k] + nstend[i, k] * deltat, 0.0)

            # End of vertical loop

        # End of sub-step loop

        # Calculate 1st order conversion rate
        for k in range(top_lev, pver):
            if qcsum_rate1ord[k] > 0.0:
                rate1ord_cw2pr_st[i, k] = qcsinksum_rate1ord[k] / qcsum_rate1ord[k] / iter
            else:
                rate1ord_cw2pr_st[i, k] = 0.0

    # End of column loop

    # Return output variables
    return (tlat, qvlat, qctend, qitend, qnitend, qrtend, nctend, nitend, nrtend, nstend,
            prect, preci, sflx, rflx, qniic, qric, nsic, nric, rainrt, effc, effc_fn, effi,
            sadice, sadsnow, pragt, psacws, pracg, psacwg, pgsacw, pgracs, prdg, qmultg,
            qmultrg, psacr, npracg, npsacwg, nscng, ngracs, nmultg, nmultrg, nsubg, eva_mpg,
            pra_mpg, acrr_mpg, mnuccr_mpg, msacr_mpg, pwsub_mpg, pgsub_mpg, evapmssnow_mpg,
            pgmlt_mpg, npra_mpg, nsubr_mpg, prc_mpg, nprc_mpg, npccn_mpg, npsacws_mpg,
            nsubg_mpg, psacws_mpg, pracg_mpg, mnuccc_mpg, mnuccr_mpg, msacwi_mpg, psacwi_mpg,
            npsacwi_mpg, nsubr_mpg, pgsacw_mpg, pmltg_mpg, psacr_mpg, nmultg_mpg, nmultr_mpg,
            qmultg_mpg, qmultr_mpg, pracs_mpg, npracs_mpg, psacws_mpg, npsacws_mpg, qcreso_mpg,
            prcio_mpg, praio_mpg, qireso_mpg, prai_mpg, qnireso_mpg, prci_mpg, prai_mpg)