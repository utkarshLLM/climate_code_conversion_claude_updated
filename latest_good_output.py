import numpy as np
import math

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
    ac = 3e7
    bc = 2.0
    as_ = 11.72
    bs = 0.41
    ai = 700.0
    bi = 1.0
    ar = 841.99667
    br = 0.8
    pi = 3.14159265358979323846

    # Initialize output variables
    errstring = ''

    # Get physics options
    do_clubb_sgs = phys_getopts(do_clubb_sgs_out=True)

    # Initialize output fields for number concentration and ice nucleation
    ncai[:ncol, :pver] = 0.0
    ncal[:ncol, :pver] = 0.0

    # Initialize rain size
    rercld[:ncol, :pver] = 0.0
    arcld[:ncol, :pver] = 0.0

    # Initialize radiation output variables
    pgamrad[:ncol, :pver] = 0.0
    lamcrad[:ncol, :pver] = 0.0
    deffi[:ncol, :pver] = 0.0

    # Initialize water vapor tendency term output
    qcsevap[:ncol, :pver] = 0.0
    qisevap[:ncol, :pver] = 0.0
    qvres[:ncol, :pver] = 0.0
    cmeiout[:ncol, :pver] = 0.0
    vtrmc[:ncol, :pver] = 0.0
    vtrmi[:ncol, :pver] = 0.0
    qcsedten[:ncol, :pver] = 0.0
    qisedten[:ncol, :pver] = 0.0

    prao[:ncol, :pver] = 0.0
    prco[:ncol, :pver] = 0.0
    mnuccco[:ncol, :pver] = 0.0
    mnuccto[:ncol, :pver] = 0.0
    msacwio[:ncol, :pver] = 0.0
    psacwso[:ncol, :pver] = 0.0
    bergso[:ncol, :pver] = 0.0
    bergo[:ncol, :pver] = 0.0
    melto[:ncol, :pver] = 0.0
    homoo[:ncol, :pver] = 0.0
    qcreso[:ncol, :pver] = 0.0
    prcio[:ncol, :pver] = 0.0
    praio[:ncol, :pver] = 0.0
    qireso[:ncol, :pver] = 0.0
    mnuccro[:ncol, :pver] = 0.0
    pracso[:ncol, :pver] = 0.0
    meltsdt[:ncol, :pver] = 0.0
    frzrdt[:ncol, :pver] = 0.0
    mnuccdo[:ncol, :pver] = 0.0

    rflx[:, :] = 0.0
    sflx[:, :] = 0.0
    effc[:, :] = 0.0
    effc_fn[:, :] = 0.0
    effi[:, :] = 0.0

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
    rho = np.zeros((ncol, pver))
    dv = np.zeros((ncol, pver))
    mu = np.zeros((ncol, pver))
    sc = np.zeros((ncol, pver))
    kap = np.zeros((ncol, pver))
    rhof = np.zeros((ncol, pver))
    arn = np.zeros((ncol, pver))
    asn = np.zeros((ncol, pver))
    acn = np.zeros((ncol, pver))
    ain = np.zeros((ncol, pver))
    dz = np.zeros((ncol, pver))

    for k in range(pver):
        for i in range(ncol):
            rho[i, k] = p[i, k] / (r * t[i, k])
            dv[i, k] = 8.794e-5 * t[i, k]**1.81 / p[i, k]
            mu[i, k] = 1.496e-6 * t[i, k]**1.5 / (t[i, k] + 120.0)
            sc[i, k] = mu[i, k] / (rho[i, k] * dv[i, k])
            kap[i, k] = 1.414e3 * 1.496e-6 * t[i, k]**1.5 / (t[i, k] + 120.0)

            # Air density adjustment for fallspeed parameters
            rhof[i, k] = (85000.0 / (r * tmelt) / rho[i, k])**0.54

            arn[i, k] = ar * rhof[i, k]
            asn[i, k] = as_ * rhof[i, k]
            acn[i, k] = ac * rhof[i, k]
            ain[i, k] = ai * rhof[i, k]

            # Get dz from dp and hydrostatic approximation
            dz[i, k] = pdel[i, k] / (rho[i, k] * g)

    # Initialization
    qc[:ncol, :top_lev-1] = 0.0
    qi[:ncol, :top_lev-1] = 0.0
    nc[:ncol, :top_lev-1] = 0.0
    ni[:ncol, :top_lev-1] = 0.0
    t1 = t.copy()
    q1 = q.copy()
    qc1 = qc.copy()
    qi1 = qi.copy()
    nc1 = nc.copy()
    ni1 = ni.copy()

    # Initialize tendencies to zero
    tlat1 = np.zeros((ncol, pver))
    qvlat1 = np.zeros((ncol, pver))
    qctend1 = np.zeros((ncol, pver))
    qitend1 = np.zeros((ncol, pver))
    nctend1 = np.zeros((ncol, pver))
    nitend1 = np.zeros((ncol, pver))

    # Initialize precip output
    qrout[:ncol, :pver] = 0.0
    qsout[:ncol, :pver] = 0.0
    nrout[:ncol, :pver] = 0.0
    nsout[:ncol, :pver] = 0.0
    dsout[:ncol, :pver] = 0.0

    drout[:ncol, :pver] = 0.0

    reff_rain[:ncol, :pver] = 0.0
    reff_snow[:ncol, :pver] = 0.0

    # Initialize variables for trop_mozart
    nevapr[:ncol, :pver] = 0.0
    nevapr2 = np.zeros((ncol, pver))
    evapsnow[:ncol, :pver] = 0.0
    prain[:ncol, :pver] = 0.0
    prodsnow[:ncol, :pver] = 0.0
    cmeout[:ncol, :pver] = 0.0

    am_evp_st[:ncol, :pver] = 0.0

    # For refl calc
    rainrt1 = np.zeros((ncol, pver))

    # Initialize precip fraction and output tendencies
    cldmax = np.full((ncol, pver), mincld)

    # Initialize aerosol number
    dum2l = np.zeros((ncol, pver))
    dum2i = np.zeros((ncol, pver))

    # Initialize avg precip rate
    prect1 = np.zeros(ncol)
    preci1 = np.zeros(ncol)

    # Get humidity and saturation vapor pressures
    es = np.zeros(ncol)
    qs = np.zeros(ncol)
    esl = np.zeros((ncol, pver))
    esi = np.zeros((ncol, pver))
    relhum = np.zeros((ncol, pver))

    for k in range(top_lev-1, pver):
        for i in range(ncol):
            es[i] = svp_water(t[i, k])
            qs[i] = svp_to_qsat(es[i], p[i, k])

            # Prevents negative values
            if qs[i] < 0.0:
                qs[i] = 1.0
                es[i] = p[i, k]

            esl[i, k] = svp_water(t[i, k])
            esi[i, k] = svp_ice(t[i, k])

            # Make sure when above freezing that esi=esl, not active yet
            if t[i, k] > tmelt:
                esi[i, k] = esl[i, k]

            relhum[i, k] = q[i, k] / qs[i]

            # Get cloud fraction, check for minimum
            cldm = max(cldn[i, k], mincld)
            cldmw = max(cldn[i, k], mincld)

            icldm = max(icecldf[i, k], mincld)
            lcldm = max(liqcldf[i, k], mincld)

            # Subcolumns, set cloud fraction variables to one
            if microp_uniform:
                cldm = mincld
                cldmw = mincld
                icldm = mincld
                lcldm = mincld

                if qc[i, k] >= qsmall:
                    lcldm = 1.0
                    cldm = 1.0
                    cldmw = 1.0

                if qi[i, k] >= qsmall:
                    cldm = 1.0
                    icldm = 1.0

            # Calculate nfice based on liquid and ice mmr (no rain and snow mmr available yet)
            nfice[i, k] = 0.0
            dumfice = qc[i, k] + qi[i, k]
            if dumfice > qsmall and qi[i, k] > qsmall:
                nfice[i, k] = qi[i, k] / dumfice

            if do_cldice and (t[i, k] < tmelt - 5.0):
                # If aerosols interact with ice set number of activated ice nuclei
                dum2 = naai[i, k]

                dumnnuc = (dum2 - ni[i, k] / icldm) / deltat * icldm
                dumnnuc = max(dumnnuc, 0.0)
                # Get provisional ni and qi after nucleation in order to calculate
                # Bergeron process below
                ninew = ni[i, k] + dumnnuc * deltat
                qinew = qi[i, k] + dumnnuc * deltat * mi0

            else:
                ninew = ni[i, k]
                qinew = qi[i, k]

            # Initialize CME components
            cme = 0.0
            cmei = 0.0

            # Bergeron process
            # Make sure to initialize bergeron process to zero
            berg = 0.0
            prd = 0.0

            # Condensation loop
            # Get in-cloud qi and ni after nucleation
            if icldm > 0.0:
                qiic = qinew / icldm
                niic = ninew / icldm
            else:
                qiic = 0.0
                niic = 0.0

            if nicons:
                niic = ninst / rho[i, k]

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
                    if qiic >= qsmall:
                        lami = (cons1 * ci * niic / qiic)**(1.0 / di)
                        n0i = niic * lami

                        # Check for slope
                        # Adjust vars
                        if lami < lammini:
                            lami = lammini
                            n0i = lami**(di + 1.0) * qiic / (ci * cons1)
                        elif lami > lammaxi:
                            lami = lammaxi
                            n0i = lami**(di + 1.0) * qiic / (ci * cons1)

                        epsi = 2.0 * pi * n0i * rho[i, k] * dv[i, k] / (lami * lami)

                        # If liquid exists
                        if qc[i, k] > qsmall:
                            # Begin bergeron process
                            # Calculate Bergeron process
                            prd = epsi * (qvl - qvi) / abi
                        else:
                            prd = 0.0

                        # Multiply by cloud fraction
                        prd = prd * min(icldm, lcldm)

                        # Transfer of existing cloud liquid to ice
                        berg = max(0.0, prd)

                    if berg > 0.0:
                        bergtsf = max(0.0, (qc[i, k] / berg) / deltat)

                        if bergtsf < 1.0:
                            berg = max(0.0, qc[i, k] / deltat)

                if bergtsf < 1.0 or icldm > lcldm:
                    if qiic >= qsmall:
                        # First case is for case when liquid water is present, but is completely depleted
                        # in time step, i.e., bergrsf > 0 but < 1
                        if qc[i, k] >= qsmall:
                            rhin = (1.0 + relhum[i, k]) / 2.0
                            if (rhin * esl[i, k] / esi[i, k]) > 1.0:
                                prd = epsi * (rhin * qvl - qvi) / abi

                                # Multiply by cloud fraction assuming liquid/ice maximum overlap
                                prd = prd * min(icldm, lcldm)

                                # Add to cmei
                                cmei += (prd * (1.0 - bergtsf))

                        # Second case is for pure ice cloud, either no liquid, or icldm > lcldm
                        if qc[i, k] < qsmall or icldm > lcldm:
                            # Note: for case of no liquid, need to set liquid cloud fraction to zero
                            # Store liquid cloud fraction in 'dum'
                            if qc[i, k] < qsmall:
                                dum = 0.0
                            else:
                                dum = lcldm

                            # Set RH to grid-mean value for pure ice cloud
                            rhin = relhum[i, k]

                            if (rhin * esl[i, k] / esi[i, k]) > 1.0:
                                prd = epsi * (rhin * qvl - qvi) / abi

                                # Multiply by relevant cloud fraction for pure ice cloud
                                # Assuming maximum overlap of liquid/ice
                                prd = prd * max((icldm - dum), 0.0)
                                cmei += prd

                # If deposition, it should not reduce grid mean rhi below 1.0
                if cmei > 0.0 and (relhum[i, k] * esl[i, k] / esi[i, k]) > 1.0:
                    cmei = min(cmei, (q[i, k] - qs[i] * esi[i, k] / esl[i, k]) / abi / deltat)

            # Evaporation should not exceed available water
            if (-berg) < -qc[i, k] / deltat:
                berg = max(qc[i, k] / deltat, 0.0)

            # Sublimation process
            if do_cldice and ((relhum[i, k] * esl[i, k] / esi[i, k]) < 1.0 and qiic >= qsmall):
                qvi = svp_to_qsat(esi[i, k], p[i, k])
                qvl = svp_to_qsat(esl[i, k], p[i, k])
                dqsidt = xxls * qvi / (rv * t[i, k]**2)
                abi = 1.0 + dqsidt * xxls / cpp

                # Get ice size distribution parameters
                lami = (cons1 * ci * niic / qiic)**(1.0 / di)
                n0i = niic * lami

                # Check for slope
                # Adjust vars
                if lami < lammini:
                    lami = lammini
                    n0i = lami**(di + 1.0) * qiic / (ci * cons1)
                elif lami > lammaxi:
                    lami = lammaxi
                    n0i = lami**(di + 1.0) * qiic / (ci * cons1)

                epsi = 2.0 * pi * n0i * rho[i, k] * dv[i, k] / (lami * lami)

                # Modify for ice fraction below
                prd = epsi * (relhum[i, k] * qvl - qvi) / abi * icldm
                cmei = min(prd, 0.0)

            # Sublimation should not exceed available ice
            if cmei < -qi[i, k] / deltat:
                cmei = -qi[i, k] / deltat

            # Sublimation should not increase grid mean rhi above 1.0
            if cmei < 0.0 and (relhum[i, k] * esl[i, k] / esi[i, k]) < 1.0:
                cmei = min(0.0, max(cmei, (q[i, k] - qs[i] * esi[i, k] / esl[i, k]) / abi / deltat))

            # Limit cmei due for roundoff error
            cmei *= omsm

            # Conditional for ice nucleation
            if do_cldice and (t[i, k] < (tmelt - 5.0)):
                # Using Liu et al. (2007) ice nucleation with hooks into simulated aerosol
                # Ice nucleation rate (dum2) has already been calculated and read in (naai)
                dum2i[i, k] = naai[i, k]
            else:
                dum2i[i, k] = 0.0

    # Initialize sub-step precip flux variables
    rflx1 = np.zeros((ncol, pver + 1))
    sflx1 = np.zeros((ncol, pver + 1))

    # Initialize final precip flux variables
    rflx[:, 0] = 0.0
    sflx[:, 0] = 0.0
    rflx[:, 1:] = 0.0
    sflx[:, 1:] = 0.0

    ltrue = np.zeros(ncol, dtype=bool)
    for i in range(ncol):
        for k in range(top_lev - 1, pver):
            if qc[i, k] >= qsmall or qi[i, k] >= qsmall or cmei[i, k] >= qsmall:
                ltrue[i] = True

    # Assign number of sub-steps to iter
    # Use 2 sub-steps, following tests described in MG2008
    iter = 2

    # Get sub-step time step
    deltat /= iter

    # Since activation/nucleation processes are fast, need to take into account
    # factor mtime = mixing timescale in cloud / model time step
    # Mixing time can be interpreted as cloud depth divided by sub-grid vertical velocity
    # For now mixing timescale is assumed to be 1 timestep for modal aerosols, 20 min bulk

    # Note: mtime for bulk aerosols was set to: mtime=deltat/1200.0
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

        qcsinksum_rate1ord = np.zeros(pver)
        qcsum_rate1ord = np.zeros(pver)

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

            # Begin new i,k loop, calculate new cldmax after adjustment to cldm above

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

            for k in range(top_lev - 1, pver):
                qcvar = relvar[i, k]
                cons2 = math.gamma(qcvar + 2.47)
                cons3 = math.gamma(qcvar)
                cons9 = math.gamma(qcvar + 2.0)
                cons10 = math.gamma(qcvar + 1.0)
                cons12 = math.gamma(qcvar + 1.15)
                cons15 = math.gamma(qcvar + bc / 3.0)
                cons18 = qcvar**2.47
                cons19 = qcvar**2
                cons20 = qcvar**1.15

                # Set cwml and cwmi to current qc and qi
                cwml = qc[i, k]
                cwmi = qi[i, k]

                # Initialize precip fallspeeds to zero
                ums = 0.0
                uns = 0.0
                umr = 0.0
                unr = 0.0

                # Calculate precip fraction based on maximum overlap assumption
                # For sub-columns cldm has already been set to 1 if cloud
                # water or ice is present, so cldmax will be correctly set below
                # and nothing extra needs to be done here

                if k == top_lev - 1:
                    cldmax[i, k] = cldm
                else:
                    # If rain or snow mix ratio is smaller than
                    # threshold, then set cldmax to cloud fraction at current level
                    if do_clubb_sgs:
                        if qc[i, k] >= qsmall or qi[i, k] >= qsmall:
                            cldmax[i, k] = cldm
                        else:
                            cldmax[i, k] = cldmax[i, k - 1]
                    else:
                        if qric[i, k - 1] >= qsmall or qniic[i, k - 1] >= qsmall:
                            cldmax[i, k] = max(cldmax[i, k - 1], cldm)
                        else:
                            cldmax[i, k] = cldm

                # Decrease in number concentration due to sublimation/evap
                # Divide by cloud fraction to get in-cloud decrease
                # Don't reduce Nc due to bergeron process

                if cmei[i, k] < 0.0 and qi[i, k] > qsmall and cldm > mincld:
                    nsubi = cmei[i, k] / qi[i, k] * ni[i, k] / cldm
                else:
                    nsubi = 0.0
                nsubc = 0.0

                # Ice nucleation if activated nuclei exist at t<-5C AND rhmini + 5%
                if do_cldice and dum2i[i, k] > 0.0 and t[i, k] < (tmelt - 5.0) and \
                   relhum[i, k] * esl[i, k] / esi[i, k] > rhmini + 0.05:

                    # If NCAI > 0. then set numice = ncai (as before)
                    # Note: this is gridbox averaged
                    nnuccd = (dum2i[i, k] - ni[i, k] / icldm) / deltat * icldm
                    nnuccd = max(nnuccd, 0.0)
                    nimax = dum2i[i, k] * icldm

                    # Calc mass of new particles using new crystal mass...
                    # Also this will be multiplied by mtime as nnuccd is...
                    mnuccd = nnuccd * mi0

                    # Add mnuccd to cmei....
                    cmei[i, k] += mnuccd * mtime

                    # Limit cmei
                    qvi = svp_to_qsat(esi[i, k], p[i, k])
                    dqsidt = xxls * qvi / (rv * t[i, k]**2)
                    abi = 1.0 + dqsidt * xxls / cpp
                    cmei[i, k] = min(cmei[i, k], (q[i, k] - qvi) / abi / deltat)

                    # Limit for roundoff error
                    cmei[i, k] *= omsm

                else:
                    nnuccd = 0.0
                    nimax = 0.0
                    mnuccd = 0.0

                # Obtain in-cloud values of cloud water/ice mixing ratios and number concentrations
                # for microphysical process calculations
                # Units are kg/kg for mixing ratio, 1/kg for number conc

                # Limit in-cloud values to 0.005 kg/kg
                qcic = min(cwml / lcldm, 5e-3)
                qiic = min(cwmi / icldm, 5e-3)
                ncic = max(nc[i, k] / lcldm, 0.0)
                niic = max(ni[i, k] / icldm, 0.0)


                if nccons:
                    ncic = ncnst / rho[i, k]
                if nicons:
                    niic = ninst / rho[i, k]

                if qc[i, k] - berg[i, k] * deltat < qsmall:
                    qcic = 0.0
                    ncic = 0.0
                    if qc[i, k] - berg[i, k] * deltat < 0.0:
                        berg[i, k] = qc[i, k] / deltat * omsm

                if do_cldice and qi[i, k] + (cmei[i, k] + berg[i, k]) * deltat < qsmall:
                    qiic = 0.0
                    niic = 0.0
                    if qi[i, k] + (cmei[i, k] + berg[i, k]) * deltat < 0.0:
                        cmei[i, k] = (-qi[i, k] / deltat - berg[i, k]) * omsm

                # Add to cme output
                cmeout[i, k] += cmei[i, k]

                # Droplet activation
                # Calculate potential for droplet activation if cloud water is present
                # Formulation from Abdul-Razzak and Ghan (2000) and Abdul-Razzak et al. (1998), AR98
                # Number tendency (npccnin) is read in from companion routine

                # Assume aerosols already activated are equal to number of existing droplets for simplicity
                # Multiply by cloud fraction to obtain grid-average tendency

                if qcic >= qsmall:
                    npccn = max(0.0, npccnin[i, k])
                    dum2l[i, k] = (nc[i, k] + npccn * deltat * mtime) / lcldm
                    dum2l[i, k] = max(dum2l[i, k], cdnl / rho[i, k])  # sghan minimum in #/cm3
                    ncmax = dum2l[i, k] * lcldm
                else:
                    npccn = 0.0
                    dum2l[i, k] = 0.0
                    ncmax = 0.0

                # Get size distribution parameters based on in-cloud cloud water/ice
                # These calculations also ensure consistency between number and mixing ratio

                # Cloud ice
                if qiic >= qsmall:
                    # Add upper limit to in-cloud number concentration to prevent numerical error
                    niic = min(niic, qiic * 1e20)

                    lami = (cons1 * ci * niic / qiic)**(1.0 / di)
                    n0i = niic * lami

                    # Check for slope
                    # Adjust vars
                    if lami < lammini:
                        lami = lammini
                        n0i = lami**(di + 1.0) * qiic / (ci * cons1)
                        niic = n0i / lami
                    elif lami > lammaxi:
                        lami = lammaxi
                        n0i = lami**(di + 1.0) * qiic / (ci * cons1)
                        niic = n0i / lami
                else:
                    lami = 0.0
                    n0i = 0.0

                if qcic >= qsmall:
                    # Add upper limit to in-cloud number concentration to prevent numerical error
                    ncic = min(ncic, qcic * 1e20)

                    ncic = max(ncic, cdnl / rho[i, k])  # sghan minimum in #/cm

                    # Get pgam from fit to observations of martin et al. 1994
                    pgam = 0.0005714 * (ncic / 1e6 * rho[i, k]) + 0.2714
                    pgam = 1.0 / (pgam**2) - 1.0
                    pgam = max(pgam, 2.0)
                    pgam = min(pgam, 15.0)

                    # Calculate lamc
                    lamc = (pi / 6.0 * rhow * ncic * math.gamma(pgam + 4.0) /
                            (qcic * math.gamma(pgam + 1.0)))**(1.0 / 3.0)

                    # lammin, 50 micron diameter max mean size
                    lammin = (pgam + 1.0) / 50e-6
                    lammax = (pgam + 1.0) / 2e-6

                    if lamc < lammin:
                        lamc = lammin
                        ncic = 6.0 * lamc**3 * qcic * math.gamma(pgam + 1.0) / (pi * rhow * math.gamma(pgam + 4.0))
                    elif lamc > lammax:
                        lamc = lammax
                        ncic = 6.0 * lamc**3 * qcic * math.gamma(pgam + 1.0) / (pi * rhow * math.gamma(pgam + 4.0))

                    # Parameter to calculate droplet freezing
                    cdist1 = ncic / math.gamma(pgam + 1.0)

                else:
                    lamc = 0.0
                    cdist1 = 0.0

                # Begin microphysical process calculations
                # Autoconversion of cloud liquid water to rain
                # Formula from Khrouditnov and Kogan (2000), modified for sub-grid distribution of qc
                # Minimum qc of 1 x 10^-8 prevents floating point error

                if qcic >= 1e-8:
                    # nprc is increase in rain number conc due to autoconversion
                    # nprc1 is decrease in cloud droplet conc due to autoconversion

                    # Assume exponential sub-grid distribution of qc, resulting in additional
                    # factor related to qcvar below

                    # hm switch for sub-columns, don't include sub-grid qc
                    if microp_uniform:
                        prc = 1350.0 * qcic**2.47 * (ncic / 1e6 * rho[i, k])**(-1.79)
                        nprc = prc / (4.0 / 3.0 * pi * rhow * (25e-6)**3)
                        nprc1 = prc / (qcic / ncic)
                    else:
                        prc = cons2 / (cons3 * cons18) * 1350.0 * qcic**2.47 * (ncic / 1e6 * rho[i, k])**(-1.79)
                        nprc = prc / cons22
                        nprc1 = prc / (qcic / ncic)
                else:
                    prc = 0.0
                    nprc = 0.0
                    nprc1 = 0.0

                # Add autoconversion to precip from above to get provisional rain mixing ratio
                # and number concentration (qric and nric)

                # 0.45 m/s is fallspeed of new rain drop (80 micron diameter)
                dum = 0.45
                dum1 = 0.45

                if k == top_lev - 1:
                    qric[i, k] = prc * lcldm * dz[i, k] / cldmax[i, k] / dum
                    nric[i, k] = nprc * lcldm * dz[i, k] / cldmax[i, k] / dum
                else:
                    if qric[i, k - 1] >= qsmall:
                        dum = umr[k - 1]
                        dum1 = unr[k - 1]

                    # No autoconversion of rain number if rain/snow falling from above
                    # This assumes that new drizzle drops formed by autoconversion are rapidly collected
                    # by the existing rain/snow particles from above
                    if qric[i, k - 1] >= 1e-9 or qniic[i, k - 1] >= 1e-9:
                        nprc = 0.0

                    qric[i, k] = (rho[i, k - 1] * umr[k - 1] * qric[i, k - 1] * cldmax[i, k - 1] +
                                  (rho[i, k] * dz[i, k] * ((pra[k - 1] + prc) * lcldm +
                                   (pre[k - 1] - pracs[k - 1] - mnuccr[k - 1]) * cldmax[i, k]))) / \
                                 (dum * rho[i, k] * cldmax[i, k])
                    nric[i, k] = (rho[i, k - 1] * unr[k - 1] * nric[i, k - 1] * cldmax[i, k - 1] +
                                  (rho[i, k] * dz[i, k] * (nprc * lcldm +
                                   (nsubr[k - 1] - npracs[k - 1] - nnuccr[k - 1] + nragg[k - 1]) * cldmax[i, k]))) / \
                                 (dum1 * rho[i, k] * cldmax[i, k])

                # Autoconversion of cloud ice to snow
                # Similar to Ferrier (1994)
                if do_cldice:
                    if t[i, k] <= 273.15 and qiic >= qsmall:
                        # Note: assumes autoconversion timescale of 180 sec
                        nprci = n0i / (lami * 180.0) * np.exp(-lami * dcs)

                        prci = pi * rhoi * n0i / (6.0 * 180.0) * \
                               (cons23 / lami + 3.0 * cons24 / lami**2 +
                                6.0 * dcs / lami**3 + 6.0 / lami**4) * np.exp(-lami * dcs)
                    else:
                        prci = 0.0
                        nprci = 0.0
                else:
                    # Add in the particles that we have already converted to snow, and
                    # don't do any further autoconversion of ice.
                    prci = tnd_qsnow[i, k] / cldm
                    nprci = tnd_nsnow[i, k] / cldm

                # Add autoconversion to flux from level above to get provisional snow mixing ratio
                # and number concentration (qniic and nsic)
                dum = (asn[i, k] * cons25)
                dum1 = (asn[i, k] * cons25)

                if k == top_lev - 1:
                    qniic[i, k] = prci * icldm * dz[i, k] / cldmax[i, k] / dum
                    nsic[i, k] = nprci * icldm * dz[i, k] / cldmax[i, k] / dum
                else:
                    if qniic[i, k - 1] >= qsmall:
                        dum = ums[k - 1]
                        dum1 = uns[k - 1]

                    qniic[i, k] = (rho[i, k - 1] * ums[k - 1] * qniic[i, k - 1] * cldmax[i, k - 1] +
                                   (rho[i, k] * dz[i, k] * ((prci + prai[k - 1] + psacws[k - 1] + bergs[k - 1]) * icldm +
                                    (prds[k - 1] + pracs[k - 1] + mnuccr[k - 1]) * cldmax[i, k]))) / \
                                  (dum * rho[i, k] * cldmax[i, k])

                    nsic[i, k] = (rho[i, k - 1] * uns[k - 1] * nsic[i, k - 1] * cldmax[i, k - 1] +
                                  (rho[i, k] * dz[i, k] * (nprci * icldm +
                                   (nsubs[k - 1] + nsagg[k - 1] + nnuccr[k - 1]) * cldmax[i, k]))) / \
                                 (dum1 * rho[i, k] * cldmax[i, k])

                # If precip mix ratio is zero so should number concentration
                if qniic[i, k] < qsmall:
                    qniic[i, k] = 0.0
                    nsic[i, k] = 0.0

                if qric[i, k] < qsmall:
                    qric[i, k] = 0.0
                    nric[i, k] = 0.0

                # Make sure number concentration is a positive number to avoid
                # taking root of negative later
                nric[i, k] = max(nric[i, k], 0.0)
                nsic[i, k] = max(nsic[i, k], 0.0)

                # Get size distribution parameters for precip
                # Rain
                if qric[i, k] >= qsmall:
                    lamr = (pi * rhow * nric[i, k] / qric[i, k])**(1.0 / 3.0)
                    n0r = nric[i, k] * lamr

                    # Check for slope
                    # Adjust vars
                    if lamr < lamminr:
                        lamr = lamminr
                        n0r = lamr**4 * qric[i, k] / (pi * rhow)
                        nric[i, k] = n0r / lamr
                    elif lamr > lammaxr:
                        lamr = lammaxr
                        n0r = lamr**4 * qric[i, k] / (pi * rhow)
                        nric[i, k] = n0r / lamr

                    # Provisional rain number and mass weighted mean fallspeed (m/s)
                    unr = min(arn[i, k] * cons4 / lamr**br, 9.1 * rhof[i, k])
                    umr = min(arn[i, k] * cons5 / (6.0 * lamr**br), 9.1 * rhof[i, k])

                else:
                    lamr = 0.0
                    n0r = 0.0
                    umr = 0.0
                    unr = 0.0

                # Snow
                if qniic[i, k] >= qsmall:
                    lams = (cons6 * cs * nsic[i, k] / qniic[i, k])**(1.0 / ds)
                    n0s = nsic[i, k] * lams

                    # Check for slope
                    # Adjust vars
                    if lams < lammins:
                        lams = lammins
                        n0s = lams**(ds + 1.0) * qniic[i, k] / (cs * cons6)
                        nsic[i, k] = n0s / lams
                    elif lams > lammaxs:
                        lams = lammaxs
                        n0s = lams**(ds + 1.0) * qniic[i, k] / (cs * cons6)
                        nsic[i, k] = n0s / lams

                    # Provisional snow number and mass weighted mean fallspeed (m/s)
                    ums = min(asn[i, k] * cons8 / (6.0 * lams**bs), 1.2 * rhof[i, k])
                    uns = min(asn[i, k] * cons7 / lams**bs, 1.2 * rhof[i, k])

                else:
                    lams = 0.0
                    n0s = 0.0
                    ums = 0.0
                    uns = 0.0

                # Heterogeneous freezing of cloud water
                if not use_hetfrz_classnuc:
                    if do_cldice and qcic >= qsmall and t[i, k] < 269.15:
                        # Immersion freezing (Bigg, 1953)
                        # Subcolumns
                        if microp_uniform:
                            mnuccc = (pi * pi / 36.0 * rhow *
                                      cdist1 * math.gamma(7.0 + pgam) *
                                      bimm * (np.exp(aimm * (273.15 - t[i, k])) - 1.0) /
                                      lamc**3 / lamc**3)

                            nnuccc = (pi / 6.0 * cdist1 * math.gamma(pgam + 4.0) *
                                      bimm *
                                      (np.exp(aimm * (273.15 - t[i, k])) - 1.0) / lamc**3)
                        else:
                            mnuccc = cons9 / (cons3 * cons19) * \
                                     (pi * pi / 36.0 * rhow *
                                      cdist1 * math.gamma(7.0 + pgam) *
                                      bimm * (np.exp(aimm * (273.15 - t[i, k])) - 1.0) /
                                      lamc**3 / lamc**3)

                            nnuccc = cons10 / (cons3 * qcvar) * \
                                     (pi / 6.0 * cdist1 * math.gamma(pgam + 4.0) *
                                      bimm *
                                      (np.exp(aimm * (273.15 - t[i, k])) - 1.0) / lamc**3)

                        # Contact freezing (-40<T<-3 C) (Young, 1974) with hooks into simulated dust
                        # Dust size and number in 4 bins are read in from companion routine
                        tcnt = (270.16 - t[i, k])**1.3
                        viscosity = 1.8e-5 * (t[i, k] / 298.0)**0.85  # Viscosity (kg/m/s)
                        mfp = 2.0 * viscosity / (p[i, k] *
                                                 np.sqrt(8.0 * 28.96e-3 / (pi * 8.314409 * t[i, k])))  # Mean free path (m)

                        nslip1 = 1.0 + (mfp / rndst[i, k, 0]) * (1.257 + (0.4 * np.exp(-(1.1 * rndst[i, k, 0] / mfp))))
                        nslip2 = 1.0 + (mfp / rndst[i, k, 1]) * (1.257 + (0.4 * np.exp(-(1.1 * rndst[i, k, 1] / mfp))))
                        nslip3 = 1.0 + (mfp / rndst[i, k, 2]) * (1.257 + (0.4 * np.exp(-(1.1 * rndst[i, k, 2] / mfp))))
                        nslip4 = 1.0 + (mfp / rndst[i, k, 3]) * (1.257 + (0.4 * np.exp(-(1.1 * rndst[i, k, 3] / mfp))))

                        ndfaer1 = 1.381e-23 * t[i, k] * nslip1 / (6.0 * pi * viscosity * rndst[i, k, 0])  # Aerosol diffusivity (m2/s)
                        ndfaer2 = 1.381e-23 * t[i, k] * nslip2 / (6.0 * pi * viscosity * rndst[i, k, 1])
                        ndfaer3 = 1.381e-23 * t[i, k] * nslip3 / (6.0 * pi * viscosity * rndst[i, k, 2])
                        ndfaer4 = 1.381e-23 * t[i, k] * nslip4 / (6.0 * pi * viscosity * rndst[i, k, 3])

                        if microp_uniform:
                            mnucct = ((ndfaer1 * (nacon[i, k, 0] * tcnt) + ndfaer2 * (nacon[i, k, 1] * tcnt) +
                                       ndfaer3 * (nacon[i, k, 2] * tcnt) + ndfaer4 * (nacon[i, k, 3] * tcnt)) *
                                      pi * pi / 3.0 * rhow *
                                      cdist1 * math.gamma(pgam + 5.0) / lamc**4)

                            nnucct = ((ndfaer1 * (nacon[i, k, 0] * tcnt) + ndfaer2 * (nacon[i, k, 1] * tcnt) +
                                       ndfaer3 * (nacon[i, k, 2] * tcnt) + ndfaer4 * (nacon[i, k, 3] * tcnt)) *
                                      2.0 * pi *
                                      cdist1 * math.gamma(pgam + 2.0) / lamc)
                        else:
                            mnucct = math.gamma(qcvar + 4.0 / 3.0) / (cons3 * qcvar**(4.0 / 3.0)) * \
                                     ((ndfaer1 * (nacon[i, k, 0] * tcnt) + ndfaer2 * (nacon[i, k, 1] * tcnt) +
                                       ndfaer3 * (nacon[i, k, 2] * tcnt) + ndfaer4 * (nacon[i, k, 3] * tcnt)) *
                                      pi * pi / 3.0 * rhow *
                                      cdist1 * math.gamma(pgam + 5.0) / lamc**4)

                            nnucct = math.gamma(qcvar + 1.0 / 3.0) / (cons3 * qcvar**(1.0 / 3.0)) * \
                                     ((ndfaer1 * (nacon[i, k, 0] * tcnt) + ndfaer2 * (nacon[i, k, 1] * tcnt) +
                                       ndfaer3 * (nacon[i, k, 2] * tcnt) + ndfaer4 * (nacon[i, k, 3] * tcnt)) *
                                      2.0 * pi *
                                      cdist1 * math.gamma(pgam + 2.0) / lamc)

                        # Make sure number of droplets frozen does not exceed available ice nuclei concentration
                        # This prevents 'runaway' droplet freezing
                        if nnuccc * lcldm > nnuccd:
                            dum = (nnuccd / (nnuccc * lcldm))
                            # Scale mixing ratio of droplet freezing with limit
                            mnuccc *= dum
                            nnuccc = nnuccd / lcldm

                    else:
                        mnuccc = 0.0
                        nnuccc = 0.0
                        mnucct = 0.0
                        nnucct = 0.0

                else:
                    if do_cldice and qcic >= qsmall:
                        con1 = 1.0 / (1.333 * pi)**0.333
                        r3lx = con1 * (rho[i, k] * qcic / (rhow * max(ncic * rho[i, k], 1.0e6)))**0.333  # in m
                        r3lx = max(4.e-6, r3lx)
                        mi0l = 4.0 / 3.0 * pi * rhow * r3lx**3

                        nnuccc = frzimm[i, k] * 1.0e6 / rho[i, k]
                        mnuccc = nnuccc * mi0l

                        nnucct = frzcnt[i, k] * 1.0e6 / rho[i, k]
                        mnucct = nnucct * mi0l

                        nnudep = frzdep[i, k] * 1.0e6 / rho[i, k]
                        mnudep = nnudep * mi0
                    else:
                        nnuccc = 0.0
                        mnuccc = 0.0

                        nnucct = 0.0
                        mnucct = 0.0

                        nnudep = 0.0
                        mnudep = 0.0

                # Snow self-aggregation from passarelli, 1978, used by reisner, 1998
                # This is hard-wired for bs = 0.4 for now
                # Ignore self-collection of cloud ice
                if qniic[i, k] >= qsmall and t[i, k] <= 273.15:
                    nsagg = -1108.0 * asn[i, k] * Eii * \
                            pi**((1.0 - bs) / 3.0) * rhosn**((-2.0 - bs) / 3.0) * rho[i, k]** \
                            ((2.0 + bs) / 3.0) * qniic[i, k]**((2.0 + bs) / 3.0) * \
                            (nsic[i, k] * rho[i, k])**((4.0 - bs) / 3.0) / \
                            (4.0 * 720.0 * rho[i, k])
                else:
                    nsagg = 0.0

                # Accretion of cloud droplets onto snow/graupel
                # Here use continuous collection equation with
                # simple gravitational collection kernel
                # Ignore collisions between droplets/cloud ice
                # since minimum size ice particle for accretion is 50 - 150 micron

                # Ignore collision of snow with droplets above freezing
                if qniic[i, k] >= qsmall and t[i, k] <= tmelt and qcic[i, k] >= qsmall:
                    # Put in size dependent collection efficiency
                    # Mean diameter of snow is area-weighted, since
                    # accretion is function of crystal geometric area
                    # Collection efficiency is approximation based on stoke's law (Thompson et al. 2004)
                    dc0 = (pgam + 1.0) / lamc
                    ds0 = 1.0 / lams
                    dum = dc0 * dc0 * uns * rhow / (9.0 * mu[i, k] * ds0)
                    eci = dum * dum / ((dum + 0.4) * (dum + 0.4))

                    eci = max(eci, 0.0)
                    eci = min(eci, 1.0)

                    # No impact of sub-grid distribution of qc since psacws
                    # is linear in qc
                    psacws = pi / 4.0 * asn[i, k] * qcic * rho[i, k] * \
                             n0s * eci * cons11 / \
                             lams**(bs + 3.0)
                    npsacws = pi / 4.0 * asn[i, k] * ncic * rho[i, k] * \
                              n0s * eci * cons11 / \
                              lams**(bs + 3.0)
                else:
                    psacws = 0.0
                    npsacws = 0.0

                # Add secondary ice production due to accretion of droplets by snow
                # (Hallet-Mossop process) (from Cotton et al., 1986)
                if not do_cldice:
                    ni_secp = 0.0
                    nsacwi = 0.0
                    msacwi = 0.0
                elif (t[i, k] < 270.16) and (t[i, k] >= 268.16):
                    ni_secp = 3.5e8 * (270.16 - t[i, k]) / 2.0 * psacws
                    nsacwi = ni_secp
                    msacwi = min(ni_secp * mi0, psacws)
                elif (t[i, k] < 268.16) and (t[i, k] >= 265.16):
                    ni_secp = 3.5e8 * (t[i, k] - 265.16) / 3.0 * psacws
                    nsacwi = ni_secp
                    msacwi = min(ni_secp * mi0, psacws)
                else:
                    ni_secp = 0.0
                    nsacwi = 0.0
                    msacwi = 0.0
                psacws = max(0.0, psacws - ni_secp * mi0)

                # Accretion of rain water by snow
                # Formula from ikawa and saito, 1991, used by reisner et al., 1998
                if qric[i, k] >= 1e-8 and qniic[i, k] >= 1e-8 and t[i, k] <= 273.15:
                    pracs = pi * pi * ecr * (((1.2 * umr - 0.95 * ums)**2 +
                                              0.08 * ums * umr)**0.5 * rhow * rho[i, k] *
                                             n0r * n0s *
                                             (5.0 / (lamr**6 * lams) +
                                              2.0 / (lamr**5 * lams**2) +
                                              0.5 / (lamr**4 * lams**3)))

                    npracs = pi / 2.0 * rho[i, k] * ecr * (1.7 * (unr - uns)**2 +
                                                           0.3 * unr * uns)**0.5 * n0r * n0s * \
                             (1.0 / (lamr**3 * lams) +
                              1.0 / (lamr**2 * lams**2) +
                              1.0 / (lamr * lams**3))
                else:
                    pracs = 0.0
                    npracs = 0.0

                # Heterogeneous freezing of rain drops
                # Follows from Bigg (1953)
                if t[i, k] < 269.15 and qric[i, k] >= qsmall:
                    mnuccr = 20.0 * pi * pi * rhow * nric[i, k] * bimm * \
                             (np.exp(aimm * (273.15 - t[i, k])) - 1.0) / lamr**3 / \
                             lamr**3

                    nnuccr = pi * nric[i, k] * bimm * \
                             (np.exp(aimm * (273.15 - t[i, k])) - 1.0) / lamr**3
                else:
                    mnuccr = 0.0
                    nnuccr = 0.0

                # Accretion of cloud liquid water by rain
                # Formula from Khrouditnov and Kogan (2000)
                # Gravitational collection kernel, droplet fall speed neglected
                if qric[i, k] >= qsmall and qcic[i, k] >= qsmall:
                    # Include sub-grid distribution of cloud water
                    # Add sub-column switch
                    if microp_uniform:
                        pra = 67.0 * (qcic * qric[i, k])**1.15
                        npra = pra / (qcic / ncic)
                    else:
                        pra = accre_enhan[i, k] * (cons12 / (cons3 * cons20) * 67.0 * (qcic * qric[i, k])**1.15)
                        npra = pra / (qcic / ncic)
                else:
                    pra = 0.0
                    npra = 0.0

                # Self-collection of rain drops
                # From Beheng(1994)
                if qric[i, k] >= qsmall:
                    nragg = -8.0 * nric[i, k] * qric[i, k] * rho[i, k]
                else:
                    nragg = 0.0

                # Accretion of cloud ice by snow
                # For this calculation, it is assumed that the Vs >> Vi
                # and Ds >> Di for continuous collection
                if do_cldice and qniic[i, k] >= qsmall and qiic >= qsmall and t[i, k] <= 273.15:
                    prai = pi / 4.0 * asn[i, k] * qiic * rho[i, k] * \
                           n0s * Eii * cons11 / \
                           lams**(bs + 3.0)
                    nprai = pi / 4.0 * asn[i, k] * niic * \
                            rho[i, k] * n0s * Eii * cons11 / \
                            lams**(bs + 3.0)
                else:
                    prai = 0.0
                    nprai = 0.0

                # Calculate evaporation/sublimation of rain and snow
                # Note: evaporation/sublimation occurs only in cloud-free portion of grid cell
                # In-cloud condensation/deposition of rain and snow is neglected
                # except for transfer of cloud water to snow through bergeron process

                # Initialize evap/sub tendencies
                pre = 0.0
                prds = 0.0

                # Evaporation of rain
                # Only calculate if there is some precip fraction > cloud fraction
                if qcic + qiic < 1e-6 or cldmax[i, k] > lcldm:
                    # Set temporary cloud fraction to zero if cloud water + ice is very small
                    # This will ensure that evaporation/sublimation of precip occurs over
                    # entire grid cell, since min cloud fraction is specified otherwise
                    if qcic + qiic < 1e-6:
                        dum = 0.0
                    else:
                        dum = lcldm

                    # Saturation vapor pressure
                    esn = svp_water(t[i, k])
                    qsn = svp_to_qsat(esn, p[i, k])

                    # Recalculate saturation vapor pressure for liquid and ice
                    esl[i, k] = esn
                    esi[i, k] = svp_ice(t[i, k])
                    # hm fix, make sure when above freezing that esi=esl, not active yet
                    if t[i, k] > tmelt:
                        esi[i, k] = esl[i, k]

                    # Calculate q for out-of-cloud region
                    qclr = (q[i, k] - dum * qsn) / (1.0 - dum)

                    if qric[i, k] >= qsmall:
                        qvs = svp_to_qsat(esl[i, k], p[i, k])
                        dqsdt = xxlv * qvs / (rv * t[i, k]**2)
                        ab = 1.0 + dqsdt * xxlv / cpp
                        epsr = 2.0 * pi * n0r * rho[i, k] * Dv[i, k] * \
                               (f1r / (lamr * lamr) +
                                f2r * (arn[i, k] * rho[i, k] / mu[i, k])**0.5 *
                                sc[i, k]**(1.0 / 3.0) * cons13 /
                                (lamr**(5.0 / 2.0 + br / 2.0)))

                        pre = epsr * (qclr - qvs) / ab

                        # Only evaporate in out-of-cloud region
                        # and distribute across cldmax
                        pre = min(pre * (cldmax[i, k] - dum), 0.0)
                        pre = pre / cldmax[i, k]
                        am_evp_st[i, k] = max(cldmax[i, k] - dum, 0.0)

                    # Sublimation of snow
                    if qniic[i, k] >= qsmall:
                        qvi = svp_to_qsat(esi[i, k], p[i, k])
                        dqsidt = xxls * qvi / (rv * t[i, k]**2)
                        abi = 1.0 + dqsidt * xxls / cpp
                        epss = 2.0 * pi * n0s * rho[i, k] * Dv[i, k] * \
                               (f1s / (lams * lams) +
                                f2s * (asn[i, k] * rho[i, k] / mu[i, k])**0.5 *
                                sc[i, k]**(1.0 / 3.0) * cons14 /
                                (lams**(5.0 / 2.0 + bs / 2.0)))
                        prds = epss * (qclr - qvi) / abi

                        # Only sublimate in out-of-cloud region and distribute over cldmax
                        prds = min(prds * (cldmax[i, k] - dum), 0.0)
                        prds = prds / cldmax[i, k]
                        am_evp_st[i, k] = max(cldmax[i, k] - dum, 0.0)

                    # Make sure RH not pushed above 100% due to rain evaporation/snow sublimation
                    # Get updated RH at end of time step based on cloud water/ice condensation/evap
                    qtmp = q[i, k] - (cmei[i, k] + (pre + prds) * cldmax[i, k]) * deltat
                    ttmp = t[i, k] + ((pre * cldmax[i, k]) * xxlv +
                                      (cmei[i, k] + prds * cldmax[i, k]) * xxls) * deltat / cpp

                    # Limit range of temperatures!
                    ttmp = max(180.0, min(ttmp, 323.0))

                    esn = svp_water(ttmp)  # Use rhw to allow ice supersaturation
                    qsn = svp_to_qsat(esn, p[i, k])

                    # Modify precip evaporation rate if q > qsat
                    if qtmp > qsn:
                        if pre + prds < -1e-20:
                            dum1 = pre / (pre + prds)
                            # Recalculate q and t after cloud water cond but without precip evap
                            qtmp = q[i, k] - (cmei[i, k]) * deltat
                            ttmp = t[i, k] + (cmei[i, k] * xxls) * deltat / cpp
                            esn = svp_water(ttmp)  # Use rhw to allow ice supersaturation
                            qsn = svp_to_qsat(esn, p[i, k])
                            dum = (qtmp - qsn) / (1.0 + cons27 * qsn / (cpp * rv * ttmp**2))
                            dum = min(dum, 0.0)

                            # Modify rates if needed, divide by cldmax to get local (in-precip) value
                            pre = dum * dum1 / deltat / cldmax[i, k]

                            # Do separately using RHI for prds....
                            esn = svp_ice(ttmp)  # Use rhi to allow ice supersaturation
                            qsn = svp_to_qsat(esn, p[i, k])
                            dum = (qtmp - qsn) / (1.0 + cons28 * qsn / (cpp * rv * ttmp**2))
                            dum = min(dum, 0.0)

                            # Modify rates if needed, divide by cldmax to get local (in-precip) value
                            prds = dum * (1.0 - dum1) / deltat / cldmax[i, k]

                # Bergeron process - evaporation of droplets and deposition onto snow
                if qniic[i, k] >= qsmall and qcic >= qsmall and t[i, k] < tmelt:
                    qvi = svp_to_qsat(esi[i, k], p[i, k])
                    qvs = svp_to_qsat(esl[i, k], p[i, k])
                    dqsidt = xxls * qvi / (rv * t[i, k]**2)
                    abi = 1.0 + dqsidt * xxls / cpp
                    epss = 2.0 * pi * n0s * rho[i, k] * Dv[i, k] * \
                           (f1s / (lams * lams) +
                            f2s * (asn[i, k] * rho[i, k] / mu[i, k])**0.5 *
                            sc[i, k]**(1.0 / 3.0) * cons14 /
                            (lams**(5.0 / 2.0 + bs / 2.0)))
                    bergs = epss * (qvs - qvi) / abi
                else:
                    bergs = 0.0

                # Conservation to ensure no negative values of cloud water/precipitation
                # In case microphysical process rates are large

                # Make sure and use end-of-time step values for cloud water, ice, due
                # condensation/deposition

                # Note: for check on conservation, processes are multiplied by omsm
                # to prevent problems due to round off error

                # Include mixing timescale (mtime)

                qce = (qc[i, k] - berg[i, k] * deltat)
                nce = (nc[i, k] + npccn * deltat * mtime)
                qie = (qi[i, k] + (cmei[i, k] + berg[i, k]) * deltat)
                nie = (ni[i, k] + nnuccd * deltat * mtime)

                # Conservation of qc
                dum = (prc + pra + mnuccc + mnucct + msacwi +
                       psacws + bergs) * lcldm * deltat

                if dum > qce:
                    ratio = qce / deltat / lcldm / (prc + pra + mnuccc + mnucct + msacwi + psacws + bergs) * omsm

                    prc *= ratio
                    pra *= ratio
                    mnuccc *= ratio
                    mnucct *= ratio
                    msacwi *= ratio
                    psacws *= ratio
                    bergs *= ratio

                # Conservation of nc
                dum = (nprc1 + npra + nnuccc + nnucct +
                       npsacws - nsubc) * lcldm * deltat

                if dum > nce:
                    ratio = nce / deltat / ((nprc1 + npra + nnuccc + nnucct +
                                             npsacws - nsubc) * lcldm) * omsm

                    nprc1 *= ratio
                    npra *= ratio
                    nnuccc *= ratio
                    nnucct *= ratio
                    npsacws *= ratio
                    nsubc *= ratio

                # Conservation of qi
                if do_cldice:
                    frztmp = -mnuccc - mnucct - msacwi
                    if use_hetfrz_classnuc:
                        frztmp = -mnuccc - mnucct - mnudep - msacwi
                    dum = (frztmp * lcldm + (prci + prai) * icldm) * deltat

                    if dum > qie:
                        frztmp = mnuccc + mnucct + msacwi
                        if use_hetfrz_classnuc:
                            frztmp = mnuccc + mnucct + mnudep + msacwi
                        ratio = (qie / deltat + frztmp * lcldm) / ((prci + prai) * icldm) * omsm
                        prci *= ratio
                        prai *= ratio

                    # Conservation of ni
                    frztmp = -nnucct - nsacwi
                    if use_hetfrz_classnuc:
                        frztmp = -nnucct - nnuccc - nnudep - nsacwi
                    dum = (frztmp * lcldm + (nprci + nprai - nsubi) * icldm) * deltat

                    if dum > nie:
                        frztmp = nnucct + nsacwi
                        if use_hetfrz_classnuc:
                            frztmp = nnucct + nnuccc + nnudep + nsacwi
                        ratio = (nie / deltat + frztmp * lcldm) / \
                                ((nprci + nprai - nsubi) * icldm) * omsm
                        nprci *= ratio
                        nprai *= ratio
                        nsubi *= ratio

                # For precipitation conservation, use logic that vertical integral
                # of tendency from current level to top of model (i.e., qrtot) cannot be negative

                # Conservation of rain mixing rat
                if ((prc + pra) * lcldm + (-mnuccr + pre - pracs) *
                    cldmax[i, k]) * dz[i, k] * rho[i, k] + qrtot < 0.0:

                    if -pre + pracs + mnuccr >= qsmall:
                        ratio = (qrtot / (dz[i, k] * rho[i, k]) + (prc + pra) * lcldm) / \
                                ((-pre + pracs + mnuccr) * cldmax[i, k]) * omsm

                        pre *= ratio
                        pracs *= ratio
                        mnuccr *= ratio

                # Conservation of nr
                # For now neglect evaporation of nr
                nsubr = 0.0

                if (nprc * lcldm + (-nnuccr + nsubr - npracs + nragg) * cldmax[i, k]) * dz[i, k] * rho[i, k] + nrtot < 0.0:
                    if -nsubr - nragg + npracs + nnuccr >= qsmall:
                        ratio = (nrtot / (dz[i, k] * rho[i, k]) + nprc * lcldm) / \
                                ((-nsubr - nragg + npracs + nnuccr) * cldmax[i, k]) * omsm

                        nsubr *= ratio
                        npracs *= ratio
                        nnuccr *= ratio
                        nragg *= ratio

                # Conservation of snow mix ratio
                if ((bergs + psacws) * lcldm + (prai + prci) * icldm + (pracs +
                    mnuccr + prds) * cldmax[i, k]) * dz[i, k] * rho[i, k] + qstot < 0.0:

                    if -prds >= qsmall:
                        ratio = (qstot / (dz[i, k] * rho[i, k]) + (bergs + psacws) * lcldm + (prai + prci) * icldm +
                                 (pracs + mnuccr) * cldmax[i, k]) / (-prds * cldmax[i, k]) * omsm

                        prds *= ratio

                # Conservation of ns
                # Calculate loss of number due to sublimation
                # For now neglect sublimation of ns
                nsubs = 0.0

                if (nprci * icldm + (nnuccr + nsubs + nsagg) * cldmax[i, k]) * \
                   dz[i, k] * rho[i, k] + nstot < 0.0:

                    if -nsubs - nsagg >= qsmall:
                        ratio = (nstot / (dz[i, k] * rho[i, k]) + nprci * icldm +
                                 nnuccr * cldmax[i, k]) / ((-nsubs - nsagg) * cldmax[i, k]) * omsm

                        nsubs *= ratio
                        nsagg *= ratio

                # Get tendencies due to microphysical conversion processes
                # Note: tendencies are multiplied by appropriate cloud/precip
                # fraction to get grid-scale values
                # Note: cmei is already grid-average values

                qvlat[i, k] += -(pre + prds) * cldmax[i, k] - cmei[i, k]

                tlat[i, k] += ((pre * cldmax[i, k]) * xxlv +
                               (prds * cldmax[i, k] + cmei[i, k]) * xxls +
                               ((bergs + psacws + mnuccc + mnucct + msacwi) * lcldm + (mnuccr +
                                pracs) * cldmax[i, k] + berg[i, k]) * xlf)

                qctend[i, k] += (-pra - prc - mnuccc - mnucct - msacwi -
                                 psacws - bergs) * lcldm - berg[i, k]

                if do_cldice:
                    frztmp = mnuccc + mnucct + msacwi
                    if use_hetfrz_classnuc:
                        frztmp = mnuccc + mnucct + mnudep + msacwi
                    qitend[i, k] += frztmp * lcldm + \
                                    (-prci - prai) * icldm + cmei[i, k] + berg[i, k]

                qrtend[i, k] += (pra + prc) * lcldm + (pre - pracs -
                                                       mnuccr) * cldmax[i, k]

                qnitend[i, k] += (prai + prci) * icldm + (psacws + bergs) * lcldm + (prds +
                                                                                     pracs + mnuccr) * cldmax[i, k]

                # Add output for cmei (accumulate)
                cmeiout[i, k] += cmei[i, k]

                # Assign variables for trop_mozart, these are grid-average
                # Evaporation/sublimation is stored here as positive term
                evapsnow[i, k] += -prds * cldmax[i, k]
                nevapr[i, k] += -pre * cldmax[i, k]
                nevapr2[i, k] += -pre * cldmax[i, k]

                # Change to make sure prain is positive: do not remove snow from
                # prain used for wet deposition
                prain[i, k] += (pra + prc) * lcldm + (-pracs -
                                                      mnuccr) * cldmax[i, k]
                prodsnow[i, k] += (prai + prci) * icldm + (psacws + bergs) * lcldm + (
                    pracs + mnuccr) * cldmax[i, k]

                # Following are used to calculate 1st order conversion rate of cloud water
                # to rain and snow (1/s), for later use in aerosol wet removal routine
                # Previously, wetdepa used (prain/qc) for this, and the qc in wetdepa may be smaller than the qc
                # used to calculate pra, prc, ... in this routine
                # qcsinksum_rate1ord = sum over iterations{ rate of direct transfer of cloud water to rain & snow }
                #                      (no cloud ice or bergeron terms)
                # qcsum_rate1ord     = sum over iterations{ qc used in calculation of the transfer terms }

                qcsinksum_rate1ord[k] += (pra + prc + psacws) * lcldm
                qcsum_rate1ord[k] += qc[i, k]

                # Microphysics output, note this is grid-averaged
                prao[i, k] += pra * lcldm
                prco[i, k] += prc * lcldm
                mnuccco[i, k] += mnuccc * lcldm
                mnuccto[i, k] += mnucct * lcldm
                msacwio[i, k] += msacwi * lcldm
                psacwso[i, k] += psacws * lcldm
                bergso[i, k] += bergs * lcldm
                bergo[i, k] += berg[i, k]
                prcio[i, k] += prci * icldm
                praio[i, k] += prai * icldm
                mnuccro[i, k] += mnuccr * cldmax[i, k]
                pracso[i, k] += pracs * cldmax[i, k]

                # Multiply activation/nucleation by mtime to account for fast timescale
                nctend[i, k] += npccn * mtime + \
                                (-nnuccc - nnucct - npsacws + nsubc -
                                 npra - nprc1) * lcldm

                if do_cldice:
                    frztmp = nnucct + nsacwi
                    if use_hetfrz_classnuc:
                        frztmp = nnucct + nnuccc + nnudep + nsacwi
                    nitend[i, k] += nnuccd * mtime + \
                                    frztmp * lcldm + (nsubi - nprci - nprai) * icldm

                nstend[i, k] += (nsubs +
                                 nsagg + nnuccr) * cldmax[i, k] + nprci * icldm

                nrtend[i, k] += nprc * lcldm + (nsubr - npracs - nnuccr +
                                                nragg) * cldmax[i, k]

                # Make sure that nc and ni at advanced time step do not exceed
                # maximum (existing N + source terms*dt), which is possible due to
                # fast nucleation timescale

                if nctend[i, k] > 0.0 and nc[i, k] + nctend[i, k] * deltat > ncmax:
                    nctend[i, k] = max(0.0, (ncmax - nc[i, k]) / deltat)

                if do_cldice and nitend[i, k] > 0.0 and ni[i, k] + nitend[i, k] * deltat > nimax:
                    nitend[i, k] = max(0.0, (nimax - ni[i, k]) / deltat)

                # Get final values for precipitation q and N, based on
                # flux of precip from above, source/sink term, and terminal fallspeed
                # See eq. 15-16 in MG2008

                # Rain
                if qric[i, k] >= qsmall:
                    if k == top_lev - 1:
                        qric[i, k] = qrtend[i, k] * dz[i, k] / cldmax[i, k] / umr
                        nric[i, k] = nrtend[i, k] * dz[i, k] / cldmax[i, k] / unr
                    else:
                        qric[i, k] = (rho[i, k - 1] * umr * qric[i, k - 1] * cldmax[i, k - 1] +
                                      (rho[i, k] * dz[i, k] * qrtend[i, k])) / (umr * rho[i, k] * cldmax[i, k])
                        nric[i, k] = (rho[i, k - 1] * unr * nric[i, k - 1] * cldmax[i, k - 1] +
                                      (rho[i, k] * dz[i, k] * nrtend[i, k])) / (unr * rho[i, k] * cldmax[i, k])
                else:
                    qric[i, k] = 0.0
                    nric[i, k] = 0.0

                # Snow
                if qniic[i, k] >= qsmall:
                    if k == top_lev - 1:
                        qniic[i, k] = qnitend[i, k] * dz[i, k] / cldmax[i, k] / ums
                        nsic[i, k] = nstend[i, k] * dz[i, k] / cldmax[i, k] / uns
                    else:
                        qniic[i, k] = (rho[i, k - 1] * ums * qniic[i, k - 1] * cldmax[i, k - 1] +
                                       (rho[i, k] * dz[i, k] * qnitend[i, k])) / (ums * rho[i, k] * cldmax[i, k])
                        nsic[i, k] = (rho[i, k - 1] * uns * nsic[i, k - 1] * cldmax[i, k - 1] +
                                      (rho[i, k] * dz[i, k] * nstend[i, k])) / (uns * rho[i, k] * cldmax[i, k])
                else:
                    qniic[i, k] = 0.0
                    nsic[i, k] = 0.0

                # Calculate precipitation flux at surface
                # Divide by density of water to get units of m/s
                prect[i] += (qrtend[i, k] * dz[i, k] * rho[i, k] +
                             qnitend[i, k] * dz[i, k] * rho[i, k]) / rhow
                preci[i] += qnitend[i, k] * dz[i, k] * rho[i, k] / rhow

                # Convert rain rate from m/s to mm/hr
                rainrt[i, k] = qric[i, k] * rho[i, k] * umr / rhow * 3600.0 * 1000.0

                # Vertically-integrated precip source/sink terms (note: grid-averaged)
                qrtot = max(qrtot + qrtend[i, k] * dz[i, k] * rho[i, k], 0.0)
                qstot = max(qstot + qnitend[i, k] * dz[i, k] * rho[i, k], 0.0)
                nrtot = max(nrtot + nrtend[i, k] * dz[i, k] * rho[i, k], 0.0)
                nstot = max(nstot + nstend[i, k] * dz[i, k] * rho[i, k], 0.0)

                # Calculate melting and freezing of precip

                # Melt snow at +2 C
                if t[i, k] + tlat[i, k] / cpp * deltat > 275.15:
                    if qstot > 0.0:
                        # Make sure melting snow doesn't reduce temperature below threshold
                        dum = -xlf / cpp * qstot / (dz[i, k] * rho[i, k])
                        if t[i, k] + tlat[i, k] / cpp * deltat + dum < 275.15:
                            dum = (t[i, k] + tlat[i, k] / cpp * deltat - 275.15) * cpp / xlf
                            dum = dum / (xlf / cpp * qstot / (dz[i, k] * rho[i, k]))
                            dum = max(0.0, dum)
                            dum = min(1.0, dum)
                        else:
                            dum = 1.0

                        qric[i, k] += dum * qniic[i, k]
                        nric[i, k] += dum * nsic[i, k]
                        qniic[i, k] = (1.0 - dum) * qniic[i, k]
                        nsic[i, k] = (1.0 - dum) * nsic[i, k]
                        # Heating tendency
                        tmp = -xlf * dum * qstot / (dz[i, k] * rho[i, k])
                        meltsdt[i, k] += tmp

                        tlat[i, k] += tmp
                        qrtot += dum * qstot
                        nrtot += dum * nstot
                        qstot = (1.0 - dum) * qstot
                        nstot = (1.0 - dum) * nstot
                        preci[i] = (1.0 - dum) * preci[i]

                # Freeze all rain at -5C for Arctic
                if t[i, k] + tlat[i, k] / cpp * deltat < (tmelt - 5.0):
                    if qrtot > 0.0:
                        # Make sure freezing rain doesn't increase temperature above threshold
                        dum = xlf / cpp * qrtot / (dz[i, k] * rho[i, k])
                        if t[i, k] + tlat[i, k] / cpp * deltat + dum > (tmelt - 5.0):
                            dum = -(t[i, k] + tlat[i, k] / cpp * deltat - (tmelt




