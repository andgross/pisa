#! /usr/bin/env python
## IMPORTS ##
import os,sys
import numpy as np
import logging
from argparse import ArgumentParser, RawTextHelpFormatter
#from utils.utils import set_verbosity,get_smoothed_map,get_osc_probLT_dict_hdf5
#from utils.json import from_json, to_json
#from flux.HondaFlux import get_flux_maps,HondaFlux

from icecube import dataio, icetray, dataclasses

# Until python2.6, default json is very slow.
try:
    import simplejson as json
except ImportError, e:
    import json

def hadronic_factor_function(hadronic_factor, leptonic_energy):
    import numpy as np
    E0 = 0.18791678;
    m  = 0.16267529;
    f0 = 0.30974123;

    return(pow(hadronic_factor,-m)*(hadronic_factor-1)+pow(leptonic_energy/E0,-m)*(1-f0));

def hadronicFactor(leptonic_energy):
    import numpy as np

    precision = 1.e-10;
    E0 = 0.18791678;
    m  = 0.16267529;
    f0 = 0.30974123;

    minlimit_hadronic_factor = 1-pow(np.exp(1./E0),-m)*(1-f0);
    min_hadronic_factor = 0.5
    max_hadronic_factor = 1.0

    min_function = hadronic_factor_function(min_hadronic_factor,leptonic_energy)
    max_function = hadronic_factor_function(max_hadronic_factor,leptonic_energy)
    if(max_function*min_function>0):
        return minlimit_hadronic_factor
    
    while(max_hadronic_factor-min_hadronic_factor > precision):
        cur_factor = (min_hadronic_factor+max_hadronic_factor)/2
        cur_function = hadronic_factor_function(cur_factor,leptonic_energy)
        if(cur_function*min_function > 0):
            min_hadronic_factor = cur_factor
            min_function = cur_function
        else: 
            max_hadronic_factor = cur_factor
            max_function = cur_function
            
  
    hadronic_factor = (min_hadronic_factor+max_hadronic_factor)/2;
    if (hadronic_factor < minlimit_hadronic_factor):
           return minlimit_hadronic_factor
    return(hadronic_factor)

def run_show_starge3():
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    czbins=np.linspace(-1,0,21)
    ebins=np.logspace(0,2,81)
    maps=get_oscillated_event_count_map_direct_from_MC_using_neutrinoflux_HDF( ebins, czbins,'/data/PINGU/multiNestProcessed/atmoweights/numu_v3cuts_newfiles_2*_*xx_stdproc.hdf5','/data/PINGU/multiNestProcessed/atmoweights/nue_v3cuts_newfiles_2*_*xx_stdproc.hdf5','/data/PINGU/multiNestProcessed/atmoweights/nutau_v3cuts_newfiles_2*_*xx_stdproc.hdf5',0.5,0.0024)

    e=maps['nue']
    ebar=maps['nuebar']
    mu=maps['numu']
    mubar=maps['numubar']
    tau=maps['nutau']
    taubar=maps['nutaubar']

    f1, axarr1=plt.subplots(3, 2,squeeze=False)

    axarr1[0][0].pcolormesh(mu['czbins'], mu['ebins'],mu['map'].transpose())
    axarr1[0][0].semilogy()
    axarr1[0][1].pcolormesh(mubar['czbins'], mubar['ebins'],mubar['map'].transpose())
    axarr1[0][1].semilogy()
    axarr1[1][0].pcolormesh(e['czbins'], e['ebins'],e['map'].transpose())
    axarr1[1][0].semilogy()
    axarr1[1][1].pcolormesh(ebar['czbins'], ebar['ebins'],ebar['map'].transpose())
    axarr1[1][1].semilogy()
    axarr1[2][0].pcolormesh(tau['czbins'], tau['ebins'],tau['map'].transpose())
    axarr1[2][0].semilogy()
    axarr1[2][1].pcolormesh(taubar['czbins'], taubar['ebins'],taubar['map'].transpose())
    axarr1[2][1].semilogy()


def run_compare_starge3():
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    czbins=np.linspace(-1,0,21)
    ebins=np.logspace(0,2,81)
    map1=get_oscillated_event_count_map_direct_from_MC_using_neutrinoflux_HDF( ebins, czbins,'/data/PINGU/multiNestProcessed/atmoweights/numu_v3cuts_newfiles_2*_*xx_stdproc.hdf5','/data/PINGU/multiNestProcessed/atmoweights/nue_v3cuts_newfiles_20*_*xx_stdproc.hdf5','/data/PINGU/multiNestProcessed/atmoweights/nutau_v3cuts_newfiles_218_*xx_stdproc.hdf5',0.39,0.0024)

    map2=get_oscillated_event_count_map_direct_from_MC_using_neutrinoflux_HDF( ebins, czbins,'/data/PINGU/multiNestProcessed/atmoweights/numu_v3cuts_newfiles_2*_*xx_stdproc.hdf5','/data/PINGU/multiNestProcessed/atmoweights/nue_v3cuts_newfiles_20*_*xx_stdproc.hdf5','/data/PINGU/multiNestProcessed/atmoweights/nutau_v3cuts_newfiles_218_*xx_stdproc.hdf5',0.39,-0.0024)
  
    e1=map1['nue']
    ebar1=map1['nuebar']
    mu1=map1['numu']
    mubar1=map1['numubar']
    tau1=map1['nutau']
    taubar1=map1['nutaubar']

    e2=map2['nue']
    ebar2=map2['nuebar']
    mu2=map2['numu']
    mubar2=map2['numubar']
    tau2=map2['nutau']
    taubar2=map2['nutaubar']


    f1, axarr1=plt.subplots(3, 2,squeeze=False)

    axarr1[0][0].pcolormesh(mu1['czbins'], mu1['ebins'],mu1['map'].transpose())
    axarr1[0][0].semilogy()
    axarr1[0][1].pcolormesh(mubar1['czbins'], mubar1['ebins'],mubar1['map'].transpose())
    axarr1[0][1].semilogy()
    axarr1[1][0].pcolormesh(mu2['czbins'], mu2['ebins'],mu2['map'].transpose())
    axarr1[1][0].semilogy()
    axarr1[1][1].pcolormesh(mubar2['czbins'], mubar2['ebins'],mubar2['map'].transpose())
    axarr1[1][1].semilogy()
    axarr1[1][0].pcolormesh(mu2['czbins'], mu2['ebins'],mu2['map'].transpose())
    axarr1[1][0].semilogy()
    axarr1[2][0].pcolormesh(mu2['czbins'], mu2['ebins'],(mu2['map'].transpose()+mubar2['map'].transpose()-mu1['map'].transpose()-mubar1['map'].transpose())/(np.sqrt(mu1['map'].transpose()+mubar1['map'].transpose())),vmin=-1., vmax=1.)
    axarr1[2][0].semilogy()
    axarr1[2][1].pcolormesh(e2['czbins'], e2['ebins'],(e2['map'].transpose()+ebar2['map'].transpose()-e1['map'].transpose()-ebar1['map'].transpose())/(np.sqrt(e1['map'].transpose()+ebar1['map'].transpose())),vmin=-1., vmax=1.)
    axarr1[2][1].semilogy()
    print (mu2['map'].transpose()+mubar2['map'].transpose()-mu1['map'].transpose()-mubar1['map'].transpose())/(np.sqrt(mu1['map'].transpose()+mubar1['map'].transpose()))
    print (e2['map'].transpose() +ebar2['map'].transpose() -e1['map'].transpose() -ebar1['map'].transpose())/ (np.sqrt(e1['map'].transpose() +ebar1['map'].transpose()))


def get_oscillated_event_count_map_direct_from_MC_using_neutrinoflux_HDF( ebins, czbins, numufiles='/data/PINGU/multiNestProcessed/atmoweights/numu_v3cuts_newfiles_201_*xx_stdproc.hdf5', nuefiles='/data/PINGU/multiNestProcessed/atmoweights/nue_v3cuts_newfiles_204_*xx_stdproc.hdf5', nutaufiles='/data/PINGU/multiNestProcessed/atmoweights/nutau_v3cuts_newfiles_218_*xx_stdproc.hdf5',sinsq_theta23=0.39, deltamsq_13=0.00216, **params):
    """
Primary module function. Produces a map of event counts in true coszen/ebins
for each of nue, nue_bar, numu, numu_bar, nutau, nutau_bar (6 maps total).
Inputs:
-- MCFiles for nue, numu, nutau
"""
    from numpy import cos
    from icecube import neutrinoflux
    import numpy as np
    from glob import glob
    import tables

    v_hadronicFactor=np.vectorize(hadronicFactor)

    mufilelist =  glob(numufiles) 
    mufilelist.sort()
    flux_numu=flux_nue=p_nue_numu=p_numu_numu=cz=cz_reco=energy=energy_reco=nt=NEvents=ow=np.array([])
    nfiles=0  
    for filename in mufilelist:
        nfiles+=100
        h5file=tables.openFile(filename)

        mc = h5file.root.I3MCWeightDict
        ow_tmp = mc.col('OneWeight')
        NEvents_tmp = mc.col('NEvents')
   
        neutrino=h5file.root.PrimaryNu
        nt_tmp=neutrino.col('type')
        energy_tmp=neutrino.col('energy')
        cz_tmp=cos(neutrino.col('zenith'))

        mn_reco=h5file.root.MultiNest_BestFitParticle
        cz_reco_tmp=cos(mn_reco.col('zenith'))
        energy_reco_tmp=mn_reco.col('energy')/v_hadronicFactor(mn_reco.col('energy'))+mn_reco.col('length')/4.5

        nm=h5file.root.numuflux
        numuflux_tmp=nm.col('value')
        ne=h5file.root.nueflux
        nueflux_tmp=ne.col('value')

        if (deltamsq_13>0):
            oscparams="/OscillationProbabilities_sth23_"+str(sinsq_theta23)+"_sth13_0.024_dm13_"+str(deltamsq_13)
        else:
            oscparams="/OscillationProbabilities_sth23_"+str(sinsq_theta23)+"_sth13_0.024_dm13_negative_"+str(-deltamsq_13)
        op=h5file.getNode(oscparams)
        p_numu_numu_tmp=op.col('PNuMuNuMu')
        p_nue_numu_tmp=op.col('PNuENuMu')
       

        ow=np.concatenate([ow,ow_tmp])
        NEvents=np.concatenate([NEvents,NEvents_tmp])
        
        nt=np.concatenate([nt,nt_tmp])
        energy=np.concatenate([energy,energy_tmp])
        cz=np.concatenate([cz,cz_tmp])
        p_numu_numu=np.concatenate([p_numu_numu, p_numu_numu_tmp])
        p_nue_numu=np.concatenate([p_nue_numu, p_nue_numu_tmp])
        
        energy_reco=np.concatenate([energy_reco,energy_reco_tmp])
        cz_reco=np.concatenate([cz_reco,cz_reco_tmp])

        flux_numu=np.concatenate([flux_numu, numuflux_tmp])
        flux_nue=np.concatenate([flux_nue, nueflux_tmp])

        print len(energy)
        
    print "now the weights..."
    print nt
    atmw_osc=86400.*365.*(p_numu_numu*flux_numu*ow*1./(NEvents/2.0)+p_nue_numu*flux_nue*ow*1./(NEvents/2.0))
    atmw_osc/=nfiles
    atmw_noosc=86400.*365.*flux_numu*ow*1./(NEvents/2.0)
    atmw_noosc/=nfiles
    print sum(atmw_osc)
    H_numu, czbins, ebins= np.histogram2d( cz_reco[nt==int(dataclasses.I3Particle.NuMu)], energy_reco[nt==int(dataclasses.I3Particle.NuMu)], bins=(czbins, ebins), weights=atmw_osc[nt==int(dataclasses.I3Particle.NuMu)])
    H_numu_noosc,czbins, ebins= np.histogram2d( cz_reco[nt==int(dataclasses.I3Particle.NuMu)], energy_reco[nt==int(dataclasses.I3Particle.NuMu)], bins=(czbins, ebins), weights=atmw_noosc[nt==int(dataclasses.I3Particle.NuMu)])
    H_numubar,czbins, ebins= np.histogram2d( cz_reco[nt==int(dataclasses.I3Particle.NuMuBar)], energy_reco[nt==int(dataclasses.I3Particle.NuMuBar)], bins=(czbins, ebins), weights=atmw_osc[nt==int(dataclasses.I3Particle.NuMuBar)])
    H_numubar_noosc,czbins, ebins= np.histogram2d( cz_reco[nt==int(dataclasses.I3Particle.NuMuBar)], energy_reco[nt==int(dataclasses.I3Particle.NuMuBar)], bins=(czbins, ebins), weights=atmw_noosc[nt==int(dataclasses.I3Particle.NuMuBar)])
 

    H2,czbins,czbins= np.histogram2d( cz, cz_reco, bins=(czbins,czbins), weights=atmw_osc)


    efilelist =  glob(nuefiles) 
    efilelist.sort()

    flux_numu=flux_nue=p_nue_nue=p_numu_nue=cz=cz_reco=energy=energy_reco=nt=NEvents=ow=np.array([])
    nfiles=0
  
    for filename in efilelist:
        nfiles+=100
        h5file=tables.openFile(filename)

        mc = h5file.root.I3MCWeightDict
        ow_tmp = mc.col('OneWeight')
        NEvents_tmp = mc.col('NEvents')
   
        neutrino=h5file.root.PrimaryNu
        nt_tmp=neutrino.col('type')
        energy_tmp=neutrino.col('energy')
        cz_tmp=cos(neutrino.col('zenith'))

        mn_reco=h5file.root.MultiNest_BestFitParticle
        cz_reco_tmp=cos(mn_reco.col('zenith'))
        energy_reco_tmp=mn_reco.col('energy')/v_hadronicFactor(mn_reco.col('energy'))+mn_reco.col('length')/4.5

        nm=h5file.root.numuflux
        numuflux_tmp=nm.col('value')
        ne=h5file.root.nueflux
        nueflux_tmp=ne.col('value')

        if (deltamsq_13>0):
            oscparams="/OscillationProbabilities_sth23_"+str(sinsq_theta23)+"_sth13_0.024_dm13_"+str(deltamsq_13)
        else:
            oscparams="/OscillationProbabilities_sth23_"+str(sinsq_theta23)+"_sth13_0.024_dm13_negative_"+str(-deltamsq_13)
        op=h5file.getNode(oscparams)
        p_numu_nue_tmp=op.col('PNuMuNuE')
        p_nue_nue_tmp=op.col('PNuENuE')
       

        ow=np.concatenate([ow,ow_tmp])
        NEvents=np.concatenate([NEvents,NEvents_tmp])
        
        nt=np.concatenate([nt,nt_tmp])
        energy=np.concatenate([energy,energy_tmp])
        cz=np.concatenate([cz,cz_tmp])
        p_numu_nue=np.concatenate([p_numu_nue, p_numu_nue_tmp])
        p_nue_nue=np.concatenate([p_nue_nue, p_nue_nue_tmp])
        
        energy_reco=np.concatenate([energy_reco,energy_reco_tmp])
        cz_reco=np.concatenate([cz_reco,cz_reco_tmp])

        flux_numu=np.concatenate([flux_numu, numuflux_tmp])
        flux_nue=np.concatenate([flux_nue, nueflux_tmp])

        print len(energy)
        
    print "now the weights..."
    
    atmw_osc=86400.*365.*(p_numu_nue*flux_numu*ow*1./(NEvents/2.0)+p_nue_nue*flux_nue*ow*1./(NEvents/2.0))
    atmw_osc/=nfiles
    atmw_noosc=86400.*365.*flux_nue*ow*1./(NEvents/2.0)
    atmw_noosc/=nfiles

    H_nue,czbins, ebins= np.histogram2d( cz_reco[nt==int(dataclasses.I3Particle.NuE)], energy_reco[nt==int(dataclasses.I3Particle.NuE)], bins=(czbins, ebins), weights=atmw_osc[nt==int(dataclasses.I3Particle.NuE)])
    H_nue_noosc,czbins, ebins= np.histogram2d( cz_reco[ nt == int(dataclasses.I3Particle.NuE)], energy_reco[nt==int(dataclasses.I3Particle.NuE)], bins=(czbins, ebins), weights=atmw_noosc[nt==int(dataclasses.I3Particle.NuE)])
    H_nuebar,czbins, ebins= np.histogram2d( cz_reco[nt==int(dataclasses.I3Particle.NuEBar)], energy_reco[nt==int(dataclasses.I3Particle.NuEBar)], bins=(czbins, ebins), weights=atmw_osc[nt==int(dataclasses.I3Particle.NuEBar)])
    H_nuebar_noosc,czbins, ebins= np.histogram2d( cz_reco[nt==int(dataclasses.I3Particle.NuEBar)], energy_reco[nt==int(dataclasses.I3Particle.NuEBar)], bins=(czbins, ebins), weights=atmw_noosc[nt==int(dataclasses.I3Particle.NuEBar)])


    taufilelist =  glob(nutaufiles) 
    taufilelist.sort()

    flux_numu=flux_nue=p_nue_nutau=p_numu_nutau=cz=cz_reco=energy=energy_reco=nt=NEvents=ow=np.array([])
    nfiles=0
  
    for filename in taufilelist:
        nfiles+=100
        h5file=tables.openFile(filename)

        mc = h5file.root.I3MCWeightDict
        ow_tmp = mc.col('OneWeight')
        NEvents_tmp = mc.col('NEvents')
   
        neutrino=h5file.root.MCNeutrino
        nt_tmp=neutrino.col('type')
        energy_tmp=neutrino.col('energy')
        cz_tmp=cos(neutrino.col('zenith'))

        mn_reco=h5file.root.MultiNest_BestFitParticle
        cz_reco_tmp=cos(mn_reco.col('zenith'))
        energy_reco_tmp=mn_reco.col('energy')/v_hadronicFactor(mn_reco.col('energy'))+mn_reco.col('length')/4.5

        nm=h5file.root.numuflux
        numuflux_tmp=nm.col('value')
        ne=h5file.root.nueflux
        nueflux_tmp=ne.col('value')

        if (deltamsq_13>0):
            oscparams="/OscillationProbabilities_sth23_"+str(sinsq_theta23)+"_sth13_0.024_dm13_"+str(deltamsq_13)
        else:
            oscparams="/OscillationProbabilities_sth23_"+str(sinsq_theta23)+"_sth13_0.024_dm13_negative_"+str(-deltamsq_13)
        op=h5file.getNode(oscparams)
        p_numu_nutau_tmp=op.col('PNuMuNuTau')
        p_nue_nutau_tmp=op.col('PNuENuTau')
       

        ow=np.concatenate([ow,ow_tmp])
        NEvents=np.concatenate([NEvents,NEvents_tmp])
        
        nt=np.concatenate([nt,nt_tmp])
        energy=np.concatenate([energy,energy_tmp])
        cz=np.concatenate([cz,cz_tmp])
        p_numu_nutau=np.concatenate([p_numu_nutau, p_numu_nutau_tmp])
        p_nue_nutau=np.concatenate([p_nue_nutau, p_nue_nutau_tmp])
        
        energy_reco=np.concatenate([energy_reco,energy_reco_tmp])
        cz_reco=np.concatenate([cz_reco,cz_reco_tmp])

        flux_numu=np.concatenate([flux_numu, numuflux_tmp])
        flux_nue=np.concatenate([flux_nue, nueflux_tmp])

        print len(energy)
        
    print "now the weights..."
    print nt
    atmw_osc=86400.*365.*(p_numu_nutau*flux_numu*ow*1./(NEvents/2.0)+p_nue_nutau*flux_nue*ow*1./(NEvents/2.0))
    atmw_osc/=nfiles

    H_nutau,czbins, ebins= np.histogram2d( cz_reco[nt==int(dataclasses.I3Particle.NuTau)], energy_reco[nt==int(dataclasses.I3Particle.NuTau)], bins=(czbins, ebins), weights=atmw_osc[nt==int(dataclasses.I3Particle.NuTau)])
    H_nutaubar,czbins, ebins= np.histogram2d( cz_reco[nt==int(dataclasses.I3Particle.NuTauBar)], energy_reco[nt==int(dataclasses.I3Particle.NuTauBar)], bins=(czbins, ebins), weights=atmw_osc[nt==int(dataclasses.I3Particle.NuTauBar)])


    #print H1
        
    maps = {}

    maps['numu']= { 'ebins' : ebins,
                        'czbins': czbins,
                        'map': H_numu}

    maps['numubar']= { 'ebins' : ebins,
                       'czbins': czbins,
                       'map': H_numubar}       

    
    maps['nue']= { 'ebins' : ebins,
                   'czbins': czbins,
                   'map': H_nue}

    maps['nuebar']= { 'ebins' : ebins,
                      'czbins': czbins,
                      'map': H_nuebar}   

    maps['nutau']= { 'ebins' : ebins,
                     'czbins': czbins,
                     'map': H_nutau}

    maps['nutaubar']= { 'ebins' : ebins,
                        'czbins': czbins,
                        'map': H_nutaubar}   

    return maps
        
