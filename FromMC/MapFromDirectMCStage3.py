#! /usr/bin/env python
## IMPORTS ##
import os,sys
import numpy as np
import logging
from argparse import ArgumentParser, RawTextHelpFormatter
from utils.utils import set_verbosity
from utils.json import from_json, to_json


# Until python2.6, default json is very slow.
try:
    import simplejson as json
except ImportError, e:
    import json

def get_oscillated_event_count_map_direct_from_MC_using_neutrinoflux_HDF( ebins, czbins, numufiles, nuefiles, nutaufiles,sinsq_theta23=0.39, deltamsq_13=0.00216, **params):
    """
Primary module function. Produces a map of event counts in true coszen/ebins
for each of nue, nue_bar, numu, numu_bar, nutau, nutau_bar (6 maps total).
Inputs:
-- MCFiles for nue, numu, nutau, oscillation parameters
"""
    from numpy import cos
    from icecube import neutrinoflux, dataclasses
    import numpy as np
    from glob import glob
    import tables
    import HadronicFactor
    v_hadronicFactor=np.vectorize(HadronicFactor.hadronicFactor)

    mufilelist =  glob(numufiles) 
    mufilelist.sort()
    flux_numu=flux_nue=p_nue_numu=p_numu_numu=cz=cz_reco=energy=energy_reco=nt=NEvents=ow=np.array([])
    nfiles=0  
    for filename in mufilelist:
        logging.info("Reading in file %s", filename)
        nfiles+=100
        h5file=tables.openFile(filename)

        ow_tmp = h5file.root.I3MCWeightDict.col('OneWeight')
        NEvents_tmp = h5file.root.I3MCWeightDict.col('NEvents')
   
        nt_tmp=h5file.root.PrimaryNu.col('type')
        energy_tmp=h5file.root.PrimaryNu.col('energy')
        cz_tmp=cos(h5file.root.PrimaryNu.col('zenith'))

        cz_reco_tmp=cos(h5file.root.MultiNest_BestFitParticle.col('zenith'))
        energy_reco_tmp=h5file.root.MultiNest_BestFitParticle.col('energy')/v_hadronicFactor(h5file.root.MultiNest_BestFitParticle.col('energy'))+h5file.root.MultiNest_BestFitParticle.col('length')/4.5

        numuflux_tmp=h5file.root.numuflux.col('value')
        nueflux_tmp=h5file.root.nueflux.col('value')

        if (deltamsq_13>0):
            oscparams="/OscillationProbabilities_sth23_"+str(sinsq_theta23)+"_sth13_0.024_dm13_"+str(deltamsq_13)
        else:
            oscparams="/OscillationProbabilities_sth23_"+str(sinsq_theta23)+"_sth13_0.024_dm13_negative_"+str(-deltamsq_13)
        try:
            op=h5file.getNode(oscparams)
        except:
            logging.error("Unable to find oscillation parameters \'%s\'"%oscparams)
            sys.exit(1)

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

        
    atmw_osc=86400.*365.*(p_numu_numu*flux_numu*ow*1./(NEvents/2.0)+p_nue_numu*flux_nue*ow*1./(NEvents/2.0))
    atmw_osc/=nfiles

    H_numu, czbins, ebins= np.histogram2d( cz_reco[nt==int(dataclasses.I3Particle.NuMu)], energy_reco[nt==int(dataclasses.I3Particle.NuMu)], bins=(czbins, ebins), weights=atmw_osc[nt==int(dataclasses.I3Particle.NuMu)])
    H_numubar,czbins, ebins= np.histogram2d( cz_reco[nt==int(dataclasses.I3Particle.NuMuBar)], energy_reco[nt==int(dataclasses.I3Particle.NuMuBar)], bins=(czbins, ebins), weights=atmw_osc[nt==int(dataclasses.I3Particle.NuMuBar)])


    efilelist =  glob(nuefiles) 
    efilelist.sort()

    flux_numu=flux_nue=p_nue_nue=p_numu_nue=cz=cz_reco=energy=energy_reco=nt=NEvents=ow=np.array([])
    nfiles=0
  
    for filename in efilelist:
        logging.info("Reading in file %s", filename)
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
        try:
            op=h5file.getNode(oscparams)
        except:
            logging.error("Unable to find oscillation parameters \'%s\'"%oscparams)
            sys.exit(1)
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

    
    atmw_osc=86400.*365.*(p_numu_nue*flux_numu*ow*1./(NEvents/2.0)+p_nue_nue*flux_nue*ow*1./(NEvents/2.0))
    atmw_osc/=nfiles

    H_nue,czbins, ebins= np.histogram2d( cz_reco[nt==int(dataclasses.I3Particle.NuE)], energy_reco[nt==int(dataclasses.I3Particle.NuE)], bins=(czbins, ebins), weights=atmw_osc[nt==int(dataclasses.I3Particle.NuE)])
    H_nuebar,czbins, ebins= np.histogram2d( cz_reco[nt==int(dataclasses.I3Particle.NuEBar)], energy_reco[nt==int(dataclasses.I3Particle.NuEBar)], bins=(czbins, ebins), weights=atmw_osc[nt==int(dataclasses.I3Particle.NuEBar)])

    taufilelist =  glob(nutaufiles) 
    taufilelist.sort()

    flux_numu=flux_nue=p_nue_nutau=p_numu_nutau=cz=cz_reco=energy=energy_reco=nt=NEvents=ow=np.array([])
    nfiles=0
  
    for filename in taufilelist:
        logging.info("Reading in file %s", filename)
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
        try:
            op=h5file.getNode(oscparams)
        except:
            logging.error("Unable to find oscillation parameters \'%s\'"%oscparams)
            sys.exit(1)
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

                
    atmw_osc=86400.*365.*(p_numu_nutau*flux_numu*ow*1./(NEvents/2.0)+p_nue_nutau*flux_nue*ow*1./(NEvents/2.0))
    atmw_osc/=nfiles

    H_nutau,czbins, ebins= np.histogram2d( cz_reco[nt==int(dataclasses.I3Particle.NuTau)], energy_reco[nt==int(dataclasses.I3Particle.NuTau)], bins=(czbins, ebins), weights=atmw_osc[nt==int(dataclasses.I3Particle.NuTau)])
    H_nutaubar,czbins, ebins= np.histogram2d( cz_reco[nt==int(dataclasses.I3Particle.NuTauBar)], energy_reco[nt==int(dataclasses.I3Particle.NuTauBar)], bins=(czbins, ebins), weights=atmw_osc[nt==int(dataclasses.I3Particle.NuTauBar)])

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
        
if __name__ == '__main__':    
#def main():
    import numpy as np
    #from utils.json import from_json, to_json
    
    parser = ArgumentParser(description='Take hdf MC files'
                            'as input and write out a set of oscillated event counts. ',
                            formatter_class=RawTextHelpFormatter)
    parser.add_argument('-o', '--outfile', dest='outfile', metavar='FILE', type=str,
                        action='store',default="event_rate.json",
                        help='file to store the output')
    parser.add_argument('-i', '--indir', dest='indir', metavar='FILE', type=str,
                        action='store',default="/data/PINGU/multiNestProcessed/atmoweights/",
                        help='file to store the output')
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help='set verbosity level')
    args = parser.parse_args()

    #Set verbosity level
    set_verbosity(args.verbose)
    
    #default bins
    czbins=np.linspace(-1,0,21)
    ebins=np.logspace(0,2,81)

    #get a map
    event_rate_maps=get_oscillated_event_count_map_direct_from_MC_using_neutrinoflux_HDF( ebins, czbins, args.indir+'/numu_v3cuts_newfiles_*xx_stdproc.hdf5', args.indir+'/nue_v3cuts_newfiles_*_*xx_stdproc.hdf5', args.indir+'/nutau_v3cuts_newfiles_*_*xx_stdproc.hdf5',0.5,0.0024)
    
    #write the map
    to_json(event_rate_maps,args.outfile)
