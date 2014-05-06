#! /usr/bin/env python
## IMPORTS ##
import os,sys
import numpy as np
import logging
from argparse import ArgumentParser, RawTextHelpFormatter
from utils.utils import set_verbosity
from utils.json import from_json, to_json
import tables

# Until python2.6, default json is very slow.
try:
    import simplejson as json
except ImportError, e:
    import json

class InfoFromMc(object):
    def __init__(self):
        self.numu_flux=np.array([])
        self.nue_flux=np.array([])
        self.p_from_nue=np.array([])
        self.p_from_numu=np.array([])
        self.cz_true=np.array([])
        self.cz_reco=np.array([])
        self.energy_true=np.array([])
        self.energy_reco=np.array([])
        self.neutrino_type=np.array([])
        self.n_events=np.array([])
        self.ow=np.array([])


def run_show_starge3():
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    czbins=np.linspace(-1,0,21)
    ebins=np.logspace(0,2,81)
    maps=get_oscillated_event_count_map_direct_from_MC_using_neutrinoflux_HDF( ebins, czbins,'/data/PINGU/multiNestProcessed/atmoweights/numu_v3cuts_newfiles_2*_*xx_stdproc.hdf5','/data/PINGU/multiNestProcessed/atmoweights/nue_v3cuts_newfiles_2*_*xx_stdproc.hdf5','/data/PINGU/multiNestProcessed/atmoweights/nutau_v3cuts_newfiles_2*_*xx_stdproc.hdf5',0.5,0.0024)

    print maps
    e=maps['nue']
    ebar=maps['nuebar']
    mu=maps['numu']
    mubar=maps['numubar']
    tau=maps['nutau']
    taubar=maps['nutaubar']

    f1, axarr1=plt.subplots(3, 2,squeeze=False)
    print mu
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
    

def append_hdf_file( mcdata, filename, flavor,sinsq_theta23, deltamsq_13, primary="PrimaryNu"):
    import HadronicFactor
    v_hadronicFactor=np.vectorize(HadronicFactor.hadronicFactor)
    print filename
    h5file=tables.openFile(filename)
    numu_flux = h5file.root.numuflux.col('value')
    nue_flux  = h5file.root.nueflux.col('value')

    if (deltamsq_13>0):
        oscparams="/OscillationProbabilities_sth23_"+str(sinsq_theta23)+"_sth13_0.024_dm13_"+str(deltamsq_13)
    else:
        oscparams="/OscillationProbabilities_sth23_"+str(sinsq_theta23)+"_sth13_0.024_dm13_negative_"+str(-deltamsq_13)
    try:
        op=h5file.getNode(oscparams)
    except:
        logging.error("Unable to find oscillation parameters \'%s\'"%oscparams)
        sys.exit(1)


    if flavor=="NuE":
        p_from_nue =op.col('PNuENuE')
        p_from_numu=op.col('PNuMuNuE')
    elif flavor=="NuMu":
        p_from_nue =op.col('PNuENuMu')
        p_from_numu=op.col('PNuMuNuMu')
    elif flavor=="NuTau":
        p_from_nue =op.col('PNuENuTau')
        p_from_numu=op.col('PNuMuNuTau')
    else:
        logging.error("Unsupported neutrino flavor \'%s\'"%flavor)
        sys.exit(1)

    if primary=="PrimaryNu":
        cz_true=np.cos(h5file.root.PrimaryNu.col('zenith'))
        energy_true=h5file.root.PrimaryNu.col('energy')
    else:
        cz_true=np.cos(h5file.root.MCNeutrino.col('zenith'))
        energy_true=h5file.root.MCNeutrino.col('energy')

    cz_reco=np.cos(h5file.root.MultiNest_BestFitParticle.col('zenith'))

    energy_reco=h5file.root.MultiNest_BestFitParticle.col('energy')/v_hadronicFactor(h5file.root.MultiNest_BestFitParticle.col('energy'))+h5file.root.MultiNest_BestFitParticle.col('length')/4.5

    neutrino_type=h5file.root.PrimaryNu.col('type')
    n_events=h5file.root.I3MCWeightDict.col('NEvents')
    ow=h5file.root.I3MCWeightDict.col('OneWeight')

    mcdata.numu_flux  =np.concatenate([mcdata.numu_flux,numu_flux])
    mcdata.nue_flux   =np.concatenate([mcdata.nue_flux,nue_flux])
    mcdata.p_from_nue =np.concatenate([mcdata.p_from_nue,p_from_nue])
    mcdata.p_from_numu=np.concatenate([mcdata.p_from_numu,p_from_numu])
    mcdata.cz_true    =np.concatenate([mcdata.cz_true,cz_true])
    mcdata.cz_reco    =np.concatenate([mcdata.cz_reco,cz_reco])
    mcdata.energy_true=np.concatenate([mcdata.energy_true,energy_true])
    mcdata.energy_reco=np.concatenate([mcdata.energy_reco,energy_reco])
    mcdata.neutrino_type=np.concatenate([mcdata.neutrino_type,neutrino_type])
    mcdata.n_events=np.concatenate([mcdata.n_events,n_events])
    mcdata.ow=np.concatenate([mcdata.ow,ow])
    h5file.close()
    return mcdata



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

    numuevents=InfoFromMc()
    nfiles=0  
    for filename in mufilelist:
        numuevents=append_hdf_file( numuevents, filename, "NuMu",sinsq_theta23, deltamsq_13)
        logging.info("Reading in file %s", filename)
        nfiles+=100
        
        
    print len(numuevents.p_from_numu), len(numuevents.numu_flux), len(numuevents.ow), len(numuevents.n_events)
    atmw_osc=86400.*365.*(numuevents.p_from_numu*numuevents.numu_flux*numuevents.ow*1./(numuevents.n_events/2.0)+numuevents.p_from_nue*numuevents.nue_flux*numuevents.ow*1./(numuevents.n_events/2.0))
    atmw_osc/=nfiles

    print numuevents.cz_reco, numuevents.neutrino_type, numuevents.energy_reco, atmw_osc
    H_numu, czbins, ebins= np.histogram2d( numuevents.cz_reco[numuevents.neutrino_type==int(dataclasses.I3Particle.NuMu)], numuevents.energy_reco[numuevents.neutrino_type==int(dataclasses.I3Particle.NuMu)], bins=(czbins, ebins), weights=atmw_osc[numuevents.neutrino_type==int(dataclasses.I3Particle.NuMu)])
    H_numubar, czbins, ebins= np.histogram2d( numuevents.cz_reco[numuevents.neutrino_type==int(dataclasses.I3Particle.NuMuBar)], numuevents.energy_reco[numuevents.neutrino_type==int(dataclasses.I3Particle.NuMuBar)], bins=(czbins, ebins), weights=atmw_osc[numuevents.neutrino_type==int(dataclasses.I3Particle.NuMuBar)])

    print H_numu

    efilelist =  glob(nuefiles) 
    efilelist.sort()
    nueevents=InfoFromMc()
    nfiles=0
  
    for filename in efilelist:
        nueevents=append_hdf_file( nueevents, filename, "NuE",sinsq_theta23, deltamsq_13)
        logging.info("Reading in file %s", filename)
        nfiles+=100


    atmw_osc=86400.*365.*(nueevents.p_from_numu*nueevents.numu_flux*nueevents.ow*1./(nueevents.n_events/2.0)+nueevents.p_from_nue*nueevents.nue_flux*nueevents.ow*1./(nueevents.n_events/2.0))
    atmw_osc/=nfiles

    H_nue, czbins, ebins= np.histogram2d( nueevents.cz_reco[nueevents.neutrino_type==int(dataclasses.I3Particle.NuE)], nueevents.energy_reco[nueevents.neutrino_type==int(dataclasses.I3Particle.NuE)], bins=(czbins, ebins), weights=atmw_osc[nueevents.neutrino_type==int(dataclasses.I3Particle.NuE)])
    H_nuebar, czbins, ebins= np.histogram2d( nueevents.cz_reco[nueevents.neutrino_type==int(dataclasses.I3Particle.NuEBar)], nueevents.energy_reco[nueevents.neutrino_type==int(dataclasses.I3Particle.NuEBar)], bins=(czbins, ebins), weights=atmw_osc[nueevents.neutrino_type==int(dataclasses.I3Particle.NuEBar)])

    taufilelist =  glob(nutaufiles) 
    taufilelist.sort()

    flux_numu=flux_nue=p_nue_nutau=p_numu_nutau=cz=cz_reco=energy=energy_reco=nt=NEvents=ow=np.array([])
    nfiles=0
    nutauevents=InfoFromMc()
    for filename in taufilelist:
        nutauevents=append_hdf_file( nutauevents, filename, "NuTau",sinsq_theta23, deltamsq_13,"MCNeutrino")
        logging.info("Reading in file %s", filename)
        nfiles+=100


                
    atmw_osc=86400.*365.*(nutauevents.p_from_numu*nutauevents.numu_flux*nutauevents.ow*1./(nutauevents.n_events/2.0)+nutauevents.p_from_nue*nutauevents.nue_flux*nutauevents.ow*1./(nutauevents.n_events/2.0))
    atmw_osc/=nfiles

    H_nutau, czbins, ebins= np.histogram2d( nutauevents.cz_reco[nutauevents.neutrino_type==int(dataclasses.I3Particle.NuTau)], nutauevents.energy_reco[nutauevents.neutrino_type==int(dataclasses.I3Particle.NuTau)], bins=(czbins, ebins), weights=atmw_osc[nutauevents.neutrino_type==int(dataclasses.I3Particle.NuTau)])
    H_nutaubar, czbins, ebins= np.histogram2d( nutauevents.cz_reco[nutauevents.neutrino_type==int(dataclasses.I3Particle.NuTauBar)], nutauevents.energy_reco[nutauevents.neutrino_type==int(dataclasses.I3Particle.NuTauBar)], bins=(czbins, ebins), weights=atmw_osc[nutauevents.neutrino_type==int(dataclasses.I3Particle.NuTauBar)])

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

