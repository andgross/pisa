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
