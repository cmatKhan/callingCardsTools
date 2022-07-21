import scipy.stats as scistat

def compute_cumulative_poisson(exp_hops_region,bg_hops_region,total_exp_hops,total_bg_hops,pseudocounts):
    #usage
    #scistat.poisson.cdf(x,mu)
    #scales sample with more hops down to sample with less hops
    #tested 6/14/17
    if total_bg_hops >= total_exp_hops:
        auc = scistat.poisson.cdf((exp_hops_region+pseudocounts), 
                                  bg_hops_region * 
                                  (float(total_exp_hops) / float(total_bg_hops)) + 
                                  pseudocounts)
    else:
        auc = scistat.poisson.cdf(((exp_hops_region * (float(total_bg_hops)/float(total_exp_hops)) ) + pseudocounts),
                                  bg_hops_region + pseudocounts)
    
    return(1-auc)
