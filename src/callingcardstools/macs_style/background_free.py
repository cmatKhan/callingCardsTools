# local
from callingcardstools import macs_style as ms
# third party
import numpy as np
import pandas as pd
import scipy.stats as scistat

def call_peaks_and_annotate(expfile,  outputfile, TTAA_file, annotation_file, 
                               pvalue_cutoff      = 1e-4, 
                               peak_pvalue_cutoff = 1e-3, 
                               window_size        = 1000, 
                               lam_win_size       = 100000, 
                               step_size          = 500, 
                               pseudocounts       = 0.2):
    """This function is the wrapper that calls the various subfunctions
    TTAA_file = '/scratch/ref/rmlab/calling_card_ref/mouse/TTAA_mm10.txt'
    
    annotation_file = '/scratch/ref/rmlab/calling_card_ref/mouse/refGene.mm10.Sorted.bed'

    

    Args:
        expfile (_type_): _description_
        outputfile (_type_): _description_
        TTAA_file (_type_): _description_
        annotation_file (_type_): _description_
        pvalue_cutoff (_type_, optional): _description_. Defaults to 1e-4.
        peak_pvalue_cutoff (_type_, optional): _description_. Defaults to 1e-3.
        window_size (int, optional): _description_. Defaults to 1000.
        lam_win_size (int, optional): _description_. Defaults to 100000.
        step_size (int, optional): _description_. Defaults to 500.
        pseudocounts (float, optional): _description_. Defaults to 0.2.
    """

    expframe = pd.read_csv(expfile,delimiter="\t",header=None)
    expframe.columns = ["Chr","Start","End","Reads","Strand","Barcode"]

    TTAAframe = pd.read_csv(TTAA_file,delimiter="\t",header=None)
    TTAAframe.columns = ["Chr","Start","End"]

    peaks_frame = find_peaks_bf(expframe, 
                                TTAAframe, 
                                peak_pvalue_cutoff, 
                                window_size,
                                lam_win_size,
                                step_size,
                                pseudocounts)

    peaks_frame = ms.annotate.annotate_peaks_frame(peaks_frame,annotation_file)
    # filter rows based on pvalue_cutoff
    peaks_frame = peaks_frame[peaks_frame["poisson_pvalue"] <= pvalue_cutoff]
    # sort by pvalue
    peaks_frame = peaks_frame.sort_values(["poisson_pvalue"])
    # write to csv
    peaks_frame.to_csv(outputfile,sep="\t",index=False)
    peaksbed_frame = ms.annotate.make_peaksbed(peaks_frame)
    peaksbed_frame.columns = ['Chr','Start','End','Col4']

    peaksbed_frame_with_description = pd.merge(peaksbed_frame,
                                               peaks_frame,
                                               how = "left",
                                               on = ['Chr','Start','End'])

    peaksbed_frame_with_description['tph_experiment'] = \
        'tph_experiment=' + \
            peaksbed_frame_with_description['tph_experiment'].astype(str)
    peaksbed_frame_with_description['poisson_pvalue'] = \
        'poisson_pvalue=' + \
            peaksbed_frame_with_description['poisson_pvalue'].astype(str)
    peaksbed_frame_with_description['lambda'] = \
        'lambda=' + \
            peaksbed_frame_with_description['lambda'].astype(str)
    peaksbed_frame_with_description['lambda_type'] = \
        'lambda_type=' + \
            peaksbed_frame_with_description['lambda_type'].astype(str)

    peaksbed_frame_with_description['Description'] = \
        list(peaksbed_frame_with_description.loc[:,['tph_experiment', 
                                                    'poisson_pvalue',
                                                    'lambda',
                                                    'lambda_type']].values)

    final_peaksbed_frame = \
        peaksbed_frame_with_description.loc[:,['Chr', 
                                               'Start',
                                               'End',
                                               'Col4',
                                               'Description']]

    peaksbed_frame = final_peaksbed_frame

    pattern = "^(\S+)\."
    group = re.search(pattern,outputfile)
    bedfilename = group.groups(0)[0]+".tmp.peaks.bed"
    peaksbed_frame.to_csv(bedfilename,sep="\t",index=False,header=False)

    df = pd.read_csv(bedfilename,sep = "\t",header = None)
    df.columns = ['Chr','Start','End','Col4','Description']
    df['Description'] = df['Description'].str.replace('\n','')

    bedfilename = group.groups(0)[0]+".peaks.bed"
    df.to_csv(bedfilename,sep = "\t",header = None,index = False)

def find_peaks_bf(experiment_frame, TTAA_frame, pvalue_cutoff = 1e-3, 
                  window_size = 1000, lam_win_size=100000, step_size = 500, 
                  pseudocounts = 0.2):
    """This function is the same as find_peaks but it calls peaks without using
    a background distribution.  To do so, it scans through the genome
    with a window of window_size and step size of step_size and looks for regions that
    have signficantly more experiment hops in the window than expected from the number of
    experimenta hops in a window of size lam_win_size.
    It merges consecutively enriched windows and computes the center of the peak.  Next 
    it computes lambda, the number of insertions per TTAA expected from the experimental
    distribution where the expected number of hops are
    estimated from the number of hops in a window of size lam_win_size around the peak center.  It then computes
    a p-value based on the expected number of hops = lamda * number of TTAAs in peak * number of hops
    in peak.  Finally, it returns a frame that has Chr,Start,End,Center,experiment_hops,
    fraction_experiment lambda_type, lambda,poisson_pvalue

    :param experiment_frame:
    :param TTAA_frame:
    :param pvalue_cutoff:
    :param window_size:
    :param lam_win_size:
    :param step-size:
    :param pseudocounts:

    :returns:

    """
    peaks_frame_columns = ["Chr", "Start", "End", "Center", 
                           "experiment_hops", "fraction_experiment", 
                           "tph_experiment", "lambda_type", "lambda", 
                           "poisson_pvalue"]

    peaks_frame = pd.DataFrame(columns = peaks_frame_columns)

    
    experiment_gnashy_dict = {}
    experiment_dict_of_trees = {}
    total_experiment_hops = len(experiment_frame)
    TTAA_frame_gbChr_dict = {} 
    TTAA_dict_of_trees = {}
    list_of_l_names = lam_win_size

    #group by chromosome and populate interval tree with TTAA positions
    for name,group in TTAA_frame.groupby('Chr'):
         
        TTAA_frame_gbChr_dict[name] = group
        TTAA_frame_gbChr_dict[name].index = TTAA_frame_gbChr_dict[name]["Start"]
        #initialize tree
        TTAA_dict_of_trees[name] = Intersecter() 
        #populate tree with position as interval
        for idx, row in TTAA_frame_gbChr_dict[name].iterrows():    
            TTAA_dict_of_trees[name].add_interval(Interval(int(idx),int(idx+3)))

    #group by chromosome and populate interval tree with positions of experiment hops
    for name,group in experiment_frame.groupby('Chr'):
        experiment_gnashy_dict[name] = group
        experiment_gnashy_dict[name].index = experiment_gnashy_dict[name]["Start"]
        #initialize tree
        experiment_dict_of_trees[name] = Intersecter() 
        #populate tree with position as interval
        for idx, row in experiment_gnashy_dict[name].iterrows():    
            experiment_dict_of_trees[name].add_interval(Interval(int(idx),int(idx)+3)) 

    #these will eventually be the columns in the peaks frame that will be returned.
    chr_list = []
    start_list = []
    end_list = []
    center_list = []
    num_exp_hops_list = []
    frac_exp_list = []
    tph_exp_list = []
    lambda_type_list =[]
    lambda_list = []
    pvalue_list = []
    l = []
    
    #group experiment gnashyfile by chomosome
    for name,group in experiment_frame.groupby('Chr'):
        
        max_pos = max(group["End"])
        sig_start = 0
        sig_end = 0
        sig_flag = 0
        
        window_range = range(lam_win_size/2,max_pos + 
            window_size+(lam_win_size/2+1),
            step_size)

        for window_start in window_range:
            overlap = experiment_dict_of_trees[name]\
                .find(window_start, window_start + window_size - 1)
            num_exp_hops = len(overlap)
            overlap_lam_win_size = experiment_dict_of_trees[name]\
                .find(window_start - (lam_win_size/2-1), window_start + 
                                                         window_size +
                                                         (lam_win_size/2) - 1)

            num_lam_win_size_hops = len(overlap_lam_win_size)
            num_TTAAs_window = len(TTAA_dict_of_trees[name]\
                .find(window_start,window_start+window_size - 1))
            num_TTAAs_lam_win_size = len(TTAA_dict_of_trees[name]\
                .find(window_start - (lam_win_size/2-1), window_start + 
                                                         window_size + 
                                                         lam_win_size/2 - 1))
            
            # expected number of hops per TTAA
            lambda_lam_win_size = \
                (num_lam_win_size_hops/(max(num_TTAAs_lam_win_size,1))) 

            #is this window significant?
            pvalue = 1 - scistat.poisson.cdf(
                (num_exp_hops + pseudocounts),
                lambda_lam_win_size * max(num_TTAAs_window,1) + pseudocounts)

            if pvalue < pvalue_cutoff:
                #was last window significant?
                if sig_flag:
                    #if so, extend end of windows
                    sig_end = window_start+window_size-1
                else:
                    #otherwise, define new start and end and set flag
                    sig_start = window_start
                    sig_end = window_start+window_size-1
                    sig_flag = 1

            else:
                #current window not significant.  Was last window significant?
                if sig_flag:
                    #add full sig window to the frame of peaks 
                    
                    #add chr, peak start, peak end
                    chr_list.append(name) #add chr to frame
                    start_list.append(sig_start) #add peak start to frame
                    end_list.append(sig_end) #add peak end to frame
                
                    #compute peak center and add to frame
                    overlap = experiment_dict_of_trees[name]\
                        .find(sig_start,sig_end)
                    exp_hop_pos_list = [x.start for x in overlap]
                    peak_center = np.median(exp_hop_pos_list)
                    center_list.append(peak_center) #add peak center to frame

                    #add number of experiment hops in peak to frame
                    num_exp_hops = len(overlap)
                    num_exp_hops_list.append(num_exp_hops)

                    #add fraction of experiment hops in peak to frame
                    frac_exp_list\
                        .append(float(num_exp_hops)/total_experiment_hops)
                    tph_exp_list\
                        .append(float(num_exp_hops) * 
                                100000/total_experiment_hops)

                    num_TTAAs_peak = \
                        len(TTAA_dict_of_trees[name].find(sig_start,sig_end))

                    #compute lambda in lam_win_size
                    num_exp_hops_lam_win_size = \
                        len(experiment_dict_of_trees[name]\
                            .find(peak_center - (lam_win_size/2-1), 
                                  peak_center+(lam_win_size/2)))

                    num_TTAAs_lam_win_size = \
                        len(TTAA_dict_of_trees[name]\
                            .find(peak_center - (lam_win_size/2-1), 
                                  peak_center+(lam_win_size/2)))
                                  
                    lambda_lam_win_size = \
                        float(num_exp_hops_lam_win_size)/\
                            (max(num_TTAAs_lam_win_size,1))


                    lambda_f = lambda_lam_win_size


                    #record type of lambda used
                    lambda_type_list.append(list_of_l_names)
                    #record lambda
                    lambda_list.append(lambda_f)
                    #compute pvalue and record it

                    pvalue = 1-scistat.poisson.cdf(
                        (num_exp_hops + pseudocounts), 
                        lambda_f * max(num_TTAAs_peak,1) + pseudocounts)

                    pvalue_list.append(pvalue)              
                    #number of hops that are a user-defined distance from peak center
                    sig_flag = 0
        
                                        
                    lambdatype = lam_win_size
                    l = [pvalue,float(num_exp_hops) * 100000/\
                        total_experiment_hops,lambda_f, lambdatype]

                #else do nothing.
                
                
    #make frame from all of the lists
    peaks_frame["Chr"] = chr_list
    peaks_frame["Start"] = start_list
    peaks_frame["End"] = end_list
    peaks_frame["Center"] = center_list
    peaks_frame["experiment_hops"] = num_exp_hops_list 
    peaks_frame["fraction_experiment"] = frac_exp_list 
    peaks_frame["tph_experiment"]= tph_exp_list
    peaks_frame["lambda_type"] = lambda_type_list
    peaks_frame["lambda"] = lambda_list
    peaks_frame["poisson_pvalue"] = pvalue_list
    return peaks_frame
