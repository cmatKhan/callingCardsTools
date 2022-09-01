"""This module annotates peaks frames.
""" 

import pandas as pd
import pybedtools.bedtool as bt 
import random

def peaks_frame_to_bed(peaks_frame, bedfilename):
    """bladdity blah

    Args:
        peaks_frame (DataFrame): _description_
        bedfilename (Str): _description_
    """
    bed_frame = peaks_frame[["Chr","Start","End"]].copy()
    #start coords of bed files are 0 indexed while ends are 1 indexed
    # TODO The above comment was in the codebase. 
    # This needs to be checked. General note: +/- should be handled once, 
    # somewhere obvious, preferrably. 20220708
    bed_frame.loc[:,"Start"] = bed_frame["Start"]-1 
    bed_frame.to_csv(bedfilename, sep = "\t", header = None, index=None)

def make_peaksbed(peaks_frame):
    """This function converts peaks from to a bed file that can be used to
    display the peak locations on the EPCC browser.  The minimum column 
    requirement for the peaks frame is [Chr,Start,End].  This function
    will return a bed frame with the following columns:  
    [Chr in mm10 or hg19 format,Start,End,name,score]
    """
    # TODO resolve this: Do I need to subtract 1 from start??
    # ^^ that is in current codebase (not a comment from NF development 20220708)

    bed_frame = peaks_frame[["Chr","Start","End"]].copy()
    score_list = [1000]*len(peaks_frame)
    bed_frame["Score"] = score_list
    # generate a list of chromosomes 1-22 (note that python is [start,end) ) and 
    # the sex chromosomes
    chrlist = ["chr"+str(x) for x in range(1,23)] + ["chrX", "chrY"]
    bed_frame["Chr"] = pd.Categorical(bed_frame["Chr"], chrlist)
    bed_frame = bed_frame.sort_values(['Chr','Start','End'])
    return bed_frame

def annotate_peaks_frame(peaks_frame, refGene_filename):
    """This function annotates peaks using pybedtools

    :param peaks_frame:
    :param refGene_filename:

    refGene_filename = '/scratch/ref/rmlab/calling_card_ref/mouse/refGene.mm10.Sorted.bed'
    """
    #convert peaks frame to sorted bed
    temp_filename = "temp_peaks_"+str(random.randint(1,1000))+".bed"
    peaks_frame_to_bed(peaks_frame,temp_filename)
    peaks_bed = bt.BedTool(temp_filename)
    peaks_bed = peaks_bed.sort()    
    peaks_bed = peaks_bed.closest(refGene_filename,D="ref",t="first",k=2)
    peaks_bed.saveas(temp_filename)
    #read in gene_annotation_bedfilename    
    temp_annotated_peaks_columns = ["Chr","Start","End","Feature Chr", 
                                    "Feature Start", "Feature End", 
                                    "Feature sName","Feature Name","Strand",
                                    "Distance"]
    temp_annotated_peaks = pd.read_csv(temp_filename, 
                                       sep = "\t", 
                                       names=temp_annotated_peaks_columns)


    # convert start coords back to 1 indexed coordinates
    temp_annotated_peaks.loc[:,"Start"] = temp_annotated_peaks["Start"] + 1 
    index_list = [(x,y,z) for x,y,z in zip(temp_annotated_peaks["Chr"],
                  temp_annotated_peaks["Start"],temp_annotated_peaks["End"])]
    
    try:
        temp_annotated_peaks.index = pd.MultiIndex.from_tuples(index_list)
    except TypeError as e:
        print('Error at line 72: %s' %e)

    index_list = [(x,y,z) for x,y,z in zip(peaks_frame["Chr"], 
                  peaks_frame["Start"],peaks_frame["End"])]

    try:
        peaks_frame.index = pd.MultiIndex.from_tuples(index_list)
    except TypeError as e:
        print('Error at line 80: %s' %e)

    try:
        peaks_frame = peaks_frame.sortlevel(0,axis=1)
    except AttributeError as e:
        print('Error at line 85: %s' %e)

    if "Background Hops" in peaks_frame.columns:
        temp_list = ["Chr","Start","End","Center","Experiment Hops",
                     "Fraction Experiment","TPH Experiment",
        "Background Hops","Fraction Background","TPH Background",
        "TPH Background subtracted","Poisson pvalue","Lambda","Lambda Type",
        "Feature 1 sName","Feature 1 Name",
        "Feature 1 Start","Feature 1 End","Feature 1 Strand","Feature 1 Distance", 
        "Feature 2 sName","Feature 2 Name","Feature 2 Start","Feature 2 End",
        "Feature 2 Strand","Feature 2 Distance"]
    else:
        temp_list = ["Chr","Start","End","Center","Experiment Hops",
                     "Fraction Experiment","TPH Experiment",
        "Poisson pvalue","Lambda","Lambda Type","Feature 1 sName","Feature 1 Name",
        "Feature 1 Start","Feature 1 End","Feature 1 Strand","Feature 1 Distance", 
        "Feature 2 sName","Feature 2 Name","Feature 2 Start","Feature 2 End",
        "Feature 2 Strand", "Feature 2 Distance"]

    peaks_frame = peaks_frame.reindex(columns=temp_list,fill_value=0)

    for idx,row in temp_annotated_peaks.iterrows():
        if not(peaks_frame.loc[idx,"Feature 1 sName"]):
            peaks_frame.loc[idx,["Feature 1 sName","Feature 1 Name",
                                 "Feature 1 Start","Feature 1 End",
                                 "Feature 1 Strand", "Feature 1 Distance"]] = \
                                    list(row[["Feature sName","Feature Name",
                                    "Feature Start", "Feature End","Strand",
                                    "Distance"]])    
        else:
            peaks_frame.loc[idx,["Feature 2 sName","Feature 2 Name",
                                 "Feature 2 Start","Feature 2 End",
                                 "Feature 2 Strand", "Feature 2 Distance"]] = \
                                    list(row[["Feature sName","Feature Name",
                                    "Feature Start", "Feature End","Strand",
                                    "Distance"]])
    return peaks_frame