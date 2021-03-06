#standard lib
import os

# imports from local repo
from callingcardstools import macs_style as ms
from callingcardstools import utils
# third party
import numpy as np
import pandas as pd
from bx.intervals.intersection import Intersecter, Interval

def call_peaks_and_annotate(expframe, bgframe, TTAAframe, outputfile,
                            annotation_file, peak_pvalue_cutoff=1e-4, 
							window_size = 1000, step_size = 500, 
							pseudocounts = 0.2):

	# This function is the wrapper that calls the various subfunctions

	# Args:
	# 	expfile (pandas DataFrame): CCF generated by experiment. columns must be 
	# 	  named the following -- ["Chr","Start","End","Reads","Strand","Barcode"]
	# 	bgfile (pandas DataFrame): Background CCF. columns must be 
	# 	  named the following -- ["Chr","Start","End","Reads","Strand","Barcode"]
	# 	outputfile (_type_): _description_
	# 	TTAA_file (str, optional): _description_. original default to '/scratch/ref/rmlab/calling_card_ref/mouse/TTAA_mm10.txt'.
	# 	annotation_file (str, optional): _description_. original default to '/scratch/ref/rmlab/calling_card_ref/mouse/refGene.mm10.Sorted.bed'.
	# 	pvalue_cutoff (_type_, optional): _description_. Defaults to 1e-4.
	# 	peak_pvalue_cutoff (_type_, optional): _description_. Defaults to 1e-3.
	# 	window_size (int, optional): _description_. Defaults to 1000.
	# 	step_size (int, optional): _description_. Defaults to 500.
	# 	pseudocounts (float, optional): _description_. Defaults to 0.2.
	



	peaks_frame = find_peaks(expframe,
	                         bgframe,
							 TTAAframe,
							 peak_pvalue_cutoff,
							 window_size,
							 step_size,
							 pseudocounts)

    # add feature names, etc to peaks_frame
	peaks_frame = \
		ms.annotate.annotate_peaks_frame(peaks_frame, annotation_file)

    # row filter and sort based on poisson_pvalue
	peaks_frame = \
		peaks_frame[peaks_frame["poisson_pvalue"] <= peak_pvalue_cutoff]
	peaks_frame = peaks_frame.sort_values(["poisson_pvalue"])

	peaksbed_frame = create_peaksbed_frame(peaks_frame)

    # TODO figure out what is going on in these files -- the bed is definitely 
	# not a bed file

	# write peaks_frame out to csv
	peaks_frame.to_csv(outputfile, sep="\t", index=False)
    # write peaksbed to csv
	bedfilename = os.path.splitext(outputfile)[0]+"_peaks.bed"
	peaksbed_frame.to_csv(bedfilename,
	                      sep="\t",
						  index=False,
						  header=False)

def create_peaksbed_frame(peaks_frame):
	"""_summary_

	Args:
		peaks_frame (_type_): _description_

	Returns:
		_type_: _description_
	"""

	peaksbed_frame = ms.annotate.make_peaksbed(peaks_frame)
	peaksbed_frame.columns = ['Chr','Start','End','Col4']

	peaksbed_frame_with_description = pd.merge(peaksbed_frame,
	                                           peaks_frame, 
											   how = "left", 
											   on = ['Chr','Start','End'])

	peaksbed_frame_with_description['tph_experiment'] = 'tph_experiment=' + \
		peaksbed_frame_with_description['tph_experiment'].astype(str)

	peaksbed_frame_with_description['TPH Background'] = 'TPH Background=' + \
		peaksbed_frame_with_description['TPH Background'].astype(str)

	peaksbed_frame_with_description['TPH Background subtracted'] = \
		'TPH Background subtracted=' + \
			peaksbed_frame_with_description['TPH Background subtracted'].astype(str)

	peaksbed_frame_with_description['poisson_pvalue'] = 'poisson_pvalue=' + \
		peaksbed_frame_with_description['poisson_pvalue'].astype(str)

	peaksbed_frame_with_description['lambda'] = 'lambda=' + \
		peaksbed_frame_with_description['lambda'].astype(str)

	peaksbed_frame_with_description['lambda_type'] = 'lambda_type=' + \
		peaksbed_frame_with_description['lambda_type'].astype(str)

	peaksbed_frame_with_description['Description'] = \
		list(peaksbed_frame_with_description.loc[:,['tph_experiment',
		                                            'TPH Background',
													'TPH Background subtracted',
													'poisson_pvalue',
													'lambda',
													'lambda_type']].values)
	
	final_peaksbed_frame = \
		peaksbed_frame_with_description.loc[:,['Chr','Start','End',
		                                       'Col4','Description']]

	return final_peaksbed_frame



def find_peaks(experiment_frame,background_frame,TTAA_frame,
	pvalue_cutoff = 1e-3,window_size = 1000, step_size = 500,
	pseudocounts = 0.2):

	"""This function is passed an experiment frame, a background frame, 
	and an TTAA_frame,all in ccf format:
	(Chr, Start, Stop, Reads, Strand, Barcode)

	It then builds interval trees containing all of the background and 
	experiment hops and all of the TTAAs.  Next, it scans through the 
	genome with a window of window_size and step size of step_size and 
	looks for regions that have signficantly more experiment hops than
	background hops (poisson w/ pvalue_cutoff). It merges consecutively 
	enriched windows and computes the center of the peak.  Next 
	it computes lambda, the number of insertions per TTAA expected from the
	background distribution by taking the max of lambda_bg, lamda_1, lamda_5,
	lambda_10.  It then computes a p-value based on the expected number of hops 
	= lamda * number of TTAAs in peak * number of hops in peak.  Finally, it 
	returns a frame that has 
	[Chr,Start,End,Center,experiment_hops,
	fraction_experiment,Background Hops,Fraction Background,poisson_pvalue]

	"""
	peaks_frame = pd.DataFrame(
		columns = ["Chr","Start","End","Center","experiment_hops",
		"fraction_experiment","tph_experiment","Background Hops",
		"Fraction Background","TPH Background","TPH Background subtracted",
		"lambda_type","lambda","poisson_pvalue"])

	
	experiment_gnashy_dict = {}
	experiment_dict_of_trees = {}
	total_experiment_hops = len(experiment_frame)
	background_gnashy_dict = {} 
	background_dict_of_trees = {}
	total_background_hops = len(background_frame)
	TTAA_frame_gbChr_dict = {} 
	TTAA_dict_of_trees = {}
	list_of_l_names = ["bg","1k","5k","10k"]

	#group by chromosome and populate interval tree with TTAA positions
	for name,group in TTAA_frame.groupby('Chr'):
		 
		TTAA_frame_gbChr_dict[name] = group
		TTAA_frame_gbChr_dict[name].index = TTAA_frame_gbChr_dict[name]["Start"]
		#initialize tree
		TTAA_dict_of_trees[name] = Intersecter() 
		#populate tree with position as interval
		for idx, row in TTAA_frame_gbChr_dict[name].iterrows():	
			TTAA_dict_of_trees[name].add_interval(Interval(int(idx),int(idx+3)))
	#group by chromosome and populate interval tree with positions of background hops
	for name,group in background_frame.groupby('Chr'):
		background_gnashy_dict[name] = group
		background_gnashy_dict[name].index = background_gnashy_dict[name]["Start"]
		#initialize tree
		background_dict_of_trees[name] = Intersecter() 
		#populate tree with position as interval
		for idx, row in background_gnashy_dict[name].iterrows():	
			background_dict_of_trees[name].add_interval(Interval(int(idx),int(idx)+3)) 

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
	num_bg_hops_list = []
	frac_exp_list = []
	tph_exp_list = []
	frac_bg_list = []
	tph_bg_list = []
	tph_bgs_list = []
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
		print(max_pos)
		print(window_size)
		print(step_size)
		for window_start in range(1,max_pos+window_size,step_size):
			overlap = experiment_dict_of_trees[name]\
				.find(window_start,window_start+window_size - 1)
			num_exp_hops = len(overlap)
			bg_overlap = background_dict_of_trees[name]\
				.find(window_start,window_start+window_size - 1)
			num_bg_hops = len(bg_overlap)

			#is this window significant?
			if utils.compute_cumulative_poisson(num_exp_hops,
			                              num_bg_hops,
										  total_experiment_hops,
										  total_background_hops,
										  pseudocounts) < pvalue_cutoff:
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
						.append(float(num_exp_hops)*100000/total_experiment_hops)


					#add number of background hops in peak to frame
					bg_overlap = background_dict_of_trees[name]\
						.find(sig_start,sig_end)
					num_bg_hops = len(bg_overlap)
					num_bg_hops_list.append(num_bg_hops)

					frac_bg_list\
						.append(float(num_bg_hops)/total_background_hops)
					tph_bg_list\
						.append(float(num_bg_hops)*100000/total_background_hops)		
					#find lambda and compute significance of peak
					if total_background_hops >= total_experiment_hops: #scale bg hops down
						#compute lambda bg
						num_TTAAs = \
							len(TTAA_dict_of_trees[name]\
								.find(sig_start,sig_end))
						lambda_bg = \
							((num_bg_hops*(float(total_experiment_hops)/\
								total_background_hops))/max(num_TTAAs,1))


						#compute lambda 1k
						num_bg_hops_1k = \
							len(background_dict_of_trees[name]\
								.find(peak_center-499,peak_center+500))
						num_TTAAs_1k = \
							len(TTAA_dict_of_trees[name]\
								.find(peak_center-499,peak_center+500))
						lambda_1k = \
							(num_bg_hops_1k*(float(total_experiment_hops)/\
								total_background_hops))/(max(num_TTAAs_1k,1))


						#compute lambda 5k
						num_bg_hops_5k = \
							len(background_dict_of_trees[name]\
								.find(peak_center-2499,peak_center+2500))
						num_TTAAs_5k = \
							len(TTAA_dict_of_trees[name]\
								.find(peak_center-2499,peak_center+2500))
						lambda_5k = \
							(num_bg_hops_5k*(float(total_experiment_hops)/\
								total_background_hops))/(max(num_TTAAs_5k,1))


						#compute lambda 10k
						num_bg_hops_10k = \
							len(background_dict_of_trees[name]\
								.find(peak_center-4999,peak_center+5000))
						num_TTAAs_10k = \
							len(TTAA_dict_of_trees[name]\
								.find(peak_center-4999,peak_center+5000))
						lambda_10k = \
							(num_bg_hops_10k*(float(total_experiment_hops)/\
								total_background_hops))/(max(num_TTAAs_10k,1))
						lambda_f = \
							max([lambda_bg,lambda_1k,lambda_5k,lambda_10k])


						#record type of lambda used
						index = \
							[lambda_bg,lambda_1k,lambda_5k,lambda_10k]\
								.index(max([lambda_bg,lambda_1k,
								            lambda_5k,lambda_10k]))
						lambda_type_list.append(list_of_l_names[index])
						#record lambda
						lambda_list.append(lambda_f)
						#compute pvalue and record it

						pvalue = 1-scistat.poisson.cdf(
							(num_exp_hops+pseudocounts),
							lambda_f*max(num_TTAAs,1)+pseudocounts)
						pvalue_list.append(pvalue)
						

						tph_bgs = \
							float(num_exp_hops)*100000/total_experiment_hops - \
								float(num_bg_hops)*100000/total_background_hops
						tph_bgs_list.append(tph_bgs)

						index = [lambda_bg,lambda_1k,lambda_5k,lambda_10k]\
							.index(max([lambda_bg,lambda_1k,
							            lambda_5k,lambda_10k]))
						lambdatype = list_of_l_names[index]
						l = [pvalue,tph_bgs,lambda_f,lambdatype]

					else: #scale experiment hops down
						#compute lambda bg
						num_TTAAs = len(TTAA_dict_of_trees[name]\
							.find(sig_start,sig_end))
						lambda_bg = (float(num_bg_hops)/max(num_TTAAs,1)) 



						#compute lambda 1k
						num_bg_hops_1k = len(background_dict_of_trees[name]\
							.find(peak_center-499,peak_center+500))
						num_TTAAs_1k = len(TTAA_dict_of_trees[name]\
							.find(peak_center-499,peak_center+500))
						lambda_1k = (float(num_bg_hops_1k)/(max(num_TTAAs_1k,1)))



						#compute lambda 5k
						num_bg_hops_5k = len(background_dict_of_trees[name]\
							.find(peak_center-2499,peak_center+2500))
						num_TTAAs_5k = len(TTAA_dict_of_trees[name]\
							.find(peak_center-2499,peak_center+2500))
						lambda_5k = (float(num_bg_hops_5k)/(max(num_TTAAs_5k,1)))



						#compute lambda 10k
						num_bg_hops_10k = len(background_dict_of_trees[name]\
							.find(peak_center-4999,peak_center+5000))
						num_TTAAs_10k = len(TTAA_dict_of_trees[name]\
							.find(peak_center-4999,peak_center+5000))
						lambda_10k = (float(num_bg_hops_10k)/(max(num_TTAAs_10k,1)))
						lambda_f = max([lambda_bg,lambda_1k,lambda_5k,lambda_10k])



						#record type of lambda used
						index = [lambda_bg,lambda_1k,lambda_5k,lambda_10k]\
							.index(max([lambda_bg,lambda_1k,lambda_5k,lambda_10k]))
						lambda_type_list.append(list_of_l_names[index])
						#record lambda
						lambda_list.append(lambda_f)
						#compute pvalue and record it
						pvalue = 1-scistat.poisson.cdf(
							((float(total_background_hops)/total_experiment_hops)*\
								num_exp_hops+pseudocounts),
								lambda_f*max(num_TTAAs,1)+pseudocounts)
						pvalue_list.append(pvalue)

						tph_bgs = float(num_exp_hops)*100000/total_experiment_hops - \
							float(num_bg_hops)*100000/total_background_hops
						tph_bgs_list.append(tph_bgs)

						index = [lambda_bg,lambda_1k,lambda_5k,lambda_10k]\
							.index(max([lambda_bg,lambda_1k,lambda_5k,lambda_10k]))
						lambdatype = list_of_l_names[index]
						l = [pvalue,tph_bgs,lambda_f,lambdatype]

					#number of hops that are a user-defined distance from peak center
					sig_flag = 0

					

				#else do nothing.		
				
	#make frame from all of the lists
	peaks_frame["Chr"] = chr_list
	peaks_frame["Start"] = start_list
	peaks_frame["End"] = end_list
	peaks_frame["Center"] = center_list
	peaks_frame["experiment_hops"] = num_exp_hops_list 
	peaks_frame["fraction_experiment"] = frac_exp_list 
	peaks_frame["tph_experiment"] = tph_exp_list
	peaks_frame["Background Hops"] = num_bg_hops_list 
	peaks_frame["Fraction Background"] = frac_bg_list
	peaks_frame["TPH Background"] = tph_bg_list
	peaks_frame["TPH Background subtracted"] = tph_bgs_list
	peaks_frame["lambda_type"] = lambda_type_list
	peaks_frame["lambda"] = lambda_list
	peaks_frame["poisson_pvalue"] = pvalue_list

	return peaks_frame

