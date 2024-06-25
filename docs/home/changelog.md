# Change Log

# version 1.6.0

### Changes

- enrichment_vectorized has been significantly restructured such that the pseudocunt
  is removed from the final calculation of the denominator. If the denominator is
  0, enrichment is now set to -1, which signifies no enrichment (enrichment cannot
  be negative, so should clue a user into the fact that they should check
  documentation on the meaning).
- poisson_pvalue has been significantly restructured such that now no pseudocount is
  used in a calculation where the poisson is defined. Where the poisson is not defined,
  the pvalue is set to 1. This is consistent with the hypergeometric pvalue.
- Additional input checking:
    1. chrmap nrow must be > 0 and ncol > 1
    1. background data (hops) must have at least one hop. Will error if not.
    1. experiment data (hops) must have at least one hop. Will error if not.
    1. promoter set files must have at least 1 promoter region
- Old un-vectorized enrichment, poisson and hypergeometric code files removed

# #Version 1.5.2

### Changes

- added kwargs arguments to PeakCalling.yeast.call_peaks to allow user
  to pass in validation method on pyranges join, background_total_hops and
  experiment_total_hops.
- moved promoters_df to promoters_pr conversion in PealCalling.yeast.call_peaks
  from call_peaks to external function. Also corrected the `slack` in the join method
  where the overlaps are counted. Now in the conversion method, the End is incremented
  by 1 to allow hops on the right endpoint, whatever that is, to be counted.

## Version 1.5.1

### Changes

- Needed to keep `name` in the output of PeakCalling.yeast.call_peaks
- adding Analysis and PeakCalling modules to the documentation API section

## Version 1.5.0

### Changes

- overhaul of the PeakCalling/yeast module to address memory usage.
  adding pyranges as a depedency as a result. removed `consider_strand`
  and added a argument to deduplicate the experiment qbeds based on
  `chr`, `start`, `end`

## Version 1.4.1

### Changes

- chipexo_promoter_sig now checks the columns and expects the original
  yeastepigenome.org allevents `coord` column to be split into `start` 
  and `end` where `end` is simply coord + 1. This is to be consistent
  with the other bed-type files.

## Version 1.4.0

### Additions

For yeast, changing the `yeast_call_peaks` `consider_strand` functionality
to collapse read counts at the same coordinate on the forward/reverse strand
in addition to ignoring the strand with regards to the promoter.

## Version 1.3.0

### Additions

For yeast, adding peak calling and analysis functionality. This includes
the following to the cmd line:
    - `yeast_call_peaks` - call peaks on a yeast qBED file
    - `yeast_chipexo_promoter_sig` - calculate the significance of the
    promoter signal over a given set of promoter regions
    - `yeast_rank_response` - Given a promoter set and expression data, compare
    the binding and expression sets
  
### Changes

- the yeast barcode qc summary outputs the actual r1/r2 sequences as opposed to
just the edit distance equivalent classes

## Version 1.2.0

### Bug fixes

- Removing the deprecated (and removed in pandas 2.0) DataFrame.append
  function calls from Qbed.py

## Version 1.1.0

### Bug fixes

- Adding SRT annotation to mammals qBED output

### Features/not bugs

- Remove header from mammals qBed
- Fix typo in mammals barcode details -- the component identified with `lrt`
should be (and now is) `ltr`

## Version 1.0.0

- Initial release