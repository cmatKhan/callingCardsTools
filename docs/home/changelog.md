# Change Log

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