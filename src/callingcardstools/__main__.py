import sys
from importlib.metadata import version

from .entry_points.split_fastq import main as split_fastq
from .entry_points.legacy_split_fastq import main as legacy_split_fastq
from .entry_points.process_alignments import main as process_alignments

def main(args=None):

	cct_version = version('callingcardstools')

	# Stuff to print before offering a list of tools
	preamble = [
		f'callingcardstools version: {cct_version}\n',
		'Description: A set of tools and API for processing and exploring '+\
			'Calling Cards experiment data\n',
		'Available tools are:\n'
	]

    # list of available tools
	# preface each with \t
	tool_list = [
		'\tsplit_fastq: parse a (possibly multiplexed) batch of reads into\n'+\
			'\texpected barcode file(s), and a set of undetermined reads\n',
		'\tlegacy_split_fastq: yeast cc_tools 3.0 version of parse_fastq\n',
		'\tprocess_alignments: add read and alignment data to a\n'+\
			'\tsqlite database for downstream QC and analysis'
		'\n'
	]

	# Anything to print after the list of tools
	epigraph = [
		'For tool usage instructions, enter callingcardstools <tool_name> --help\n\n'+\
			'For more detailed documentation on the tools, the API, and tutorials\n'+\
				'describing interactive usage, see https://cmatkhan.github.io/callingCardsTools'
	]

	help_str = "\n".join(
		["\n".join(preamble),"\n".join(tool_list),"\n".join(epigraph)])

	try:
		tool = sys.argv[1]
	except IndexError:
		print(help_str)
		sys.exit()
	
	if tool == '--version':
		print("calligncardstools version: "+cct_version)
	elif tool == 'split_fastq':
		split_fastq(sys.argv[2:])
	elif tool == 'legacy_split_fastq':
		legacy_split_fastq(sys.argv[2:])
	elif tool == 'process_alignments':
		process_alignments(sys.argv[2:])
	else:
		print(f"{sys.argv[1]} not recognized! Try callingcardstools "\
			f"--help for a list of available tools")

if __name__ == "__main__":
    sys.exit(main())