import sys

def map_novo():
	return ('\n').join(['#!/bin/bash\n', '#SBATCH -n 8','#SBATCH -N 1', '#SBATCH --mem=30000', 
	'#SBATCH -o novo_map.log', '#SBATCH -e novo_map.err', '#SBATCH -J novo_map',	
	'\nmodule use /opt/htcf/modules', 'ml novoalign','eval $(spack load --sh samtools@1.13)\n',
    'output_dir="$2"\n','read fastq_file tf_name r1_trim_seq < <(sed -n ${SLURM_ARRAY_TASK_ID}p "$1")\n',
	'novoalign \\',
	'\t-o SAM \\',
	'\t-d /scratch/mblab/chasem/calling_cards/nf_pipeline/yeast_rob_genome/refactored_names_genome/sacCer3_plasmids_minus_Adh1.nix \\',
	'\t-f ${fastq_file} \\',
	'\t-5 ${r1_trim_seq} \\',
	'\t-n 102 2> ${output_dir}/${tf_name}_novoalign.log |\\',
	'samtools view -Sb -q 10 |\\',
	'novosort -i -o ${output_dir}/${tf_name}.bam - 2> ${tf_name}_novosort.log'])


def main(args=None):
	with open("map_novo_spack.sh","w") as f:
		f.write(map_novo())

if __name__ == "__main__":
    sys.exit(main())
