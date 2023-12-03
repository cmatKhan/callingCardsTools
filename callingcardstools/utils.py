import warnings
import os

__all__ = ['check_file_path', 'map_novo', 'remove_suffix', 'deprecated']


def check_file_path(path: str) -> str:
    """check if a file exists

    Args:
        path (str): path to file

    Raises:
        FileExistsError: if file does not exist
    """
    if os.path.isfile(path):
        return path
    else:
        raise FileExistsError(f"{path} does not exist")


def map_novo():
    novo_str = ('\n').join(
        ['#!/bin/bash\n', '#SBATCH -n 8', '#SBATCH -N 1', '#SBATCH --mem=30000',
         '#SBATCH -o novo_map.log', '#SBATCH -e novo_map.err', '#SBATCH -J novo_map',
         '\nmodule use /opt/htcf/modules', 'ml novoalign', 'eval $(spack load --sh samtools@1.13)\n',
         'output_dir="$2"\n', 'read fastq_file tf_name r1_trim_seq < <(sed -n ${SLURM_ARRAY_TASK_ID}p "$1")\n',
         'novoalign \\',
         '\t-o SAM \\',
         '\t-d /scratch/mblab/chasem/calling_cards/nf_pipeline/yeast_rob_genome/refactored_names_genome/sacCer3_plasmids_minus_Adh1.nix \\',
         '\t-f ${fastq_file} \\',
         '\t-5 ${r1_trim_seq} \\',
         '\t-n 102 2> ${output_dir}/${tf_name}_novoalign.log |\\',
         'samtools view -Sb -q 10 |\\',
         'novosort -i -o ${output_dir}/${tf_name}.bam - 2> ${tf_name}_novosort.log'])

    with open("map_novo_spack.sh", "w") as f:
        f.write(novo_str)


def remove_suffix(input_string: str, suffix: str) -> str:
    """a mimic of python 3.9 removesuffix

    Args:
            input_string (str): _description_
            suffix (str): _description_

    Returns:
            str: _description_
    """
    if suffix and input_string.endswith(suffix):
        return input_string[:-len(suffix)]
    return input_string


def deprecated(message: str):
    """cite: https://stackoverflow.com/a/48632082/9708266 

    Args:
        message (str): _description_
    """
    def deprecated_decorator(func):
        def deprecated_func(*args, **kwargs):
            warnings.warn(
                "{} is a deprecated function. {}".format(
                    func.__name__, message),
                category=DeprecationWarning,
                stacklevel=2)
            warnings.simplefilter('default', DeprecationWarning)
            return func(*args, **kwargs)
        return deprecated_func
    return deprecated_decorator
