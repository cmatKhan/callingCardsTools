# Mapping with novoalign

```bash
sbatch --array=1-5 map_novo_spack.sh map_novo_lookup.txt novo_out
```

where map_novo_lookup.txt looks like

```raw
split_fastq/E0001_RIM101_JP016_TGATA_ACCTGTTT_R1.fastq	RIM101	TGATAAATTCACTACGTCAACA
split_fastq/E0001_RTG1_JP016_TACTC_CAGAGGGG_R1.fastq	RTG1	TACTCAATTCACTACGTCAACA
split_fastq/E0001_ERT1_JP016_CCTGC_ACGCACGC_R1.fastq	ERT1	CCTGCAATTCACTACGTCAACA
split_fastq/E0001_MTH1_JP016_TCGTC_TACGGCGT_R1.fastq	MTH1	TCGTCAATTCACTACGTCAACA
split_fastq/E0001_MET32_JP016_AACGC_TTTTTCTA_R1.fastq	MET32	AACGCAATTCACTACGTCAACA
```

Where the last characters will always be:

```raw
AATTCACTACGTCAACA
```

# make ccf files

```bash
#!/bin/bash

#SBATCH --mem=20000
#SBATCH -o legacy_makeccf.log
#SBATCH -e legacy_makeccf.err
#SBATCH -J legacy_makeccf

source env/bin/activate

read bam < <(sed -n ${SLURM_ARRAY_TASK_ID}p "$1")

legacy_makeccf -s $bam -o "$2"

```

```bash
sbatch --array=1-5 legacy_makeccf_batch.sh bam_lookup.txt cctools_legacy_ccf
```