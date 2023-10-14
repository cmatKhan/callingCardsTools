# Barcode Details Json

The Barcode Details Json provides callingCardsTools with information about 
what and where to expect non-genomic sequence which are included in the raw 
reads and serve to both demultiplex, in the case of yeast data, and confirm 
that a given read represents a transposon insertion.

The mammals data barcode details is generally static -- you can use the same 
barcode details for all data. However, the yeast barcode details must be 
adjusted for each specific library.

## Yeast Barcode Details Format

It is most convenient to store the TF barcode sequences as a `tsv` and then 
use the cmd line tool `callingcardstools barcode_table_to_json` to convert 
the table to the appropriate json template.  

Save a `tsv` in the following format:
**Do not include a header**

|             |       |          |
|-------------|-------|----------|
| MTH1        | GTCCC | CAGAGGGG |
| SKN7        | TCAAG | ATCAGACC |
| HAP3        | AATGA | GGGGGTAG |

The output of 

```bash
callingcardstools barcode_table_to_json -t run_6354_bc_table.tsv -r run_6354
```

will be a json in the following format:

```json
{
    "r1": {
        "primer": {
            "trim": true,
            "index": [
                0,
                5
            ]
        },
        "transposon": {
            "trim": true,
            "index": [
                5,
                22
            ]
        }
    },
    "r2": {
        "transposon": {
            "trim": true,
            "index": [
                0,
                8
            ]
        },
        "restriction": {
            "trim": true,
            "index": [
                8,
                20
            ]
        }
    },
    "components": {
        "r1_transposon": {
            "map": [
                "AATTCACTACGTCAACA"
            ],
            "bam_tag": "RT"
        },
        "r2_restriction": {
            "map": {
                "TCGAGCGCCCGG": "Hpall",
                "TCGAGCGC": "HinP1I",
                "TCGA": "TaqAI"
            },
            "match_type": "greedy",
            "require": false,
            "bam_tag": "RS"
        },
        "tf": {
            "components": [
                "r1_primer",
                "r2_transposon"
            ],
            "map": {
                "GTCCCCAGAGGGG": "MTH1",
                "TCAAGATCAGACC": "SKN7",
                "AATGAGGGGGTAG": "HAP4"
            },
            "bam_tag": "TF"
        }
    },
    "match_allowance": {
        "r1_transposon": 0
    },
    "batch": "run_6354"
}
```

## Mammals Barcode Details Format

Since this will generally be the same for all mammals Calling Cards data, you 
can likely simply copy and paste this onto your system and use it directly:

```json
{
    "batch": "",
    "tf": "",
    "r1": {
        "pb": {"trim": true,
               "index": [0,3]},

        "ltr1": {"trim": true,
                    "index": [3,28]},
        "srt": {"trim": true,
                "index":[28,32]},

        "ltr2": {"trim": true,
                    "index": [32,38]}
    },
    "r2":{},
    "components": {

        "r1_pb":    {"map":["TAG"],
                     "match_allowance": 0,
                     "bam_tag": "PB"},

        "r1_ltr1":  {"map": ["CGTCAATTTTACGCAGACTATCTTT"],
                     "match_type": "edit_distance",
                     "match_allowance": 0,
                     "require": true,
                     "bam_tag": "L1"},
        "r1_srt":   {"map": ["CTAG", "CAAC", "CTGA", "GCAT", "GTAC", "CACA", "TGAC", "GTCA",
                             "CGAT", "CTCT", "GAAG", "TCGA", "CATG", "GTTG", "CTTC", "GCTA",
                             "GAGA", "GTGT", "CGTA", "TGGT", "GGAA", "ACAC", "TCAG", "TTGG",
                             "CAGT", "TTTT"],
                     "match_type": "edit_distance",
                     "match_allowance": 0,
                     "require": true,
                     "bam_tag": "ST",
                     "annotation": true},
        "r1_ltr2":  {"map": ["GGTTAA"],
                     "match_type": "edit_distance",
                     "match_allowance": 0,
                     "require": true,
                     "bam_tag": "L2"}
    },
    "insert_seq": ["TTAA"],
    "max_mismatch": 0
}
```