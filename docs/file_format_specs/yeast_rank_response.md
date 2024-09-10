# Yeast Rank Response Configuration

The `yeast_rank_response` tool requires a json configuration file as input.
The configuration file specifies the location of the input files and various
other settings.  

Below, you will find first an [example](#example) and next a table defining
each key. In the [definitions](#definitions) table, if a key/value pair is
not marked as [REQUIRED], then there is a default which is set in the tool.
[REQUIRED] key/value pairs must be present.

## Example

```json
{
  "binding_data_path": "path/to/binding/data",
  "binding_source": "eg, 'chipexo_1234' or cc_1234'",
  "binding_identifier_col": "identifier_column_name",
  "binding_effect_col": "effect_column_name_or_none",
  "binding_pvalue_col": "pvalue_column_name_or_none",
  "rank_by_binding_effect": false,
  "expression_data_path": "path/to/expression/data",
  "binding_source": "eg, 'mcisaac_1234' or kemmeren_1234'",
  "expression_identifier_col": "identifier_column_name",
  "expression_effect_col": "effect_column_name_or_none",
  "expression_effect_thres": 0.0,
  "expression_pvalue_col": "pvalue_column_name_or_none",
  "expression_pvalue_thres": 0.05,
  "rank_bin_size": 5,
  "normalize": false,
  "output_file": "rank_response.csv",
  "compress": false
}
```

## Definitions

| Key                          | Description |
|------------------------------|-------------|
| `binding_data_path`          | [REQUIRED] Path to the binding data file. The `binding_effect_col`, `binding_pval_col`, and 'gene_id' are required.|
| `binding_source`          | [REQUIRED] A description of where the binding data comes from. This might just be the data source, eg 'chipexo', or it might be a identifier, eg callingcards_17 or cc_17|
| `binding_identifier_col`     | [REQUIRED] Name of the feature identifier column in the binding data.|
| `binding_effect_col`         | [REQUIRED] Name of the effect column in the binding data. Set to `none` if an effect column does not exist.|
| `binding_pvalue_col`         | [REQUIRED] Name of the pvalue column in the binding data. Set to `none` if a pvalue column does not exist.|
| `rank_by_binding_effect`             | `true` or `false`. Defaults to `false` if not provided. Set to `true` to rank by the binding effect column.|
| `expression_data_path`       | [REQUIRED] Path to the expression data file.|
| `expression_source`       | [REQUIRED] A description of where the expressin data comes from. This might just be the data source, eg 'mcisaac', or it might be a identifier, eg mcisaac_17|
| `expression_identifier_col`  | [REQUIRED] Name of the feature identifier column in the expression data.|
| `expression_effect_col`      | [REQUIRED] Name of the effect column in the gene expression data. Set to `none` if an effect column does not exist.|
| `expression_effect_thres`    | [REQUIRED] Threshold for effect expression. Set to `none` if an effect column does not exist|
| `expression_pvalue_col`      | [REQUIRED] Name of the pvalue column in the gene expression data. Set to `none` if a pvalue column does not exist.|
| `expression_pvalue_thres`    | [REQUIRED] Threshold for pvalue of effect expression. Set to `none` if a pvalue column does not exist|
| `rank_bin_size`              | Defaults to 5 if not provided. Bin size for rank response.|
| `normalize`                  | This is not currently implemented -- it is a placeholder for future development when list input of binding/effect data is supported. `true` or `false`.|
| `output_file`                | Path to the output file. Defaults to `rank_response.csv`.|
| `compress`                   | Set this flag to gzip the output file. Defaults to `false` |
