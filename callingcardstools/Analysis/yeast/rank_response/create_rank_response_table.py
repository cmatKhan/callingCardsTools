import logging

logger = logging.getLogger(__name__)


def create_rank_response_tables(config_dict: dict) -> dict:
    raise NotImplementedError
    # TODO: offer normalization method
    # if normalize:
    #     min_responsive = \
    #         df[(df['effect_expression'].abs() > effect_expression_thres)
    #            & (df['p_expression'] < p_expression_thres)].shape[0]
    # else:
    #     min_responsive = int('inf')

    # try:
    #     binding_data = read_in_data(
    #         args.binding_data_path,
    #         args.identifier_col,
    #         args.binding_effect_col,
    #         args.binding_pval_col,
    #         'binding',
    #         args.binding_source)
    # except (KeyError, FileExistsError) as exc:
    #     logger.error("Error reading in binding data: %s", exc)
    #     raise

    # try:
    #     expression_data = read_in_data(
    #         args.expression_data_path,
    #         args.expression_effect_col,
    #         args.expression_pval_col,
    #         'expression',
    #         args.expression_source)
    # except (KeyError, FileExistsError) as exc:
    #     logger.error("Error reading in expression data: %s", exc)
    #     raise

    # df = expression_data.merge(binding_data[['binding_pvalue', 'gene_id']],
    #                            how='left',
    #                            left_on='gene_id',
    #                            right_on='gene_id')

    # df, random_expectation_df, rank_response_df = \
    #     rank_response_ratio_summarize(
    #         df,
    #         effect_expression_thres=args.effect_expression_thres,
    #         p_expression_thres=args.p_expression_thres,
    #         normalize=args.normalize,
    #         bin_size=args.bin_size)

    # df.to_csv(args.output_prefix + '_rank_response.tsv', sep='\t', index=False)
    # random_expectation_df.to_csv(args.output_prefix +
    #                              '_random_expectation.tsv',
    #                              sep='\t',
    #                              index=False)
    # rank_response_df.to_csv(args.output_prefix +
    #                         '_rank_response.tsv',
    #                         sep='\t',
    #                         compression='gzip' if args.compress else None,
    #                         index=True)
