from callingcardstools.BarcodeParser import BarcodeParser


class QcParser(BarcodeParser):

    def __init__(self, id_to_barcode_tsv: str,
                 barcode_details_json: str,
                 **kwargs):
        super().__init__(barcode_details_json=barcode_details_json,
                         **kwargs)


def add_batch_qc(self,
                 barcode_details: BarcodeParser,
                 id_to_bc_map_path: str = None,
                 split_char="x") -> int:
    """_summary_

    Args:
        barcode_details (BarcodeDetails): _description_
        id_to_bc_map_path (str, optional): _description_. Defaults to None.
        split_char (str, optional): _description_. Defaults to "x".

    Returns:
        int: _description_
    """
    # extract data from the barcode_details
    tfs = list(
        barcode_details.barcode_dict['components']['tf']['map'].values())
    batch = barcode_details.barcode_dict['batch']

    parsed_tfs = []
    replicates = []
    for x in tfs:
        try:
            tf = x.split(split_char)[0].replace(')', '')
            rep = x.split(split_char)[1].replace(')', '')
        except IndexError:
            tf = x
            rep = 'none'
        parsed_tfs.append(tf)
        replicates.append(rep)

    parsed_tfs.append('undetermined')
    replicates.append('none')

    # create df with nrow == len(tfs). batch is repeated for each record.
    # records is either split from the tf name, or 'none'
    df = pd.DataFrame(
        {'batch': [batch for x in parsed_tfs],
         'tf': parsed_tfs,
         'replicate': replicates})

    # by row, check if the entry already exists. if it does, then insert.
    # else, skip
    add_id_map = False
    insert_sql = "INSERT INTO batch (batch,tf,replicate) VALUES ('%s','%s','%s')"
    for idx, row in df.iterrows():
        try:
            self.get_batch_id(row['batch'], row['tf'], row['replicate'])
        except IndexError:
            add_id_map = True
            cur = self.con.cursor()
            self.db_execute(cur, insert_sql %
                            (row['batch'], row['tf'], row['replicate']))
            self.con.commit()

    # note that a trigger creates which creates default entries in all other
    # qc tables

    if id_to_bc_map_path and add_id_map:
        self.add_read_qc(barcode_details, id_to_bc_map_path, split_char)


def _update_qc_table(self, tablename: str, id_col: str, id_value: str, update_dict: dict) -> None:
    """_summary_

    Args:
        tablename (str): _description_
        id_col (str): _description_
        id_value (str): _description_
        update_dict (dict): _description_
    """
    for record in update_dict:
        sql = f"UPDATE {tablename} SET tally = {record.get('tally')} "\
            f"WHERE {id_col} = {id_value} AND "\
            f"edit_dist = {record.get('edit_dist')}"
        cur = self.con.cursor()
        self.db_execute(cur, sql)
        self.con.commit()


def _summarize_r1_primer(self, id_to_bc_df: pd.DataFrame, r1_primer: str, r2_transposon: str) -> pd.DataFrame:
    """group by a given r1_primer and return the number of r2_transposon seq within 0 to 4 edit distance 
    of the expected r2_transposon seq

    Args:
        id_to_bc_df (pd.DataFrame): _description_
        r1_primer (str): _description_
        r2_transposon (str): _description_

    Returns:
        pd.DataFrame: dataframe with 5 rows
    """

    primer_df = id_to_bc_df[id_to_bc_df.r1_primer == r1_primer]
    summarised_dict = {}
    if not primer_df.shape[0] == 0:
        primer_df = primer_df\
            .assign(edit_dist=primer_df
                    .apply(lambda x:
                           edlib.align(
                               x['r2_transposon'],
                               r2_transposon)['editDistance'],
                           axis=1))

        summarised_dict = primer_df\
            .sort_values('edit_dist')\
            .groupby('edit_dist')[['edit_dist']]\
            .count()\
            .rename(columns={'edit_dist': 'tally'})\
            .loc[:10, :]\
            .to_dict(orient='index')

    output_dict = {'edit_dist': [], 'tally': []}
    for edit_dist in range(5):
        output_dict['edit_dist'].append(edit_dist)
        output_dict['tally'].append(
            summarised_dict.get(edit_dist, {}).get('tally', 0))

    return pd.DataFrame(output_dict)


def _summarize_r2_transposon(self, id_to_bc_df: pd.DataFrame, r1_primer: str, r2_transposon: str) -> pd.DataFrame:
    """group by a given r2_transposon seq and return the number of r1_primer seqs within 0 to 4 edit distance of the 
    expected r1_primer seq

    Args:
        id_to_bc_df (pd.DataFrame): _description_
        r1_primer (str): _description_
        r2_transposon (str): _description_

    Returns:
        pd.DataFrame: dataframe with 5 rows
    """

    trans_df = id_to_bc_df[id_to_bc_df.r2_transposon == r2_transposon]
    summarised_dict = {}
    if not trans_df.shape[0] == 0:
        trans_df = trans_df\
            .assign(edit_dist=trans_df
                    .apply(lambda x:
                           edlib.align(
                               x['r1_primer'],
                               r1_primer)['editDistance'],
                           axis=1))

        summarised_dict = trans_df\
            .sort_values('edit_dist')\
            .groupby('edit_dist')[['edit_dist']]\
            .count()\
            .rename(columns={'edit_dist': 'tally'})\
            .loc[:10, :]\
            .to_dict(orient='index')

    output_dict = {'edit_dist': [], 'tally': []}
    for edit_dist in range(5):
        output_dict['edit_dist'].append(edit_dist)
        output_dict['tally'].append(
            summarised_dict.get(edit_dist, {}).get('tally', 0))

    return pd.DataFrame(output_dict)


def _summarize_tf_to_r1_transposon(self, id_to_bc_df: pd.DataFrame, r1_primer: str, r2_transposon: str, r1_transposon: str) -> pd.DataFrame:
    """_summary_

    Args:
        id_to_bc_df (pd.DataFrame): _description_
        r1_primer (str): _description_
        r2_transposon (str): _description_
        r1_transposon (str): _description_

    Returns:
        pd.DataFrame: _description_
    """

    r1_trans_df = id_to_bc_df[id_to_bc_df.r1_primer == r1_primer]
    r1_trans_df = r1_trans_df[r1_trans_df.r2_transposon == r2_transposon]

    # this is here to handle case in which there are no r1, r2 expected
    # matches
    summarised_dict = {}
    if not r1_trans_df.shape[0] == 0:
        r1_trans_df = r1_trans_df\
            .assign(edit_dist=r1_trans_df
                    .apply(lambda x:
                           edlib.align(
                               x['r1_transposon'],
                               r1_transposon)['editDistance'],
                           axis=1))

        summarised_dict = r1_trans_df\
            .sort_values('edit_dist')\
            .groupby('edit_dist')[['edit_dist']]\
            .count()\
            .rename(columns={'edit_dist': 'tally'})\
            .loc[:10, :]\
            .to_dict(orient='index')

    output_dict = {'edit_dist': [], 'tally': []}
    for edit_dist in range(5):
        output_dict['edit_dist'].append(edit_dist)
        output_dict['tally'].append(
            summarised_dict.get(edit_dist, {}).get('tally', 0))

    return pd.DataFrame(output_dict)
