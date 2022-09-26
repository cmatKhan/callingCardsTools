"""ReadParser will examine single or paired end reads for the expected read components given by a barcode details json


Todo:
    * document each function
"""
from callingcardstools.bam_parsers import BarcodeParser
from Bio import SeqRecord

class ReadParser(BarcodeParser):
    """Given either single or paired end reads, use the provided barcode details json to examine expected read components

    Depending on the entries in the barcode details json, this class will parse 
    the read(s), return the assembled components (this could include both what 
    could be construed as a barcode as well as any other components), and the 
    reads trimmed for the components labelled for trimming.
    
    Example:
    ```
    >>> rb = ReadParser('/path/to/barcode_details.json')
    >>> r1 = next(SeqIO.parse('/path/to/r1.fq', format="fastq"))
    >>> r2 = next(SeqIO.parse('/path/to/r2.fq', format="fastq"))
    >>> read_dict = rb.process_read(r1,r2)
    ```
    """

    def __init__(self, barcode_details_json:str = "") -> None:
        """_summary_

        Args:
            barcode_details_json (str, optional): _description_. Defaults to "".
        """
        if barcode_details_json:
            super().__init__(barcode_details_json)

    
    # def verify(self, barcode_details_json: str) -> None:
    #     """_summary_

    #     Args:
    #         barcode_details_json (str): _description_

    #     Raises:
    #         FileNotFoundError: _description_
    #     """
    #     if not os.path.exists(barcode_details_json):
    #         raise FileNotFoundError(f"barcode details json DNE: {barcode_details_json}")
    #     with open(barcode_details_json) as f1:
    #         barcode_dict = json.load(f1)
        
    #     if len(setdiff1d(list(barcode_dict.keys()),['r1', 'r2'])) > 0:
    #         KeyError(f"Only two top level keys allowed, 'r1' and 'r2'. Current " \
    #             f"barcode dict keys {','.join(list(barcode_dict.keys()))} are invalid")
        
    #     update_dict = {"r1": {}, "r2": {}}
    #     for read,read_bc_components in barcode_dict.items():
    #         index = [0,0]
    #         for component,details in read_bc_components.items():
    #             if details.get('indicies') or len(details.get('indicies', [])) != 2:
    #                 AttributeError(f"{read}.{component} does not have the required field 'index', or it is malformed. Update json and resubmit")
    #             elif details.get('indicies')[0] < index[0] or details.get('indicies')[1] < index[1]:
    #                 AttributeError(f"Each subsequent barcode component must be larger than the previous, starting at [0,0]. " \
    #                 f"Fields in {read}.{component} violate this requirement. Fix and resubmit")
                
    #             if details.get('trim', None) is None:
    #                 print(f"'trim' is not set for {read}.{component}. Setting "\
    #                     f"to False.")
    #                 update_dict[read][component]['trim'] = False

    #             if details.get('append', None) is None:
    #                 print(f"'append' is not set for {read}.{component}. Setting "\
    #                     f"to True.")
    #                 update_dict[read][component]['append'] = True

    #     for read,component in update_dict.items():
    #         for key,value in component.items():
    #             barcode_dict[read][component][key] = value

    #     self.json(barcode_details_json) #pylint:disable=E1102
    #     self.barcode_dict(barcode_dict) #pylint:disable=E1102

    def process_read(self, r1:SeqRecord, r2:SeqRecord=None, add_bc_to_id:bool = True, set_name_to_empty:bool = True) -> dict:
        """_summary_

        Args:
            r1 (SeqRecord): _description_
            r2 (SeqRecord, optional): _description_. Defaults to None.
            add_bc_to_id (bool, optional): _description_. Defaults to True.
            set_name_to_empty (bool, optional): _description_. Defaults to True.

        Returns:
            dict: _description_
        """

        if r2 is not None and r1.id != r2.id:
            ValueError(f"r1 ID: {r1.id} does not match r2 ID: {r2.id}")

        # add id to read_dict prior to adding barcode
        read_id = r1.id

        read_dict = {
            'r1': r1,
            'r2': r2
        }
        offset_dict = {
            'r1': 0,
            'r2': 0
        }
        bc_sequence = ""
        for end,read in read_dict.items():
            if read:
                for k,v in self.barcode_dict[end].items():
                    bc_subseq = read.seq[v['index'][0]:v['index'][1]]
                    if v['append']:
                        bc_sequence += str(bc_subseq)
                    # adjust offset for trimming
                    if v['trim']: 
                        left_index = v['index'][0] - (offset_dict[end])
                        right_index = v['index'][1] - (offset_dict[end])
                        read_dict[end] = \
                            read_dict[end][:left_index] + read_dict[end][right_index:]
                        # adjust offset
                        offset_dict[end] += v['index'][1]
        # add bc to id a la UMITools
        if add_bc_to_id:
            for k,v in read_dict.items():
                read_dict[k].id = "_".join([v.id,bc_sequence])
        # set name to empty string if flag set
        if set_name_to_empty:
            r1.name = ""
            r2.name = ""
        # add barcode to read_dict 
        read_dict['bc'] = bc_sequence
        read_dict['id'] = read_id
        return read_dict
        