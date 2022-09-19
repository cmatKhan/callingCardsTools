import os
import json

from numpy import setdiff1d

class BarcodeExtractor:
    json = ""
    barcode_dict = {}
    def __init__(self, barcode_details_json = ""):
        if barcode_details_json:
            self.verify(barcode_details_json, set = True)

    
    def verify(self, barcode_details_json, set = False):
        if not os.path.exists(barcode_details_json):
            raise FileNotFoundError(f"barcode details json DNE: {barcode_details_json}")
        with open(barcode_details_json) as f1:
            barcode_dict = json.load(f1)
        
        if len(setdiff1d(list(barcode_dict.keys()),['r1', 'r2'])) > 0:
            KeyError(f"Only two top level keys allowed, 'r1' and 'r2'. Current " \
                f"barcode dict keys {','.join(list(barcode_dict.keys()))} are invalid")
        
        update_dict = {"r1": {}, "r2": {}}
        for read,read_bc_components in barcode_dict.items():
            index = [0,0]
            for component,details in read_bc_components.items():
                if details.get('indicies') or len(details.get('indicies', [])) != 2:
                    AttributeError(f"{read}.{component} does not have the required field 'index', or it is malformed. Update json and resubmit")
                elif details.get('indicies')[0] < index[0] or details.get('indicies')[1] < index[1]:
                    AttributeError(f"Each subsequent barcode component must be larger than the previous, starting at [0,0]. " \
                    f"Fields in {read}.{component} violate this requirement. Fix and resubmit")
                
                if details.get('trim', None) is None:
                    print(f"'trim' is not set for {read}.{component}. Setting "\
                        f"to False.")
                    update_dict[read][component]['trim'] = False

                if details.get('append', None) is None:
                    print(f"'append' is not set for {read}.{component}. Setting "\
                        f"to True.")
                    update_dict[read][component]['append'] = True

        for read,component in update_dict.items():
            for key,value in component.items():
                barcode_dict[read][component][key] = value

        self.json = barcode_details_json
        self.barcode_dict = barcode_dict

    def process_read(self, read, end):
        if end not in ['r1', 'r2']:
            raise ValueError(f"{end} not valid. Must be either 'r1' or 'r2'")
        
        left_offset = 0
        bc_seq = ""

        for k,v in self.barcode_dict[end].items():
            left_index = v['index'][0] - left_offset
            right_index = v['index'][1] - left_offset
            seq = read.seq[left_index:right_index]
            if v['append']:
                bc_seq += seq
            if v.get("trim"):
                # trim seq
                read = read[:left_index] + read[right_index:]
                # adjust left_offset
                left_offset += v['index'][1]
        # set name to empty string
        return read, bc_seq
        
        


