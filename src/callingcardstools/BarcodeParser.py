"""An object which parses a barcode string extracted from a read according 
to a barcode_details json file. Note that all sequences are cast to upper"""
import os
import json
import logging 

logging.getLogger(__name__).addHandler(logging.NullHandler())

class BarcodeParser:
    """Using a json which describes acceptable values for given barcodes, 
    check the edit distance between the barcode components (substrings) and 
    the """
    # constructor --------------------------------------------------------------
    def __init__(self, barcode_details_json: str) -> None:
        """BarcodeParser Constructor

        Args:
            barcode_details_json (str): Path to the barcode details json file
        """
        # set the initialized barcode details
        # note that this also sets the match allowances
        self.barcode_details_json(barcode_details_json)
    
    # attributes ---------------------------------------------------------------
    _key_dict = {
        "indicies":"indicies",
        "restriction_enzyme": "restriction_site",
        "components":"components",
        "insert_seqs": "insert_seqs",
        "tf_map": "tf_map",
        "tf_map_keys": {"bc_components": "bc_components", "tf": "tf"},
        # this is created from tf_map, which is provided in barcode_details
        "tf_dict": "tf_dict",
        "match_allowance": "match_allowance",
        "max_mismatch_key": "max",
        "barcode_length": "length",
        "insert_length": "insert_length"
    }
    _barcode_dict = {}
    _barcode_details_json = ""
    _barcode = ""

    # getters/setters ----------------------------------------------------------
    @property
    def barcode_details_json(self):
        """path to the barcode details json file"""
        return self._barcode_details_json
    @barcode_details_json.setter
    def barcode_details_json(self, new_barcode_details_json):
        # check that barcode json exists
        if not os.path.exists(new_barcode_details_json):
            raise FileNotFoundError(f"Invalid Path: {new_barcode_details_json}")
        # open json, read in as dict
        with open(new_barcode_details_json, 'r') as f1: #pylint:disable=W1514
            barcode_dict = json.load(f1)
        # update the barcode_dict attribute
        logging.info("Updating the barcode dict to reflect the new barcode details json...")
        self.barcode_dict(barcode_dict)
        # if that works, update the barcode_details_json path
        logging.info("Success! Setting the barcode details json path")
        self._barcode_details_json = new_barcode_details_json
    
    @property
    def key_dict(self):
        """A dictionary which describes the defined fields of a valid barcode dict json file"""
        return self._key_dict
    
    @property
    def barcode_dict(self):
        """The barcode details json file verified and parsed into a python dictionary"""
        return self._barcode_dict
    @barcode_dict.setter
    def barcode_dict(self, new_barcode_dict):
        # check that the indicies and components match
        if not list(new_barcode_dict[self.key_dict['components']].keys()).sort() == \
            list(new_barcode_dict[self.key_dict['indicies']].keys()).sort():
            raise KeyError('The keys in the barcode_component json and the ' +
                'barcode_component_indicies json are not the same')
        # check that the components are the correct length according to the 
        # indicies
        barcode_components_length_check_dict = { k:v[1]-v[0] for k,v in \
                                                new_barcode_dict['indicies'].items()}
        for k,v in new_barcode_dict['components'].items():
            if isinstance(v, dict):
                # TODO address hard coding 'seq' here
                if not max([len(x) for x in v['seq']]) == \
                    barcode_components_length_check_dict[k]:
                    raise ValueError(f"There exists a component in "\
                    f"barcode_components which is not the length " \
                    f"described by the barcode_details['indicies']: "\
                    f"{k}, {barcode_components_length_check_dict[k]}")
            else:
                for comp in v:
                    if not len(comp) == barcode_components_length_check_dict[k]:
                        raise ValueError(f"There exists a component in "\
                        f"barcode_components which is not the length " \
                        f"described by the barcode_details['indicies']: "\
                        f"{comp}, {barcode_components_length_check_dict[k]}")
        
        # set barcode length
        new_barcode_dict[self.key_dict['barcode_length']] = \
            sum(barcode_components_length_check_dict.values())
        
        # if a tf_map dict exists, check that the bc_components and tf list exist
        if self.key_dict['tf_map'] in new_barcode_dict:
            for bc_comp in new_barcode_dict[self.key_dict['tf_map']]\
                [self.key_dict['tf_map_keys']['bc_components']]:
                if bc_comp not in new_barcode_dict[self.key_dict['components']]:
                    raise KeyError(f"tf barcode component {bc_comp} is not " \
                        f"in barcode details {self.key_dict['components']}")
        # if the match_allowance field exists, check that the components match 
        # those in the components field, and that the distance is an integer
        if self.key_dict['match_allowance'] in new_barcode_dict:
            for comp, dist in new_barcode_dict[self.key_dict['match_allowance']].items():
                if comp not in new_barcode_dict[self.key_dict['components']] and not comp == "max":
                    raise KeyError(f"The keys in "\
                        f"{self.barcode_dict[self.key_dict['match_allowance']]} " \
                            f"must be in {self.barcode_dict[self.key_dict['components']]}")
                if not isinstance(dist, int) or dist < 0:
                    raise ValueError("Match allowance must be a positive integer which \
                        represents the maximum edit distance from the expected \
                            barcode")
        
        # cast all values in self.key_dict['components'] to upper
        # note that the barcode is also always cast to upper
        for k,seq_list in new_barcode_dict[self.key_dict['components']].items():
            if isinstance(seq_list, dict):
                new_barcode_dict[self.key_dict['components']][k]['seq'] = \
                    [x.upper() for x in new_barcode_dict[self.key_dict['components']][k]['seq']]
            elif isinstance(seq_list, list):
                new_barcode_dict[self.key_dict['components']][k] = \
                    [x.upper() for x in seq_list]
            else:
                ValueError(f"Unrecognized element in {k}. Recognized data types are list and dict")

        # cast insert seq to upper if it exists
        if self.key_dict['insert_seqs'] in new_barcode_dict:
            upper_insert_seqs = []
            for sequence in new_barcode_dict[self.key_dict['insert_seqs']]:
                upper_insert_seqs.append(sequence.upper())
            new_barcode_dict[self.key_dict['insert_seqs']] = upper_insert_seqs
        
        # if tf_map exists, then set the tf_dict
        if self.key_dict['tf_map'] in new_barcode_dict:
            new_barcode_dict[self.key_dict['tf_dict']] = self.tf_dict(new_barcode_dict)
        
        # if the match_allowance field DNE, create one with an empty dict
        if self.key_dict['match_allowance'] not in new_barcode_dict:
            new_barcode_dict[self.key_dict['match_allowance']] = {}
        # for each key in the barcode_dict[comp], set to 0 if it does not 
        # already exist
        for comp in new_barcode_dict[self.key_dict['components']]:
            x = new_barcode_dict[self.key_dict['match_allowance']]\
                .setdefault(comp,0)
        # if a 'max' category is not set, set max to the sum of the component 
        # match allowances
        if 'max' not in new_barcode_dict[self.key_dict['match_allowance']]:
            new_barcode_dict[self.key_dict['match_allowance']]\
                [self.key_dict['max_mismatch_key']] = \
                    sum([v for v in \
                        new_barcode_dict[self.key_dict['match_allowance']].values()])
        
        # set insert_seqs if DNE
        if self.key_dict['insert_seqs'] not in new_barcode_dict:
            new_barcode_dict[self.key_dict['insert_seqs']] = ["*"]
        
        # set insertion length
        # TODO check that all insertion seqs are the same length
        new_barcode_dict[self.key_dict['insert_length']] = \
            len(new_barcode_dict[self.key_dict['insert_seqs']][0])
        
        # set barcode dict
        self._barcode_dict = new_barcode_dict
    
    @property
    def barcode_length(self):
        """Extract the barcode_length from the barcode dictionary"""
        return self.barcode_dict[self.key_dict['barcode_length']]

    @property
    def insert_length(self):
        """Extract the insertion sequence length from the barcode dictionary"""
        return self.barcode_dict[self.key_dict['insert_length']]
    
    @property
    def insert_seqs(self):
        """Getter for the insert seq sequence from the barcode details json. Returns upper case.

        Raises:
            AttributeError: Raised if the current barcode details json does 
            not have an insert seq key
        """
        if self.key_dict['insert_seqs'] in self.barcode_dict:
            return self.barcode_dict[self.key_dict['insert_seqs']]
        else:
            raise AttributeError(f'Current barcode details '\
                f'{self.barcode_details_json} does not have an ' \
                    f'insert seq component')
    
    # methods ------------------------------------------------------------------
    def component_edit_distance(self, barcode:str) -> dict:
        """Check the barcode against the expected values at each component 
        substring location.

        Args:
            barcode (str): A barcode string (the entire thing)

        Raises:
            AttributeError: Raised if self.barcode_dict is empty
            IndexError: Raised if the indices of a given component are out 
            of bounds for the barcode input string

        Returns:
            dict: A dictionary where the component names are keys and the values 
            are the minimum edit distance between the barcode component 
            substring and the values listed in the barcode_details json 
        """
        if barcode == "":
            raise AttributeError("No barcode set")

        if not self.barcode_dict:
            raise AttributeError("Barcode Dict is empty. Set a the barcode" +
                "dict by calling BarcodeChecker.set_barcode_details(json_path)")
        # Using the keys of the barcode_components dict, make a dictionary
        # with structure {component1: "TAG", component2:"TTAAGG", ...} where the keys are
        # barcode components, and the values the substring of the barcode 
        # which corresponds
        try:
            parsed_barcode_dict = \
            { k:barcode[v[0]:v[1]] for k,v in self.barcode_dict['indicies'].items()}
        except IndexError as exc:
            raise f'Barcode {barcode} does not have valid indicies for ' \
            f'all expected barcode components. {exc}' from exc

        # iterate over each part, return a dict in form {component: min_edit_dist}
        barcode_summary = {k: self.min_edit_dist(k,v) for k,v in parsed_barcode_dict.items()}

        return barcode_summary

    def barcode_check(self, barcode: str) -> dict:
        """Determine if the barcode passes (True) or fails (False) given the 
        edit distances between it and the expected components, and the allowable 
        edit distance between a given component and the actual value.

        Args:
            component_edit_dist_dict (dict): A dictionary where the keys are 
            barcode components and the values are the minimum edit distance 
            of a given barcode against the corresponding allowable components

        Returns:
            dict: A dict of structure {"pass": Boolean, True if the barcode 
            passes, "tf": Str, where the value is eithe "*" if unknown or a TF 
            string from the barcode_details} 
        """
        if barcode == "":
            raise AttributeError("No barcode set")

        allowance_dict = self.barcode_dict[self.key_dict['match_allowance']]

        tf = self.get_tf(barcode)
        # note that this is redundant -- the same code essentially is used 
        # to check edit distance below.
        # TODO implement a different parser for each component type
        restriction_enzyme = self.get_restriction_enzyme(barcode)
        component_edit_dist_dict = self.component_edit_distance(barcode)
        try:
            # iterate over the component edit dist dict and test if a given
            # barcode component matches. If the match is not exact, but is
            # within the match_allowance, increment mismatch_counter. If a
            # component has distance greater than the match_allowance, or
            # exceeds 'max', mark the barcode as failing
            barcode_passing = True
            mismatch_counter = 0
            it = iter(component_edit_dist_dict.items())
            key_value = next(it, None)
            while key_value and barcode_passing:
                if key_value[1] > allowance_dict[key_value[0]]:
                    barcode_passing = False
                if key_value[1] > 0:
                    mismatch_counter += 1
                if mismatch_counter > allowance_dict[self.key_dict['max_mismatch_key']]:
                    barcode_passing = False
                key_value = next(it, None)
        except KeyError:
            KeyError("A given component was not present in the component match dict")

        return ({"pass": barcode_passing, "tf": tf, "restriction_enzyme": restriction_enzyme})

    def min_edit_dist(self, barcode_component:str, barcode_substr:str) -> int:
        """get the minimum levenshtein edit distance between a given barcode 
        component and the allowable strings for that component

        Args:
            barcode_component (str): Name of a given barcode component
            barcode_substr (str): A substring of the read barcode to check 
            against the the barcode_details['component'] list

        Raises:
            KeyError: Raised if the barcode_component does not exist in 
            self.barcode_dict

        Returns:
            int: The minimum distance between the barcode_substr and the 
            values allowable for the barcode_component
        """
        if barcode_component not in self.barcode_dict[self.key_dict['components']]:
            raise KeyError(f"{barcode_component} is not in " \
                f"the keys of the current barcode_details dict "\
                    f"{self.barcode_details_json}")
        
        component_set = self.barcode_dict[self.key_dict['components']][barcode_component]
        
        # TODO handle hard coding
        if isinstance(component_set, dict):
            min_dist = 1
            for seq in component_set['seq']:
                if seq in barcode_substr:
                    min_dist = 0
                    break
        else:
            min_dist = min([self.edit_dist(barcode_substr, valid_component) \
                for valid_component in component_set])

        return min_dist

    def edit_dist(self, s1:str, s2:str) -> int:
        """Calculate the levenshtein distance between two strings
       
        Cite: https://stackoverflow.com/a/32558749

        Args:
            s1 (str): One of two strings to match against the other
            s2 (str): One of two strings to match against the other

        Returns:
            int: A number where 0 >= n <= longer of the two strings
        """
        if len(s1) > len(s2):
            s1, s2 = s2, s1

        distances = range(len(s1) + 1)
        for i2, c2 in enumerate(s2):
            distances_ = [i2+1]
            for i1, c1 in enumerate(s1):
                if c1 == c2:
                    distances_.append(distances[i1])
                else:
                    distances_.append(1 + min((distances[i1], 
                                      distances[i1 + 1], distances_[-1])))
            distances = distances_
        return distances[-1]
    
    def get_restriction_enzyme(self, barcode:str) -> str:
        """extract the restriction enzyme name, if one exists. Note that this is 
        set up for exact matches only, no correction/fuzzy matching

        Args:
            barcode (str): _description_

        Raises:
            AttributeError: when no barcode is set in self

        Returns:
            str: either the name of the restriction enzyme, or "*"
        """
 
        if barcode == "":
            raise AttributeError("No barcode set")

        seq_indicies = self.barcode_dict[self.key_dict['indicies']][self.key_dict['restriction_enzyme']]
        restriction_seq = barcode[seq_indicies[0]:seq_indicies[1]]
        restriction_enzyme_dict = self.barcode_dict[self.key_dict['components']][self.key_dict['restriction_enzyme']]
        
        # TODO address hard coding in seq
        for i in range(len(restriction_enzyme_dict['seq'])):
            seq = restriction_enzyme_dict['seq'][i]
            if seq in restriction_seq:
                return restriction_enzyme_dict['name'][i]
        
        return "*"

    def get_tf_barcode(self, barcode:str) -> str:
        """If the tf_map attribute is set, use the tf_bc_components to get the 
        string at the correct indicies in the input barcode.

        Args:
            barcode (str): the barcode extracted from a given read

        Raises:
            KeyError: Raised if tf_map does not exist in 
            self.barcode_dict

        Returns:
            str: The substring of barcode at the indicies given in the 
            tf_map.tf_bc_components
        """
        if barcode == "":
            raise AttributeError("No barcode set")

        if self.key_dict['tf_map'] not in self.barcode_dict.keys():
            raise KeyError('tf_map DNE')
        # extract the tf barcode from the appropriate ranges
        tf_bc = []
        barcode_components = self.barcode_dict[self.key_dict['tf_map']]\
            [self.key_dict['tf_map_keys']['bc_components']]
        for component in barcode_components:
            comp_indicies = self.barcode_dict[self.key_dict['indicies']]\
                [component]
            tf_bc.append(barcode[comp_indicies[0]:comp_indicies[1]])
        return "".join(tf_bc)

    def get_tf(self, barcode:str) -> str:
        """given a barcode string, return the tf if it exists in the tf map

        Args:
            barcode (str): the barcode extracted from the read

        Returns:
            str: If the barcode is found in the tf_dict, the tf name will 
            be returned. else "*"
        """
        if barcode == "":
            raise AttributeError("No barcode set")
        
        # if the tf_map in the barcode_details DNE return "*"
        try:
            tf_barcode = self.get_tf_barcode(barcode)
            # if there is an exact match, return the TF
            try:
                tf = self.barcode_dict[self.key_dict['tf_dict']][tf_barcode]
            # but if there isn't, return the best matching tf and the edit 
            # distance, eg MIG2_2 where the _2 is the edit distance and MIG2
            # is the best match
            except KeyError:
                curr_min = 100
                best_match = ""
                for bc in self.barcode_dict[self.key_dict['tf_dict']]:
                    dist = self.edit_dist(tf_barcode,bc)
                    if dist < curr_min:
                        curr_min = dist
                        best_match = self.barcode_dict[self.key_dict['tf_dict']][bc]
                tf = best_match + f"_{curr_min}"
        except KeyError:
            tf = "*"

        return tf

    def tf_dict(self, barcode_dict):
        """Set the tf_dict field in the barcode_dict attribute

        Raises:
            AttributeError: Raised if the number of items in the different 
            TF barcode categories differ
        """
        if self.key_dict['tf_map'] not in barcode_dict.keys():
            raise KeyError(f"{self.key_dict['tf_map']} not in barcode_dict")
        
        # extract the following simply to make the dictionary comprehension
        # easier to read
        barcode_components = barcode_dict[self.key_dict['tf_map']]\
            [self.key_dict['tf_map_keys']['bc_components']]
        components_list = barcode_dict[self.key_dict['components']]
        # extract the lengths of each of the components of the barcode -- 
        # if there are more than one, check to make sure that the lists are 
        # of equal length
        component_lengths = { k:len(components_list[k]) for k in barcode_components }
        expected_barcode_vector_length = component_lengths[list(component_lengths.keys())[0]]
        # check that the number of elements in the given category matches 
        # every other category
        for component,vector_length in component_lengths.items():
            if vector_length != expected_barcode_vector_length:
                raise AttributeError(f'Number of items in barcode ' \
                    f'component {component} components lists do not match component "\
                        f"{list(component_lengths.keys())[0]}')
        # this will be constructed in the loop below. It will have the 
        # structure {'sequence_barcode': TF_name , ... }
        tf_dict = {}
        for i in range(expected_barcode_vector_length):
            tf = barcode_dict[self.key_dict['tf_map']][self.key_dict['tf_map_keys']['tf']][i]
            tf_bc = []
            for component in component_lengths:
                tf_bc.append(barcode_dict[self.key_dict['components']][component][i])
            if tf in tf_dict:
                raise KeyError(f"The same TF appears more than once in ]"\
                    f"{self.key_dict['tf_map']}")
            # once the barcode is constructed, add it to the tf_dict
            tf_dict.update({"".join(tf_bc):tf})

        return tf_dict
