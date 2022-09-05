"""An object which parses a barcode string extracted from a read according 
to a barcode_details json file. Note that all sequences are cast to upper"""
import os
import json

class BarcodeParser:
    """Using a json which describes acceptable values for given barcodes, 
    check the edit distance between the barcode components (substrings) and 
    the """
    # set attributes
    key_dict = {
        "indicies":"indicies",
        "components":"components",
        "insert_seq": "insert_seq",
        "tf_map": "tf_map",
        "tf_map_keys": {"bc_components": "bc_components", "tf": "tf"},
        # this is created from tf_map, which is provided in barcode_details
        "tf_dict": "tf_dict",
        "match_allowance": "match_allowance"
    }
    barcode_dict = {}
    barcode_details_json = ""
    barcode = ""
    # constructor
    def __init__(self, barcode_details_json):
        # set attributes
        self.barcode_details_json = barcode_details_json
        # set the initialized barcode details
        # note that this also sets the match allowances
        self.set_barcode_details(barcode_details_json)

    def set_barcode(self, barcode):
        """setter for the barcode attribute. Note that the barcode is cast to 
        upper

        Args:
            barcode (str): barcode to store in the barcode attribute
        """
        self.barcode = barcode.upper()

    def set_barcode_details(self, barcode_details_json):
        """parse the calling cards bed file name column into components and perform
        This expects a dataframe with AT LEAST the columns (there may be more):
        ['chrom','chromStart','chromEnd','barcode', 'strand', 'insert_seq']

        :param bed_df:
        :param barcode_details:
        :param fltr_bed_output_name:

        """
        # check barcode_details path
        for input_path in [barcode_details_json]:
            if not os.path.exists(input_path):
                raise FileNotFoundError(f"Input file DNE: {input_path}")
        # open json, read in as dict
        with open(barcode_details_json) as f1:
            barcode_dict = json.load(f1)
        # verify format
        barcode_dict = self.verify_barcode_details(barcode_dict)
        # update self
        self.barcode_dict = barcode_dict
        self.set_component_match_allowances()
        self.barcode_details_json = barcode_details_json
    
    def verify_barcode_details(self, barcode_dict):
        """Verify expectations on the barcode details

        Args:
            barcode_dict (dict): A barcode_details dictionary parsed from the 
            barcode_details json file

        Raises:
            KeyError: Raised if an expected key in the barcode_details DNE
            ValueError: Raised if a value in a given field in the barcode_details 
            does not match expectations
        
        Returns:
            (dict): The verified barcode_details dictionary
        """
        # check that the indicies and components match
        if not list(barcode_dict[self.key_dict['components']].keys()).sort() == \
            list(barcode_dict[self.key_dict['indicies']].keys()).sort():
            raise KeyError('The keys in the barcode_component json and the ' +
                'barcode_component_indicies json are not the same')
        
        # check that the components are the correct length according to the 
        # indicies
        barcode_components_length_check_dict = { k:v[1]-v[0] for k,v in \
                                                barcode_dict['indicies'].items()}
        for k,v in barcode_dict['components'].items():
            for comp in v:
                if not len(comp) == barcode_components_length_check_dict[k]:
                    raise ValueError(f"There exists a component in "\
                    f"barcode_components which is not the length " \
                    f"described by the barcode_details['indicies']: "\
                    f"{comp}, {barcode_components_length_check_dict[k]}")
        
        # if a tf_map dict exists, check that the bc_components and tf list 
        # exist
        if self.key_dict['tf_map'] in barcode_dict:
            for bc_comp in barcode_dict[self.key_dict['tf_map']]\
                [self.key_dict['tf_map_keys']['bc_components']]:
                if bc_comp not in barcode_dict[self.key_dict['components']]:
                    raise KeyError(f"tf barcode component {bc_comp} is not " \
                        f"in barcode details {self.key_dict['components']}")
        # if the match_allowance field exists, check that the components match 
        # those in the components field, and that the distance is an integer
        if self.key_dict['match_allowance'] in barcode_dict:
            for comp, dist in barcode_dict[self.key_dict['match_allowance']].items():
                if comp not in barcode_dict[self.key_dict['components']]:
                    raise KeyError(f"The keys in "\
                        f"{self.barcode_dict[self.key_dict['match_allowance']]} " \
                            f"must be in {self.barcode_dict[self.key_dict['components']]}")
                if not isinstance(dist, int) or dist < 0:
                    raise ValueError("Match allowance must be a positive integer which \
                        represents the maximum edit distance from the expected \
                            barcode")
        
        # cast all values in self.key_dict['components'] to upper
        # note that the barcode is also always cast to upper
        upper_components_dict = {}
        for k,seq_list in barcode_dict[self.key_dict['components']].items():
            upper_seq_list = []
            for sequence in seq_list:
                upper_seq_list.append(sequence.upper())
            upper_components_dict.setdefault(k,upper_seq_list)
        barcode_dict[self.key_dict['components']] = upper_components_dict

        # cast insert seq to upper if it exists
        if self.key_dict['insert_seq'] in barcode_dict:
            upper_insert_seq = []
            for sequence in barcode_dict[self.key_dict['insert_seq']]:
                upper_insert_seq.append(sequence.upper())
            barcode_dict[self.key_dict['insert_seq']] = upper_insert_seq
        
        # if tf_map exists, then set the tf_dict
        if self.key_dict['tf_map'] in barcode_dict:
            barcode_dict[self.key_dict['tf_dict']] = self.tf_dict(barcode_dict)
        
        return barcode_dict

    def set_component_match_allowances(self):
        """update the barcode_dict attribute with match allowances for all 
        components which are not already set in the barcode_details json"""
        # if the match_allowance field DNE, create one with an empty dict
        x = self.barcode_dict.setdefault(self.key_dict['match_allowance'], {})
        # for each key in the barcode_dict[comp], set to 0 if it does not 
        # already exist
        for comp in self.barcode_dict[self.key_dict['components']]:
            x = self.barcode_dict[self.key_dict['match_allowance']]\
                .setdefault(comp,0)

    def component_edit_distance(self):
        """Check the barcode against the expected values at each component 
        substring location.

        Args:
            barcode (Str): A barcode string (the entire thing)

        Raises:
            AttributeError: Raised if self.barcode_dict is empty
            IndexError: Raised if the indices of a given component are out 
            of bounds for the barcode input string

        Returns:
            Dict: A dictionary where the component names are keys and the values 
            are the minimum edit distance between the barcode component 
            substring and the values listed in the barcode_details json 
        """
        if self.barcode == "":
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
            { k:self.barcode[v[0]:v[1]] for k,v in self.barcode_dict['indicies'].items()}
        except IndexError as exc:
            raise f'Barcode {self.barcode} does not have valid indicies for ' \
            f'all expected barcode components. {exc}' from exc

        # iterate over each part, return a dict in form {component: min_edit_dist}
        barcode_summary = {k: self.min_edit_dist(k,v) for k,v in parsed_barcode_dict.items()}

        return barcode_summary

    def barcode_check(self):
        """Determine if the barcode passes (True) or fails (False) given the 
        edit distances between it and the expected components, and the allowable 
        edit distance between a given component and the actual value.

        Args:
            component_edit_dist_dict (Dict): A dictionary where the keys are 
            barcode components and the values are the minimum edit distance 
            of a given barcode against the corresponding allowable components

        Returns:
            (Dict): A dict of structure {"pass": Boolean, True if the barcode 
            passes, "tf": Str, where the value is eithe "*" if unknown or a TF 
            string from the barcode_details} 
        """
        if self.barcode == "":
            raise AttributeError("No barcode set")

        tf = self.get_tf()
        component_edit_dist_dict = self.component_edit_distance()
        try:
            false_count = sum([0 if v <= \
                self.barcode_dict[self.key_dict['match_allowance']][k] else 1 \
                for k,v in component_edit_dist_dict.items()])
        except KeyError:
            KeyError("A given component was not present in the component match dict")

        bc_pass = True if false_count == 0 else False

        return ({"pass": bc_pass, "tf": tf})

    def min_edit_dist(self, barcode_component, barcode_substr):
        """get the minimum levenshtein edit distance between a given barcode 
        component and the allowable strings for that component

        Args:
            barcode_component (Str): Name of a given barcode component
            barcode_substr (Str): A substring of the read barcode to check 
            against the the barcode_details['component'] list

        Raises:
            KeyError: Raised if the barcode_component does not exist in 
            self.barcode_dict

        Returns:
            Int: The minimum distance between the barcode_substr and the 
            values allowable for the barcode_component
        """
        if barcode_component not in self.barcode_dict[self.key_dict['components']]:
            raise KeyError(f"{barcode_component} is not in " \
                f"the keys of the current barcode_details dict "\
                    f"{self.barcode_details_json}")

        min_dist = min([self.edit_dist(barcode_substr, valid_component) for \
            valid_component in self.barcode_dict[self.key_dict['components']][barcode_component]])

        return min_dist

    def edit_dist(self, s1, s2):
        """Calculate the levenshtein distance between two strings

        Cite:
            https://stackoverflow.com/a/32558749

        Args:
            s1 (Str): One of two strings to match against the other
            s2 (Str): One of two strings to match against the other

        Returns:
            Int: A number where 0 >= n <= longer of the two strings
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

    def get_insert_seqs(self):
        """Getter for the insert seq sequence from the barcode details json

        Raises:
            AttributeError: Raised if the current barcode details json does 
            not have an insert seq key 

        Returns:
            Str: The insert seq string (upper case) from the barcode details
        """
        if self.key_dict['insert_seq'] in self.barcode_dict:
            return self.barcode_dict[self.key_dict['insert_seq']]
        else:
            raise AttributeError(f'Current barcode details '\
                f'{self.barcode_details_json} does not have an ' \
                    f'insert seq component')

    def get_tf_barcode(self):
        """If the tf_map attribute is set, use the tf_bc_components to get the 
        string at the correct indicies in the input barcode.

        Args:
            barcode (Str): the barcode extracted from a given read

        Raises:
            KeyError: Raised if tf_map does not exist in 
            self.barcode_dict

        Returns:
            Str: The substring of barcode at the indicies given in the 
            tf_map.tf_bc_components
        """
        if self.barcode == "":
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
            tf_bc.append(self.barcode[comp_indicies[0]:comp_indicies[1]])
        return "".join(tf_bc)

    def get_tf(self):
        """given a barcode string, return the tf if it exists in the tf map

        Args:
            barcode (Str): the barcode extracted from the read

        Returns:
            Str: If the barcode is found in the tf_dict, the tf name will 
            be returned. else "*"
        """
        if self.barcode == "":
            raise AttributeError("No barcode set")
        
        # if the tf_map in the barcode_details DNE return "*"
        try:
            tf_barcode = self.get_tf_barcode()
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
                    f'component {k} components lists do not match component "\
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
