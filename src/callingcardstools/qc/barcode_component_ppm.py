def barcode_components_ppm(barcode_vector_dict):
    """Convert each barcode component into a position probability matrix.
    Writ the ppm for each component to file.

    :param barcode_vector_dict: a dictionary where the keys are the names of the
    components of the barcode, and the values are lists of those components
    which have been parsed out of the corresponding barcode
    """
    freq_table_dict = {k:frequency_matrix(v) for k,v in barcode_vector_dict.items()}

    for k,v in freq_table_dict.items():
        v.to_csv(k+"_ppm.tsv", sep="\t")
