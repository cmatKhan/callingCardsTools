from callingcardstools.QC.StatusFlags import StatusFlags

def test_status_flags():
	sf = StatusFlags

	status = 3

	decomposed = sf.decompose(status)

	decomposed_str = [sf(x).name for x in decomposed]

	assert decomposed == [0, 1]

	assert decomposed_str == ['BARCODE', 'MAPQ']