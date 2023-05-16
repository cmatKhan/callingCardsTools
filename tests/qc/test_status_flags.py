from callingcardstools.QC.StatusFlags import StatusFlags

def test_status_flags():
	sf = StatusFlags

	status = 3

	assert sf.decompose(status, as_str=False) == [0, 1]

	assert sf.decompose(status) == ['BARCODE', 'MAPQ']