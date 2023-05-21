
# pylint:disable=W1203,W0102
import logging
from typing import Callable

from pysam import AlignedSegment  # pylint:disable=E0611

from callingcardstools.QC.StatusFlags import StatusFlags

__all__ = ['create_status_coder']

logger = logging.getLogger(__name__)


def create_status_coder(
        insert_seqs: list = ['*'],
        mapq_threshold: int = 10,
        check_5_prime_clip: bool = False,
        check_passing: bool = True) -> Callable[[AlignedSegment], int]:
    """
    A factory function which returns a function capable of determining the
    status code of a read tagged by an AlignmentTagger object.

    Args:
        insert_seqs (list): A list of acceptable insert sequences. Defaults
            to ['*'], which will skip the insert seq check altogether.
        mapq_threshold (int): A mapq_threshold. Reads with map quality less
            than this value will be marked as failing the mapq threshold test.
            Default is 10.
        check_5_prime_clip (bool): Whether to check for 5' end clipping in
            the read. Defaults to False.
        check_passing (bool, optional): Whether to check the passing key in
            the barcode_details dict. Defaults to True.

    Returns:
        Callable[[AlignedSegment], int]: A function which given a tagged
        pysam AlignedSegment will return the status code for the read.
    """

    def coder(read_details: AlignedSegment,
              status_code: int = 0) -> int:
        """
        Returns the status code for a given read after checking for various 
        flags.

        Args:
            read_details (AlignedSegment): A pysam AlignedSegment object.
            status_code (int, optional): Initial status code. Defaults to 0.

        Raises:
            ValueError: If read_details is not a dictionary or does not
                contain expected keys.
            KeyError: If required keys are not present in read_details.
            ValueError: If the types of values in read_details do not 
                match the expected types.

        Returns:
            int: The status code for a given read.
        """
        if not isinstance(read_details, dict):
            raise ValueError('read_details must be a dictionary with '
                             'keys "read" which is a pysam.AlignedSegment and '
                             '"barcode_details" which is a dict')
        if not {'read', 'barcode_details'} - read_details.keys() == set():
            raise KeyError('"read" and "barcode_details" must be keys in'
                           'read_details')
        if not isinstance(read_details.get('read'), AlignedSegment):
            raise ValueError('read_details["read"] must be a '
                             'pysam.AlignedSegment object')
        if not isinstance(read_details.get('barcode_details'), dict):
            raise ValueError('read_details["barcode_details"] must be a '
                             'dict')
        if check_passing:
            if not isinstance(
                    read_details.get('barcode_details').get('passing', None),
                    bool):
                raise KeyError('passing must be a key in '
                               'read_details["barcode_details"]')

        # if check passing is set to false, then the passing key may not
        # exist. In this event, assume the read is passing
        if not read_details.get('barcode_details').get('passing', True):
            status_code += StatusFlags.BARCODE.flag()
        # if the read is unmapped, add the flag, but don't check
        # other alignment metrics
        if read_details.get('read').is_unmapped:
            status_code += StatusFlags.UNMAPPED.flag()
        else:
            if read_details.get('read').is_qcfail:
                status_code += StatusFlags.ALIGNER_QC_FAIL.flag()
            if read_details.get('read').is_secondary or \
                    read_details.get('read').is_supplementary:
                status_code += StatusFlags.NOT_PRIMARY.flag()
            if read_details.get('read').mapping_quality < mapq_threshold:
                status_code += StatusFlags.MAPQ.flag()
            # note: for mammals, this isn't necessary as the insert seq can
            # be checked
            if check_5_prime_clip:
                # if the read is clipped on the 5' end, flag
                if (read_details.get('read').is_forward and
                    read_details.get('read').query_alignment_start != 0) or \
                        (read_details.get('read').is_reverse and
                            read_details.get('read').query_alignment_end !=
                         read_details.get('read').infer_query_length()):
                    status_code += StatusFlags.FIVE_PRIME_CLIP.flag()
        # check the insert sequence
        try:
            if insert_seqs != ["*"]:
                if read_details.get('read').get_tag("XZ") not in insert_seqs:
                    status_code += StatusFlags.INSERT_SEQ.flag()
        except AttributeError as exc:
            logger.debug(
                f"insert sequence not found in Barcode Parser. {exc}")

        return status_code

    return coder
