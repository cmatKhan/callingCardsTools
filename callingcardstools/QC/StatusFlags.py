"""Enumerate Status bit flags which will be used to mark why a read
 passes/fails"""
from enum import IntFlag
from numpy import log2

__all__ = ['StatusFlags']


class StatusFlags(IntFlag):
    """A barcode failure is 0x0, a mapq failure is 0x1 and a insert seq failure
     is 0x2. A read that fails both barcode and mapq for instance would have
     status 3.
    """
    BARCODE = 0x0
    MAPQ = 0x1
    INSERT_SEQ = 0x2
    FIVE_PRIME_CLIP = 0x3
    UNMAPPED = 0x4
    NOT_PRIMARY = 0x5
    ALIGNER_QC_FAIL = 0x6
    RESTRICTION_ENZYME = 0x7

    def flag(self):
        return 2**self.value

    @staticmethod
    def decompose(num: int) -> list:
        """decompose a flag into the sum of powers of two

        Args:
            num (int): the flag to decompose, eg 10

        Returns:
            list: list representing the sum of the powers of two, eg 
            10 decomposes into [2,8]
        """
        # cite: https://codereview.stackexchange.com/a/201461
        powers = []
        while num != 0:
            powers.append(log2(num & -num))
            num = num & (num - 1)
        return powers
