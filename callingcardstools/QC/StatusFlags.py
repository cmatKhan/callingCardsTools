"""Enumerate Status bit flags which will be used to mark why a read
 passes/fails"""
from typing import Union, List
from enum import IntFlag

from numpy import log2, ndarray
from pandas import Series

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
    def decompose(nums: Union[int, List[int], ndarray, Series],
                  as_str: bool = True) -> List:
        """
        Decompose a number, list, ndarray, or pandas Series of numbers
            representing the sum of the powers of two into a list of those
            powers of two. Optionally, return the string representation of
            the powers of two according to the StatusFlags object.

        Args:
            nums (Union[int, List[int], ndarray, pd.Series]):
                The input number, list, ndarray, or pandas Series to decompose.
            as_str (bool, optional): Whether to return the string
                representation of the powers of two according to the
                StatusFlags object. Defaults to True.

        Returns:
            List: A list representing the sum of the powers of two if
                `as_str` is False, e.g., 10 decomposes into [2, 8].
                If `as_str` is true, then the result would be
                ['MAPQ', 'RESTRICTION_ENZYME'].

        Raises:
            TypeError: If the input type is neither int, list, numpy array,
                nor pandas Series.
            ValueError: If the input is a negative integer.
        """
        def decompose_single(num: int, as_str: bool = True) -> list:
            # check input
            if num < 0:
                raise ValueError("Invalid input, expected positive int")
            # if num is 0 and as_str is true, return NO_STATUS. otherwise, the
            # decomposed list will be empty
            if num == 0:
                if as_str:
                    return ['NO_STATUS']
            # Use list comprehension to find the powers of two that
            # compose the input number
            powers = [int(log2(1 << i)) for i
                      in range(num.bit_length()) if num & (1 << i)]

            if as_str:
                # Convert the powers of two to their string representation
                powers = [StatusFlags(int(x)).name for x in powers]

            return powers

        if isinstance(nums, (list, ndarray, Series)):
            return [decompose_single(num, as_str) for num in nums]
        elif isinstance(nums, int):
            return decompose_single(nums, as_str)
        else:
            raise TypeError(
                "Invalid input type, expected int, list, or numpy array")

    # def decompose(num: int, as_str: bool = True) -> list:
    #     """decompose a number which represents the sum of the powers of two
    #         into a list of those powers of two. Optionally,

    #     Args:
    #         num (int): the flag to decompose, eg 10
    #         as_str (bool, optional): whether to return the string
    #             representation of the powers of two according to the
    #             StatusFlag object. Defaults to True.

    #     Returns:
    #         list: list representing the sum of the powers of two if `as_str` is
    #          False, eg 10 decomposes into [2,8]. If `as_str` is true, then the
    #             result would be ['MAPQ', 'RESTRICTION_ENZYME']
    #     """
    #     # cite: https://codereview.stackexchange.com/a/201461
    #     powers = []
    #     while num != 0:
    #         powers.append(log2(num & -num))
    #         num = num & (num - 1)

    #     if as_str:
    #         powers = [StatusFlags(int(x)).name for x in powers]

    #     return powers
