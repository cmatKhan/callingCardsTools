"""Enumerate Status bit flags which will be used to mark why a read
 passes/fails"""
from typing import Union, List
from enum import IntFlag

from numpy import log2, ndarray
from pandas import Series

__all__ = ['StatusFlags']


class StatusFlags(IntFlag):
    """
    A class used to represent different status flags for read alignment.
    Each status flag corresponds to a different power of 2, allowing
    combinations of flags to be represented as a sum of these powers.

    Attributes:
        BARCODE (int): Corresponds to a barcode failure.
        MAPQ (int): Corresponds to a MAPQ failure.
        INSERT_SEQ (int): Corresponds to an insert sequence failure.
        FIVE_PRIME_CLIP (int): Corresponds to a failure due to 5' end
            clipping in the read.
        UNMAPPED (int): Corresponds to the read being unmapped.
        NOT_PRIMARY (int): Corresponds to the read not being primary.
        ALIGNER_QC_FAIL (int): Corresponds to the read failing aligner QC.
        RESTRICTION_ENZYME (int): Corresponds to a failure due to
            restriction enzyme.
    """
    BARCODE = 0x0
    MAPQ = 0x1
    INSERT_SEQ = 0x2
    FIVE_PRIME_CLIP = 0x3
    UNMAPPED = 0x4
    NOT_PRIMARY = 0x5
    ALIGNER_QC_FAIL = 0x6
    RESTRICTION_ENZYME = 0x7

    def flag(self) -> int:
        """
        Returns the power of 2 corresponding to the flag's value.

        Returns:
            int: The power of 2 corresponding to the flag's value.
        """
        return 2**self.value

    @staticmethod
    def decompose(nums: Union[int, List[int], ndarray, Series],
                  as_str: bool = True) -> List[int]:
        """
        Decomposes a number, list, ndarray, or pandas Series of numbers 
        representing the sum of the powers of two into a list of those 
        powers of two. Optionally, return the string representation of the 
        powers of two according to the StatusFlags object.

        Args:
            nums (Union[int, List[int], ndarray, pd.Series]): The input 
                number, list, ndarray, or pandas Series to decompose.
            as_str (bool, optional): Whether to return the string 
                representation of the powers of two according to the 
                StatusFlags object. Defaults to True.

        Returns:
            List: A list representing the sum of the powers of two if 
                `as_str` is False, e.g., 10 decomposes into [2, 8]. If 
                `as_str` is true, then the result would be 
                ['MAPQ', 'RESTRICTION_ENZYME'].

        Raises:
            TypeError: If the input type is neither int, list, numpy array, 
                nor pandas Series.
            ValueError: If the input is a negative integer.
        """
        def decompose_single(num: int, as_str: bool = True) -> list:
            """
            Helper function to decompose a single integer into powers of 2.

            Args:
                num (int): The integer to decompose.
                as_str (bool): Whether to return the string representation of 
                    the powers of two according to the StatusFlags object. 
                    Defaults to True.

            Returns:
                list: List of powers of two composing the input number. If 
                    `as_str` is True, the powers of two are represented by 
                    their corresponding flag names.

            Raises:
                ValueError: If the input integer is negative.
            """
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
