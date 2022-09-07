"""Enumerate Status bit flags which will be used to mark why a read passes/fails"""
from enum import IntFlag

class StatusFlags(IntFlag):
    """A barcode failure is 0x1, a mapq failure is 0x2 and a insert seq failure 
    is 0x3. A read that fails both barcode and mapq for instance would have 
    status 3.
    """
    BARCODE    = 0x0
    MAPQ       = 0x1
    INSERT_SEQ = 0x2

    def flag(self):
        return 2**self.value