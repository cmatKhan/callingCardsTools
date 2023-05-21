import importlib.resources as pkg_resources
from . import yeast, mammals

__all__ = ['Resources']


class Resources:
    """An object to provide access to package resources to the user"""

    def __init__(self) -> None:

        self._configured_organisms = ['yeast', 'mammals']

        self._yeast = {
            'barcode_details':
                pkg_resources.read_text(yeast, "barcode_details.json")
        }

        self._mammals = {
            'barcode_details':
                pkg_resources.read_text(mammals, "barcode_details.json"),
        }

    @property
    def configured_organisms(self):
        """list of organisms for which there are resources"""
        return self._configured_organisms

    @property
    def yeast(self):
        """dict of paths to resources for yeast"""
        return self._yeast

    @property
    def mammals(self):
        """dict of paths to resources for mammals"""
        return self._mammals
