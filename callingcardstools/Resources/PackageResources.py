import importlib.resources as pkg_resources
from . import yeast, mouse

__all__ = ['Resources']


class Resources:
    """An object to provide access to package resources to the user"""

    def __init__(self) -> None:

        self._configured_organisms = ['yeast', 'mouse']

        self._yeast_resources = {
            'barcode_details':
                pkg_resources.path(yeast, "barcode_details.json"),
            'yeast_db':
                pkg_resources.path(yeast, "yeast.db"),
        }

        self._mouse_resources = {
            'barcode_details':
                pkg_resources.path(mouse, "barcode_details.json"),
        }

    @property
    def configured_organisms(self):
        """list of organisms for which there are resources"""
        return self._configured_organisms

    @property
    def yeast_resources(self):
        """dict of paths to resources for yeast"""
        return self._yeast_resources

    @property
    def mouse_resources(self):
        """dict of paths to resources for mouse"""
        return self._mouse_resources
