from pathlib import PurePosixPath
from typing import Any, Dict

import fsspec
import numpy as np

from kedro.io import AbstractDataset
from kedro.io.core import get_filepath_str, get_protocol_and_path
import logging

import anndata
import mudata

logger = logging.getLogger(__name__)


class AnnDataset(AbstractDataset[anndata.AnnData, anndata.AnnData]):
    def __init__(self, filepath: str):
        """Creates a new instance of GrapeDataset to load / save a grape.Graph for given filepath.
        `
        Args:
            filepath: The location of the grapes file to load / save data.
                      Convention is: filepath is a directory with two files: `nodes.tsv` and `edges.tsv
        """
        # parse the path and protocol (e.g. file, http, s3, etc.) for nodes
        protocol, path = get_protocol_and_path(
            PurePosixPath(filepath)
        )
        self._protocol = protocol
        self._filepath = PurePosixPath(path)
        self._fs = fsspec.filesystem(self._protocol)

    def _load(self) -> np.ndarray:
        """Loads data from the image file.

        Returns:
            Data from the image file as a numpy array
        """
        # using get_filepath_str ensures that the protocol and path are appended correctly for different filesystems
        load_path = get_filepath_str(self._filepath, self._protocol)
        logger.debug(f"Loading anndata from {load_path}!")

        adata = anndata.read_h5ad(load_path)

        logger.debug("Done loading!")
        return adata

    def _save(self, adata: anndata.AnnData) -> None:
        """Saves grape.Graph to the specified filepath."""

        save_path = get_filepath_str(self._filepath, self._protocol)

        logger.debug(f"Saving anndata nodes to {save_path}")

        adata.write_h5ad(save_path)
        logger.debug(f"Saving anndata edges to {save_path}")


    def _describe(self) -> Dict[str, Any]:
        return {"summy": "dummy"}


class MuDataset(AbstractDataset[anndata.AnnData, anndata.AnnData]):
    def __init__(self, filepath: str):
        """Creates a new instance of MuDataset to load / savee for given filepath.
        `
        Args:
            filepath: The location of the grapes file to load / save data.
                      Convention is: filepath is a directory with two files: `nodes.tsv` and `edges.tsv
        """
        # parse the path and protocol (e.g. file, http, s3, etc.) for nodes
        protocol, path = get_protocol_and_path(
            PurePosixPath(filepath)
        )
        self._protocol = protocol
        self._filepath = PurePosixPath(path)
        self._fs = fsspec.filesystem(self._protocol)

    def _load(self) -> np.ndarray:
        """Loads data from the image file.

        Returns:
            Data from the image file as a numpy array
        """
        # using get_filepath_str ensures that the protocol and path are appended correctly for different filesystems
        load_path = get_filepath_str(self._filepath, self._protocol)
        logger.debug(f"Loading MuData from {load_path}!")

        # adata = anndata.read_h5ad(load_path)
        mdata = mudata.read(load_path)
        logger.debug("Done loading!")
        return mdata

    def _save(self, mdata: mudata.MuData) -> None:
        """Saves grape.Graph to the specified filepath."""

        save_path = get_filepath_str(self._filepath, self._protocol)

        logger.debug(f"Saving MuData nodes to {save_path}")

        mdata.write(save_path)
        logger.debug(f"Saving MuData edges to {save_path}")


    def _describe(self) -> Dict[str, Any]:
        return {"summy": "dummy"}
