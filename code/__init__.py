#!/usr/bin/env python

""" Chemically tag the Gaia-ESO Survey data. """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

# Standard library.
import logging

# Module-specific.
from . import (data, infer, realisations, plot)
from .data import DataSet

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

