#!/usr/bin/env python

""" Chemically tag the Gaia-ESO Survey data. """

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

# Standard library.
import logging

# Module-specific.
from . import (data, infer, realisations)

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

#formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
#ch.setFormatter(formatter)
#root.addHandler(ch)

