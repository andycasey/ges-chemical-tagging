#!/usr/bin/env python

""" Test the realisation creation """

from __future__ import absolute_import, print_function, with_statement

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

# Standard library.
import logging
import unittest

# Module specific.
import code as tagging

logger = logging.getLogger("code")

DATA_PATH = "data/GES_iDR2iDR3_WG10+WG11.fits"

class TestRealisationCreation(unittest.TestCase):

    def setUp(self):
        """ Set up the data for creating realisations. """

        # Load the data
        dataset = tagging.data.DataSet.from_fits(DATA_PATH, extension=1)

        # Assign all as field.
        dataset.data["FIELD/CLUSTER"] = "FIELD"

        # [TODO] Delete benchmarks
        clusters = ("Cha_I", "Br81", "M15", "NGC2808", "NGC6633", "IC4665", 
            "NGC104", "gamma2_Vel", "GJ880", "NGC4815", "NGC2547", "NGC5927",
            "NGC4833", "NGC1851", "NGC2243", "NGC3532", "NGC6752", "Br25", 
            "NGC4372", "NGC6705", "M67", "NGC2516", "Trumpler20")

        # Assign all as members.
        for cluster in clusters:
            members = dataset.assign_cluster_members(cluster,
                lambda row: row["TARGET"].startswith(cluster))

            # Special hack:
            if cluster == "Trumpler20":
                members += dataset.assign_cluster_members(cluster,
                    lambda row: row["TARGET"].startswith("Trumpler_20"))

        logger.info("Assigned stars to {} clusters".format(len(clusters)))
        self.dataset = dataset
        return None


    def test_something(self):
        print("done")
        return True


    def tearDown(self):

        # Close the image?
        # [TODO]

        return None


