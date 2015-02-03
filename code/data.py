#!/usr/bin/env python

""" Handle the data set for the cluster chemical tagging project """

from __future__ import absolute_import, print_function, with_statement

__all__ = ("DataSet", )
__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

# Third-party
from astropy.io import fits
from astropy.table import Column, Table

class DataSet(object):

    def __init__(self, data):
        """
        Create a DataSet class using the data table provided.

        :param data:
            The data table to use to create the set.

        """

        self.data = data
        self._add_membership_column()

    def _add_membership_column(self):
        """
        Add a FIELD/CLUSTER membership column to the data. We need this for
        creating realisations and evaluating how correct we were later on.
        """
        if "FIELD/CLUSTER" not in self.data.dtype.names:
            cluster_column = Column(
                data=[""] * len(self.data), name="FIELD/CLUSTER", dtype="|S256")
            self.data.add_column(cluster_column)
        return True

    @classmethod
    def from_fits(cls, filename, **kwargs):
        """
        Create a DataSet class from a FITS filename.

        :param filename:
            The FITS filename to load the data set from.

        :type filename:
            str
        """

        extension = kwargs.pop("extension", 0)
        data = Table(fits.open(filename, **kwargs)[extension].data)
        return cls.__init__(data)


    def writeto(self, filename, **kwargs):
        """
        Write the `DataSet` to a new file.

        :param filename:
            The file path.

        :type filename:
            str
        """

        return self.data.write(filename, **kwargs)


    def assign_cluster_candidates(self, cluster_name, star_filter):

        pass

    def assign_cluster_members(self, cluster_name, star_filter):
        pass


    def assign_field_star(self, star_filter):
        pass



# PREPARATION
# Create some data set.
# Assign field stars based on rules.
# Assign cluster candidates based on rules.
# Assign cluster members based on rules.
# Write the data to disk.

# ASSIGNING REALISATIONS
# Create realisations using the data table.
# Save the information about the realisations to disk.

# DATA ANALYSIS
# For each realisation, we need to use particular amounts of data and try and
# work out the number of clusters present.

# Scripts for running on each realisation. Each script might be using a
# different amount of data (e.g., metallicities only?)