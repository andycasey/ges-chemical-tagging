#!/usr/bin/env python

""" Handle the data set for the cluster chemical tagging project """

from __future__ import absolute_import, print_function, with_statement

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"
__all__ = ("DataSet", )

# Third-party.
import numpy as np
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
            cluster_column = Column(data=[""] * len(self.data),
                name="FIELD/CLUSTER", dtype="|S256")
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
        return cls(data)


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
        """
        Assign stars as cluster candidates if they have attributes that meet
        the filtering criteria.

        :param cluster_name:
            The cluster name to assign stars as candidates to.

        :type cluster_name:
            str

        :param star_filter:
            A callable function or evaluable string that describes whether a
            star is a candidate of `cluster_name`. If the function evaluates
            to be True, then the star is set as a cluster candidate.

        :type star_filter:
            callable or str
        """

        if cluster_name.strip().endswith("?"):
            raise ValueError("cluster name cannot end with a question mark")
        return self._assign_cluster("{}?".format(cluster_name), star_filter)


    def assign_cluster_members(self, cluster_name, star_filter):
        """
        Assign stars as cluster members if they have attributes that meet
        the filtering criteria.

        :param cluster_name:
            The cluster name to assign stars as members to.

        :type cluster_name:
            str

        :param star_filter:
            A callable function or evaluable string that describes whether a
            star is a member of `cluster_name`. If the function evaluates
            to be True, then the star is set as a cluster member.

        :type star_filter:
            callable or str
        """

        return self._assign_cluster(cluster_name, star_filter)


    def _assign_cluster(self, cluster_name, star_filter):
        """ Assign stars to clusters. """

        if not hasattr(star_filter, "__call__") \
        and not isinstance(star_filter, (str, unicode)):
            raise TypeError("star filter must be a callable function or an "\
                "evaluable string")

        # Do we need to update the FIELD/CLUSTER dtype?
        _ = self.data.dtype.names.index("FIELD/CLUSTER")
        max_len = int(self.data.dtype[_].str[2:])
        if len(cluster_name) > max_len:
            raise ValueError("cluster name '{0}' is too long; max length {1}"\
                .format(cluster_name, max_len))

        evaluate_filter = not hasattr(star_filter, "__call__")
        mask = np.zeros(len(self.data), dtype=bool)
        for i, row in enumerate(self.data):
            if evaluate_filter:
                try:
                    mask[i] = eval(star_filter, env={"row": row})
                except:
                    logger.exception("Exception evaluating '{0}' on row {1}:"\
                        .format(star_filter, i))
            else:
                mask[i] = star_filter(row)

        # Update, and return the number of affected rows.
        self.data["FIELD/CLUSTER"][mask] = cluster_name
        return int(mask.sum())




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