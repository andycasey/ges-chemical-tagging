#!/usr/bin/env python

""" Create realisations of cluster and field stars for a MCMC analysis """

from __future__ import absolute_import, print_function, with_statement

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"


import numpy as np
import pyfits


def create_realisation(data, num_clusters=np.inf, exclude_clusters=None,
    cluster_size_limits=None, field_star_fraction=0, 
    cluster_star_constraints=None, field_star_constraints=None):
    """
    Return a set of row indices that mimic the input constraints. This function
    allows the user to set constraints on what the realisation should look like.

    Possible permutations might include:

    # cluster names to exclude from realisation
    # number of clusters to include in simulation
    # minimum/maximum number of stars in each cluster
    # fraction of field stars
    # constraints on field stars
    # constraints on cluster stars

    :param data:
        The observed data table containing field and cluster stars.

    :type data:
        :class:`pyfits.hdu.table.BinTableHDU`

    :param num_clusters: [optional]
        The number of clusters to draw from. By default the realisation will
        draw from all clusters.

    :type num_clusters:
        int

    :param exclude_clusters: [optional]
        A list of cluster names to exclude in this realisation.

    :type exclude_clusters:
        tuple

    :param cluster_size_limits: [optional]
        A two-length tuple providing the lower and upper limits on how many
        stars to draw from any given cluster. If no limit is provided then the
        realisation can be dominated by whichever clusters have been observed
        more frequently. To place a single bound use (None, x) or (x, None) to
        set upper or lower limits, respectively.

    :type cluster_size_limits:
        tuple

    :param field_star_fraction: [optional]
        Specify the fraction of field stars for this realisation. This value
        must be between [0, 1]. 

    :type field_star_fraction:
        float

    :param cluster_star_constraints: [optional]
        A callable function that takes a single argument (the data row) and 
        returns true or false as to whether a cluster star should be employed 
        for this realisation or not.

    :type cluster_star_constraints:
        callable

    :param field_star_constraints: [optional]
        A callable function that takes a single argument (the data row) and
        returns true or false as to whether a field star should be employed
        for this realisation or not.

    :type field_star_constraints:
        callable

    :returns:
        Data indices for stars in this realisation, and a dictionary containing
        the names of the clusters employed (keys) and the number of stars from
        that cluster (values).

    :raises:
        A shitload of exceptions for different circumstances that I can't even
        fathom to consider yet.

    """

    # number of clusters is +infinity by default (e.g., all)
    if not np.isposinf(num_clusters):
        try:
            num_clusters = int(num_clusters)
        except (TypeError, ValueError):
            raise TypeError("number of clusters must be an integer")

        if 0 > num_clusters:
            raise ValueError("number of clusters cannot be less than zero")

    if exclude_clusters is None:
        exclude_clusters = ()
    if not isinstance(exclude_clusters, (list, tuple, np.ndarray)):
        raise TypeError("exclude clusters must be a list-like object of strings")
    # Note: If exclude_clusters contains clusters that we don't recognise, just 
    #       log a warning later on

    if cluster_size_limits is not None:
        assert len(cluster_size_limits) == 2, "Cluster size limits must be a "
            "two-length tuple"

    try:
        field_star_fraction = float(field_star_fraction)
    except (TypeError, ValueError):
        raise TypeError("field star fraction must be a float between [0, 1]")

    if not (1 >= field_star_fraction >= 0):
        raise ValueError("field star fraction must be between [0, 1]")

    if cluster_star_constraints is None:
        cluster_star_constraints = lambda _: True
    else:
        if not hasattr(cluster_star_constraints, "__call__"):
            raise TypeError("cluster star constraints must be a callable")

    if field_star_constraints is None:
        field_star_constraints = lambda _: True
    else:
        # Must be callable
        if not hasattr(field_star_constraints, "__call__"):
            raise TypeError("field star constraints must be a callable")

    # OK, now we need some logic.
    # First let's exclude stars based on our criteria:
    # - Exclude any stars based on which cluster they are in.
    # - Exclude any cluster stars based on our constraint functions
    # - Exclude any field stars based on our constraint function

    # We can presumably always add more field stars because we have more of them
    # so we should do the cluster star selections first (unless the field_star)

    # If we have a cluster but the number of stars in that cluster is below the
    # minimum limit then we should raise an exception

    # If we have been asked for num_clusters to be more than the actual clusters
    # we have (and num_clusters is not +infinity) then we should raise an
    # exception

    # If we cannot draw non-duplicating field stars (e.g., there are not enough
    # field stars to obtain the field fraction) then we should raise an
    # exception.

    # Steps:
    # (1) Randomly pick `num_cluster` names
    # (2) In each cluster, pick a number between cluster_size_limits, or if none
    #     are given, somewhere between (0, number_of_stars).
    # (3) Randomly select non-duplicating rows in each cluster to identify the
    #     actual members.
    # (4) If we have field_star_fraction > 0, we should calculate the number of
    #     field stars we need to include in this sample, based on how many
    #     cluster stars we already have.
    #     
    #     Randomly draw non-duplicating field stars until we have the required
    #     fraction.
    # (5) Return the row indices and the counter containing the number of stars
    #     in each cluster.

    raise NotImplementedError


data_indices, cluster_sizes = create_realisation(data, exclude_clusters=("NGC 2419", ),
    cluster_size_limits=(10, None), field_star_fraction=0.9,
    field_star_constraints=None,
    cluster_star_constraints=None)



