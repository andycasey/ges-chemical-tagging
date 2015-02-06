#!/usr/bin/env python

""" Create realisations of cluster and field stars """

from __future__ import absolute_import, print_function, with_statement

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"
__all__ = ("create", "perturb_abundances")

# Standard library
import logging
import random
import warnings
from hashlib import md5

# Third-party
import numpy as np
from astropy.table import Table

# Initialise logger.
logger = logging.getLogger(__name__)

# Reproducibility is key.
#random.seed(888)

# Warn the user if they asked for more clusters than what was available.
class InsufficientClustersAvailableWarning(Warning):
    pass
warnings.simplefilter("once", InsufficientClustersAvailableWarning)

def create(dataset, num_clusters=None, exclude_clusters=None,
    cluster_size_limits=None, field_star_fraction=0, 
    cluster_star_constraints=None, field_star_constraints=None, **kwargs):
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

    :param dataset:
        The observed dataset containing field and cluster stars.

    :type dataset:
        :class:`data.DataSet`

    :param num_clusters: [optional]
        The number of clusters to draw from. If `None` is provided (default) the
        realisation will pick a random (uniform) number of clusters up to the 
        limit of what is available. Use `np.inf` to draw from all clusters. If
        the number of clusters requested exceeds that available, the behaviour
        will be the same as using `np.inf` except a warning will also be made.

    :type num_clusters:
        None or int

    :param exclude_clusters: [optional]
        A list of cluster names to exclude in this realisation.

    :type exclude_clusters:
        tuple of str

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
    
    if not hasattr(dataset, "data") or not isinstance(dataset.data, Table):
        raise TypeError("data table must be a data.DataSet object")

    if num_clusters is not None and not np.isposinf(num_clusters):
        try:
            num_clusters = int(num_clusters)
        except (TypeError, ValueError):
            raise TypeError("number of clusters must be an integer")

        if 0 > num_clusters:
            raise ValueError("number of clusters cannot be less than zero")

    if exclude_clusters is None:
        exclude_clusters = ()
    if not isinstance(exclude_clusters, (list, tuple, np.ndarray)):
        raise TypeError("exclude clusters must be a list of strings")

    # Note: If exclude_clusters contains clusters that we don't recognise, just 
    #       log a warning later on. But we can't have just an empty string.
    if "" in exclude_clusters:
        raise ValueError("cannot exclude cluster names with an empty string")

    if "FIELD" in map(str.upper, exclude_clusters):
        raise ValueError("cannot exclude protected cluster name 'FIELD'")

    if cluster_size_limits is not None:
        assert len(cluster_size_limits) == 2, "Cluster size limits must be a "\
            "two-length tuple"

        minimum, maximum = cluster_size_limits
        if minimum is not None and maximum is not None and minimum > maximum:
            raise ValueError("maximum cluster size must be greater than the "
                "minimum cluster size")
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

    """
    OK, now we need some logic.

    We can presumably always add more field stars because we have more of them
    so we should do the cluster star selections first. If this assumption is
    not true then we will raise an exception later on

    If we have a cluster but the number of stars in that cluster is below the
    minimum limit then we should raise an exception

    If we have been asked for num_clusters to be more than the actual clusters
    we have (and num_clusters is not +infinity) then we should raise an
    exception

    If we cannot draw non-duplicating field stars (e.g., there are not enough
    field stars to obtain the field fraction) then we should raise an
    exception.

    Steps
    -----
    First let's exclude stars based on our criteria:
    (1) Exclude any stars based on which cluster they are in.
    (2) Exclude any cluster stars based on our constraint functions
    (3) Exclude any field stars based on our constraint function
    (4) Randomly pick `num_clusters` names
    (5) In each cluster, pick a number between cluster_size_limits, or if none
        are given, somewhere between (0, number_of_stars).
    (6) Randomly select non-duplicating rows in each cluster to identify the
        actual members.
    (7) If we have field_star_fraction > 0, we should calculate the number of
        field stars we need to include in this sample, based on how many
        cluster stars we already have.
        
        Randomly draw non-duplicating field stars until we have the required
        fraction.
    (8) Return the row indices and the counter containing the number of stars
        in each cluster.

    Notes
    -----
    Data are assumed to have some column named 'CLUSTER_NAME'. If that column
    is None or "" then the star is assumed to be a field star. Stars that are
    actually part of a cluster (after some radial velocity cut) are labelled
    as 'MEMBER_<CLUSTER_NAME>', but those that were targetted (and no velocity
    cut was made) are labelled as 'CANDIDATE_<CLUSTER_NAME>'
    """

    data = dataset.data
    FIELD_IDENTIFIER = kwargs.pop("__FIELD_IDENTIFIER", "FIELD")
    cluster_names = [each for each in \
        set(data["FIELD/CLUSTER"]).difference([FIELD_IDENTIFIER]) \
        if not each.endswith("?")]
    logging.debug("There are {0} clusters with confirmed members: {1}".format(
        len(cluster_names), ", ".join(cluster_names)))

    # Exclusions
    # (1) Exclude any stars based on what cluster they are in.
    # [TODO] Consider checking the clusters for their sizes and automatically
    # excluding those that do not meet the `cluster_size_limits` minimum
    included = np.ones(len(data), dtype=bool)
    for cluster_name in exclude_clusters:
        indices = (data["FIELD/CLUSTER"] == cluster_name)
        included[indices] *= False
        if indices.sum() == 0:
            logging.warn("Found no members or candidates of cluster {0} to "
                "exclude".format(cluster_name))
        else:
            logging.debug("Excluded {0} stars from cluster {1}".format(
                indices.sum(), cluster_name))

    # (2) Exclude any cluster stars based on our constraint function
    # (3) Exclude any field stars based on our constraint function
    for index, row in enumerate(data):
        constraint_function = field_star_constraints if \
            row["FIELD/CLUSTER"] == FIELD_IDENTIFIER else \
            cluster_star_constraints
        included[index] *= constraint_function(row)

    # (4) Randomly pick `num_clusters` names
    included_cluster_names = set(cluster_names).difference(exclude_clusters)
    if num_clusters is None: # Pro-Tip: np.isposinf(None) raises TypeError :)
        # This forces at least one cluster to be drawn from the sample.
        num_clusters = \
            int(round(random.uniform(1, len(included_cluster_names) + 1)))

    elif np.isposinf(num_clusters):
        num_clusters = len(included_cluster_names)

    elif num_clusters > len(included_cluster_names):
        warnings.warn("Number of requested clusters exceeds the available "
            "clusters after accounting for exclusions ({0} > {1}). Setting the "
            "number of clusters to be {1}".format(
                num_clusters, len(included_cluster_names)),
            InsufficientClustersAvailableWarning)
        num_clusters = len(included_cluster_names)

    # Randomly pick `num_clusters` from the subset
    included_cluster_names = random.sample(included_cluster_names,
        num_clusters)
    logging.debug("Chose {0} clusters: {1}".format(num_clusters,
        ", ".join(included_cluster_names)))
    
    # (5) In each cluster, pick a number between cluster_size_limits, or if none
    #     are given, somewhere between (0, number_of_stars).
    cluster_sizes = {}
    chosen_cluster_indices = []
    for cluster_name in included_cluster_names:
        indices = data["FIELD/CLUSTER"] == cluster_name
        num_members = indices.sum()
        default_min, default_max = 0, num_members
        if cluster_size_limits is not None:
            # Check the minimum value
            minimum, maximum = cluster_size_limits
            if minimum is not None and num_members < minimum:
                raise OverConstrainedRealisationException(
                    "the number of members in {0} does not meet the minimum "
                    "cluster size: {1} < {2}".format(cluster_name, num_members,
                        minimum))

            minimum = [minimum, default_min][minimum is None]
            maximum = [maximum, default_max][maximum is None]
            limits = (minimum, maximum)

        else:
            limits = (default_min, default_max)

        # Pick a random number of members for this cluster
        num_chosen_members = np.clip(random.randint(*limits), 0, num_members)
        cluster_sizes[cluster_name] = num_chosen_members

        # (6) Randomly select non-duplicating rows in each cluster to identify 
        #     the actual members.
        chosen_indices = random.sample(np.where(indices)[0], num_chosen_members)
        chosen_cluster_indices.extend(chosen_indices)

    # (7) If we have field_star_fraction > 0, we should calculate the number of
    #     field stars we need to include in this sample, based on how many
    #     cluster stars we already have.
    if field_star_fraction > 0:
        num_cluster_stars = len(chosen_cluster_indices)
        num_field_stars = int(round((field_star_fraction * num_cluster_stars) \
            / (1. - field_star_fraction)))
        logging.debug("Calculated number of field stars to be {0:.0f} to give a "
            "field star fraction of {1:.2f} with {2:.0f} cluster stars".format(
                num_field_stars, field_star_fraction, num_cluster_stars))

        indices = np.where(data["FIELD/CLUSTER"] == FIELD_IDENTIFIER)[0]
        if len(indices) < num_field_stars:
            raise OverConstrainedRealisationException("There are not enough "
                "field stars to yield the required field star fraction. For "
                "{0:.0f} cluster stars, {1:.0f} field stars are required to "
                "yield the requested field star fraction of {2:.2f} but only"
                " {3:.0f} field stars are available.".format(num_cluster_stars,
                    num_field_stars, field_star_fraction, len(indices)))

        if len(indices) == num_field_stars:
            logging.warn("The number of available field stars matches the number"
                "required to yield a field star fraction of {0:.2f}".format(
                    field_star_fraction))

        chosen_field_indices = random.sample(indices, num_field_stars)

    else:
        chosen_field_indices = []

    # (8) Return the row indices and the counter containing the number of stars
    #     in each cluster.
    chosen_indices = sorted(chosen_cluster_indices + chosen_field_indices)

    # Some final assertions to make sure things are as they should be.
    assert len(chosen_cluster_indices) == sum(cluster_sizes.values())
    assert len(chosen_cluster_indices) == len(set(chosen_cluster_indices))
    assert len(chosen_field_indices) == len(set(chosen_field_indices))
    assert len(chosen_indices) == len(set(chosen_indices))

    # Create some checksum of the chosen indices
    checksum = md5(" ".join(map(str, chosen_indices))).hexdigest()
    return (chosen_indices, cluster_sizes, checksum)


class OverConstrainedRealisationException(Exception):
    # This exception class is for when the inputs given to a realisation have
    # made it over constrained. 
    pass

def perturb_abundances(data_subset, abundance_columns):
    raise NotImplementedError

