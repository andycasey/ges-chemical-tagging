#!/usr/bin/env python

""" Create realisations for tests. """

from __future__ import absolute_import, print_function, with_statement

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

# Standard library.
import cPickle as pickle
import logging
import multiprocessing as mp

# Third-party.
import numpy as np

# Module-specific.
import code as tagging

logger = logging.getLogger("code")

# Load dataset.
dataset = tagging.DataSet.from_fits("data/GES_iDR2iDR3_WG10+WG11.fits",
    extension=1)

# Assign field/cluster stars.
# Unless otherwise told, this star is in the field:
dataset.data["FIELD/CLUSTER"] = "FIELD"

# Here we will define the cluster names which will be used to describe the
# candidate listing, and we will use specific rules to mark them as true
# candidates or not.
clusters = {
    "Cha_I": lambda row: row["FIELD/CLUSTER"] == "Cha_I?" and True,
    "Br81": lambda row: row["FIELD/CLUSTER"] == "Br81?" and True,
    "M15": lambda row: row["FIELD/CLUSTER"] == "M15?" and True,
    "NGC2808": lambda row: row["FIELD/CLUSTER"] == "NGC2808?" and True,
    "NGC6633": lambda row: row["FIELD/CLUSTER"] == "NGC6633?" and True,
    "IC4665": lambda row: row["FIELD/CLUSTER"] == "IC4665?" and True,
    "NGC104": lambda row: row["FIELD/CLUSTER"] == "NGC104?" and True,
    "gamma2_Vel": lambda row: row["FIELD/CLUSTER"] == "gamma2_Vel?" and True,
    "GJ880": lambda row: row["FIELD/CLUSTER"] == "GJ880?" and True,
    "NGC4815": lambda row: row["FIELD/CLUSTER"] == "NGC4815?" and True,
    "NGC2547": lambda row: row["FIELD/CLUSTER"] == "NGC2547?" and True,
    "NGC5927": lambda row: row["FIELD/CLUSTER"] == "NGC5927?" and True,
    "NGC4833": lambda row: row["FIELD/CLUSTER"] == "NGC4833?" and True,
    "NGC1851": lambda row: row["FIELD/CLUSTER"] == "NGC1851?" and True,
    "NGC2243": lambda row: row["FIELD/CLUSTER"] == "NGC2243?" and True,
    "NGC3532": lambda row: row["FIELD/CLUSTER"] == "NGC3532?" and True,
    "NGC6752": lambda row: row["FIELD/CLUSTER"] == "NGC6752?" and True,
    "Br25": lambda row: row["FIELD/CLUSTER"] == "Br25?" and True,
    "NGC4372": lambda row: row["FIELD/CLUSTER"] == "NGC4372?" and True,
    "NGC6705": lambda row: row["FIELD/CLUSTER"] == "NGC6705?" and True,
    "M67": lambda row: row["FIELD/CLUSTER"] == "M67?" and True,
    "NGC2516": lambda row: row["FIELD/CLUSTER"] == "NGC2516?" and True,
    "Trumpler20": lambda row: row["FIELD/CLUSTER"] == "Trumpler20?" and True,
}

# Assign candidates.
for cluster in clusters.keys():
    candidates = dataset.assign_cluster_candidates(cluster,
        lambda row: row["TARGET"].startswith(cluster))

    # Special hack:
    if cluster == "Trumpler20":
        candidates += dataset.assign_cluster_candidates(cluster,
            lambda row: row["TARGET"].startswith("Trumpler_20"))
    logger.info("Cluster {0} has {1} candidates".format(cluster, candidates))

# Assign members based on rules.
for cluster, membership_rule in clusters.items():
    members = dataset.assign_cluster_members(cluster, membership_rule)
    logger.info("Cluster {0} has {1} confirmed members".format(cluster, members))

# Remove benchmarks, etc.
# This will remove all stars that are not GES_* stars which were not labelled as
# being a cluster candidate or confirmed cluster member.
# Cluster stars have an assignment, MW stars will be marked as field. So we need
# to remove stars that have 'FIELD' and are not GES_MW_*, GES_CRT_*
stars_for_removal = np.array([row["FIELD/CLUSTER"] == "FIELD" and not
    (row["TARGET"].startswith("GES_MW") or row["TARGET"].startswith("GES_CRT"))\
    for row in dataset.data])
dataset.data = dataset.data[~stars_for_removal]
logger.info("Removed {0} benchmarks/misc stars from the sample. There are {1} "\
    "stars remaining".format(stars_for_removal.sum(), len(dataset.data)))

# Save the dataset.
dataset.write("data/ges-no-member-criteria.fits", overwrite=True)


# ASSIGNING REALISATIONS
# TEST 1 (TOY MODEL WITH POSITIONS AND? KINEMATICS)
#---------------------------------------------------
# Let's use the positional information (RA, DEC) to test the machinery. Set up
# many permutations using a minimum number of stars in a cluster to be 5 (or 10;
# whatever we need for it to work). Do this with GMM/AIC and VBGMM and DPGMM to
# examine the differences in the method. Set some maximum limit on the max number
# of stars in a cluster.

num_realisations = 10000
realisation_kwds = {
    "num_clusters": None, # Randomly select the number of clusters to use.
    "exclude_clusters": None, # Don't exclude any clusters
    "cluster_size_limits": (5, 50), # Select random number of stars from each cluster within these limits
    "field_star_fraction": 0, # Don't include field stars
    "cluster_star_constraints": None, # No constraints on cluster stars
    "field_star_constraints": None # No constraints on field stars
}

# Create the realisations for Test 1
logger.info("Creating {0} realisations".format(num_realisations))

def func(_):
    return tagging.realisations.create(dataset, **realisation_kwds)

pool = mp.Pool(processes=50)
all_realisations = []
def callback(_):
    global all_realisations
    all_realisations.append(_)

r = pool.map_async(func, xrange(int(num_realisations)), callback=callback)
r.wait()

# Winter is coming.
pool.close()
pool.join()

# Because Python.
all_realisations = all_realisations[0]

# Save the realisations for Test 1
with open("realisations/test_1.pickle", "wb") as fp:
    pickle.dump(all_realisations, fp, -1)
logger.info("Saved realisations to disk")