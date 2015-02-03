#!/usr/bin/env python

""" Chemically tag the Gaia-ESO Survey dataset. """

from __future__ import absolute_import, print_function, with_statement

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

# Standard library.
import logging

# Module-specific.
import code as tagging

logger = logging.getLogger("code")

# PREPARATION
# Create some data set.

# Assign field stars based on rules.
# Assign cluster candidates based on rules.
# Assign cluster members based on rules.

# Write the data to disk.
#dataset.write("data/ges-data-set.fits")

# ASSIGNING REALISATIONS
# Create realisations using the data table.
# Save the information about the realisations to disk.

# DATA ANALYSIS
# For each realisation, we need to use particular amounts of data and try and
# work out the number of clusters present.

# Scripts for running on each realisation. Each script might be using a
# different amount of data (e.g., metallicities only?)




# SHOULD WE:
# Load dataset.
dataset = tagging.DataSet.from_fits("data/GES_iDR2iDR3_WG10+WG11.fits",
    extension=1)

# Assign field/cluster stars.
# Unless otherwise told, this star is in the field:
dataset.data["FIELD/CLUSTER"] = "FIELD"

# [TODO] Delete benchmarks
clusters = ("Cha_I", "Br81", "M15", "NGC2808", "NGC6633", "IC4665", "NGC104",
    "gamma2_Vel", "GJ880", "NGC4815", "NGC2547", "NGC5927", "NGC4833",
    "NGC1851", "NGC2243", "NGC3532", "NGC6752", "Br25", "NGC4372", "NGC6705",
    "M67", "NGC2516", "Trumpler20")

# [TODO] Assign candidates?
# Assign members.
for cluster in clusters:
    members = dataset.assign_cluster_members(cluster,
        lambda row: row["TARGET"].startswith(cluster))

    # Special hack:
    if cluster == "Trumpler20":
        members += dataset.assign_cluster_members(cluster,
            lambda row: row["TARGET"].startswith("Trumpler_20"))
    print("Cluster {0} has {1} members".format(cluster, members))

# Save the dataset.
#dataset.write("data/ges-no-member-criteria.fits")


rs_indices, rs_counts, rs_hash = tagging.realisations.create(dataset,
    exclude_clusters=["Br25"], num_clusters=5)

kwds = {
    "covariance_type": "full",
    "perturb_within_uncertainties": False
}

aic_num, model, aics = tagging.infer.cluster_count(dataset.data[rs_indices], ["RA", "DEC"],
    model="GMM/AIC", **kwds)

bic_num, model, bics = tagging.infer.cluster_count(dataset.data[rs_indices], ["RA", "DEC"],
    model="GMM/BIC", **kwds)

dpgmm_num, model = tagging.infer.cluster_count(dataset.data[rs_indices], ["RA", "DEC"],
    model="DPGMM", **kwds)

vbgmm_num, model = tagging.infer.cluster_count(dataset.data[rs_indices], ["RA", "DEC"],
    model="VBGMM", **kwds)

print(aic_num, bic_num, dpgmm_num, vbgmm_num)

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot(np.arange(aics.size) + 1, aics, c="b", linestyle="-.", label="GMM/AIC")
ax.plot(np.arange(bics.size) + 1, bics, c="k", linestyle="-.", label="GMM/BIC")

raise a

"""
data, num_clusters=np.inf, exclude_clusters=None,
    cluster_size_limits=None, field_star_fraction=0, 
    cluster_star_constraints=None, field_star_constraints=None,
"""
# Create some realisation where we use all cluster stars.
# Save the realisation information to disk.

# Infer the number of clusters in the realisation data set using AIC/BIC/XGMM/VGMM/XD.
# Save the information about the realisation, what data were used, and what was inferred.

# Then create more complex cluster candidate/rules.
# Then create more complex realisations.
