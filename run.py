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

# Remove *all* benchmarks!
# Assign field stars based on rules (e.g., GES_MW_*)
# Assign cluster candidates based on rules (by TARGET)
# Assign cluster members based on rules (by TARGET + VELs)

# Write the data to disk.
#dataset.write("data/ges-data-set.fits")

# ASSIGNING REALISATIONS
# TEST 1 (TOY MODEL WITH POSITIONS AND? KINEMATICS)
#---------------------------------------------------
# Let's use the positional information (RA, DEC) to test the machinery. Set up
# many permutations using a minimum number of stars in a cluster to be 5 (or 10;
# whatever we need for it to work). Do this with GMM/AIC and VBGMM and DPGMM to
# examine the differences in the method. Set some maximum limit on the max number
# of stars in a cluster.


# TEST 2 (METALLICITIES ONLY, AS A FUNCTION OF N_CLUSTERS)
#----------------------------------------------------------
# Chose different numbers of clusters to include in the set, and compare how
# many that are recovered, using just [Fe/H]. Use the same (lower, upper) limits
# on stars/cluster as in Test 1.


# TEST 3 (METALLICITIES + ALPHAS)
#---------------------------------
# As per test 2, but use metallicities and alphas.
# (consider doing with uncertainties as well)


# TEST 4 (METALLICITIES + ALPHAS + Fe-peak)
#-------------------------------------------
# As per test 2, but use met + alpha + fe-p.
# (considering doing by perturbing uncertainties as well)


# TEST 5 (ALL ABUNDANCES)
#--------------------------------------------
# (perturb uncertainties as well)


# TEST 6 (TAKE THE BEST COMBINATION OF ABUNDANCES, NO UNCERTAINTIES, ADD FIELD STARS)
#-------------------------------------------------------------------------------------

# TEST 7 (BEST COMBINATION OF ABUNDANCES, W/ UNCERTAINTIES, ADD FIELD STARS)
#----------------------------------------------------------------------------

# TEST 8 (MINIMUM NUMBER OF CLUSTER MEMBERS?)
#---------------------------------------------
# Is there a minimum number of stars per cluster? How do the above tests
# change if we set the minimum stars in a cluster to be like 5, 10, 20, or 50?
# How many clusters do we end up actually having?

# TEST 9 (MINIMUM PRECISION REQUIRED FROM SOME TOY MODEL?)
#----------------------------------------------------------
# Take the mean abundance distributions of our clusters & approximate the covariance
# between them. Then we could say "OK, some cluster might have these properties,
# and the abundance differences within that cluster are immeasurable" -- and we
# could ask the question "how well would we need to measure these abundances to
# correctly infer their separation in amongst the field?"

# TEST 9 is pretty much doing all the other tests using faux data, based on the
# properties of the data. So we might need some code to create faux data based
# on the properties of the real data....



# Create realisations using the data table.
# Save the realisations to disk for later.

# DATA ANALYSIS
# Create functions to run on a bunch of realisations, using some method, and
# some columns of data >> save the results from each realisation to disk.

# PLOTTING
# Create way to plot hexbins/etc in some easy way.
# Create way to plot discretised Ncluster_true along the x-axis too.

# Make plots for all tests

# OPEN QUESTIONS:
# (1) Should any clusters be excluded because they are not really representative
#     of the GALAH-like sample?
# (2) How sensitive are the test results to the maximum number of stars in any
#     given cluster?
# (3) With the faux data, how do the tests scale to very large Ncluster?
# (4) Do we need to really use WAIC?


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
    "full_output": True,
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
