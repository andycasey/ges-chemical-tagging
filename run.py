#!/usr/bin/env python

""" Chemically tag the Gaia-ESO Survey dataset. """

from __future__ import absolute_import, print_function, with_statement

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

# Standard library.
import cPickle as pickle
import logging
import multiprocessing as mp


# Third-party.
import numpy as np
from astropy.table import Table

# Module-specific.
import code as tagging

# Create loggers.
logger = logging.getLogger("code")
mp_logger = mp.log_to_stderr()
mp_logger.setLevel(mp.SUBDEBUG)


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
# Transfer over .matplotlibrc file and update as necessary
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

# Load the stuff
dataset = tagging.DataSet.from_fits("data/ges-no-member-criteria.fits",
    extension=1)
with open("realisations/test_1.pickle", "rb") as fp:
    all_realisations = pickle.load(fp)

# Specify the keywords for Test 1 w/ AIC
test_1_kwds = {
    "model": "GMM/AIC", # Available: GMM/AIC, GMM/BIC, DPGMM, VBGMM
    "data_columns": ["RA", "DEC"],
    "full_output": False,
    "covariance_type": "full",
    "perturb_within_uncertainties": False,
}

output = []
pool = mp.Pool(processes=50)
def callback(_):
    output.append(_)

for i, (indices, true_clusters, unique_hash) in enumerate(all_realisations):
    kwds = test_1_kwds.copy()
    kwds["__mp_return_prefix"] = i
    pool.apply_async(tagging.infer.cluster_count, args=(dataset.data[indices], ),
        kwds=kwds, callback=callback)

# Winter is coming.
pool.close()
logger.info("Pool closed. Joining...")

try:
    pool.join()

except KeyboardInterrupt:
    # Screw it. Use as is!
    logger.warn("Assuming we have timed out. Pool is not closed!")

logger.info("Collating results together...")

# Plot stuff from the results? Nah...
results = []
for each in output:
    realisation = all_realisations[each[0]]
    # Unique hash, true clusters, inferred clusters.
    results.append([realisation[2], len(realisation[1]), each[1]])

# Create a table of the results from this realisation.
results = Table(data=np.array(results), names=("hash", "N_true", "N_AIC"),
    dtype=("S32", "i4", "i4"))

# Write the test results to disk.
results.write("results/test_1_AIC.fits", overwrite=True)

# Make plots.
fig, ax = plt.subplots()
ok = results["N_AIC"] > 0
ax.scatter(results["N_true"][ok], results["N_AIC"][ok]-results["N_true"][ok],
    facecolor="k")

raise a

"""

# Now let's try and infer the number of clusters in each realisation.
test_results = []
logger.info("Inferring cluster numbers from realisations")
times = []
for i, (indices, true_clusters, unique_hash) in enumerate(all_realisations):

    t_init = time()

    # Infer the number of clusters using AIC.
    inferred_clusters_by_aic, gmm_model_aic, aics = tagging.infer.cluster_count(
        dataset.data[indices], model="GMM/AIC", **test_1_kwds)

    # Infer the number of clusters using BIC.
    inferred_clusters_by_bic, gmm_model_bic, bics = tagging.infer.cluster_count(
        dataset.data[indices], model="GMM/BIC", **test_1_kwds)

    # Infer the number of clusters using a DPGMM.
    inferred_clusters_by_dpgmm, dpgmm_model = tagging.infer.cluster_count(
        dataset.data[indices], model="DPGMM", **test_1_kwds)

    # Infer the number of clusters using a VBGMM.
    inferred_clusters_by_vbgmm, vbgmm_model = tagging.infer.cluster_count(
        dataset.data[indices], model="VBGMM", **test_1_kwds)

    t_init2 = time()
    print("serial {}".format(t_init2 - t_init))

    results = tagging.infer.parallel_cluster_count(dataset.data[indices], **test_1_kwds)

    times.append([t_init2 - t_init, time() - t_init2])

    print(i, times[-1])
    continue
    # [TODO] Store the AICS and BICS?
    fig, ax = plt.subplots(2)
    # Take first two data columns only
    x = dataset.data[test_1_kwds["data_columns"][0]][indices]
    y = dataset.data[test_1_kwds["data_columns"][1]][indices]
    ax[0].scatter(x, y, facecolor="k")
    ax[0].set_xlabel(test_1_kwds["data_columns"][0])
    ax[0].set_ylabel(test_1_kwds["data_columns"][1])
    ax[1].plot(np.arange(len(aics)) + 1, aics, c="r", lw=2, label="AIC")
    ax[1].plot(np.arange(len(bics)) + 1, bics, c="b", lw=2, label="BIC")
    ax[1].scatter(inferred_clusters_by_aic, aics[inferred_clusters_by_aic - 1], facecolor="r", s=50)
    ax[1].scatter(inferred_clusters_by_bic, bics[inferred_clusters_by_bic - 1], facecolor="b", s=50)
    ax[1].axvline(len(true_clusters), c="k", lw=2)
    ax[1].axvline(inferred_clusters_by_dpgmm, c="k", ls=":", label="DPGMM")
    ax[1].axvline(inferred_clusters_by_vbgmm, c="K", ls="-.", label="VBGMM")
    ax[1].legend()
    ax[1].set_xlim(1, max(map(len, [aics, bics])) + 1)
    fig.savefig("figures/test1-{}.png".format(unique_hash))
    plt.close("all")

    # Store the results for this realisation.
    row = [
        unique_hash, len(true_clusters), inferred_clusters_by_aic,
        inferred_clusters_by_bic, inferred_clusters_by_dpgmm,
        inferred_clusters_by_vbgmm
    ]
    logger.debug("Results from realisation {0}: {1}".format(i, row))
    test_results.append(row)

raise a
"""



"""

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