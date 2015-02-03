#!/usr/bin/env python

""" Chemically tag the Gaia-ESO Survey data. """

from __future__ import absolute_import, print_function, with_statement

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

import code as ct

# PREPARATION
# Create some data set.

# Assign field stars based on rules.
# Assign cluster candidates based on rules.
# Assign cluster members based on rules.

# Write the data to disk.
#data.write("data/ges-data-set.fits")

# ASSIGNING REALISATIONS
# Create realisations using the data table.
# Save the information about the realisations to disk.

# DATA ANALYSIS
# For each realisation, we need to use particular amounts of data and try and
# work out the number of clusters present.

# Scripts for running on each realisation. Each script might be using a
# different amount of data (e.g., metallicities only?)




# SHOULD WE:
# Load data.
data = ct.data.DataSet.from_fits("data/GES_iDR2iDR3_WG10+WG11.fits", extension=1)

# Assign field/cluster stars.
# Unless otherwise told, this star is in the field:
data.data["FIELD/CLUSTER"] = "FIELD"

# Assign cluster stars based on rules.


# Save the data.


# Create some realisation where we use all cluster stars.
# Save the realisation information to disk.

# Infer the number of clusters in the realisation data set using AIC/BIC/XGMM/VGMM/XD.
# Save the information about the realisation, what data were used, and what was inferred.

# Then create more complex cluster candidate/rules.
# Then create more complex realisations.
