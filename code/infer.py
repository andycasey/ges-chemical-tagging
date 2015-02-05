#!/usr/bin/env python

""" Infer the number of clusters present in the given data realisation. """

from __future__ import absolute_import, print_function, with_statement

__author__ = "Andy Casey <arc@ast.cam.ac.uk>"

# Standard library.
import logging
import threading
import warnings

# Third-party.
import numpy as np
from astropy.io import fits
from astropy.table import Column, Table
from scipy import linalg
from sklearn import mixture

# Create loggers.
logger = logging.getLogger(__name__)

class InsufficientSamplesWarning(Warning):
    pass
warnings.simplefilter("once", InsufficientSamplesWarning)

def cluster_count_by_dpgmm(data, max_components=100, **kwargs):
    """
    Determine the number of clusters by fitting a Dirichlet Process Gaussian
    Mixture Model to the data.

    :param data:
        The data array.

    :type data:
        :class:`numpy.array`
    """
    logger.debug("Inferring cluster number by DPGMM")

    if max_components > len(data):
        warnings.warn("Maximum number of components exceeded the sample size. S"
            "etting the maximum number of components equal to the sample size."\
            .format(max_components, len(data)), InsufficientSamplesWarning)
        max_components = len(data)

    # Extras.
    full_output = kwargs.pop("full_output", False)
    mp_queue = kwargs.pop("__mp_queue", False)

    # Defaults.
    kwargs.setdefault("covariance_type", "full")

    # Fit stuff.
    model = mixture.DPGMM(max_components, **kwargs)
    try:
        model.fit(data)
    except ValueError:
        logger.exception("Failed to infer cluster count by DPGMM:")
        if full_output:
            if mp_queue: mp_queue.put((-1, model))
            return (-1, model)
        if mp_queue: mp_queue.put(-1)
        return -1

    # Now actually find out how many components were used
    Y, components = model.predict(data), 0
    for i, (mean, covar) in enumerate(zip(model.means_, model._get_covars())):
        v, w = linalg.eigh(covar)
        u = w[0] / linalg.norm(w[0])
        if np.any(Y == i):
            components += 1

    if full_output:
        if mp_queue: mp_queue.put((components, model))
        return (components, model)
    if mp_queue: mp_queue.put(components)
    return components


def cluster_count_by_vbgmm(data, max_components=100, **kwargs):
    """
    Determine the number of clusters by using a Variational Inference for
    Gaussian Mixture Model with the data.

    :param data:
        The data array.

    :type data:
        :class:`numpy.array`
    """
    logger.debug("Inferring cluster number by VBGMM")

    if max_components > len(data):
        warnings.warn("Maximum number of components exceeded the sample size. S"
            "etting the maximum number of components equal to the sample size."\
            .format(max_components, len(data)), InsufficientSamplesWarning)
        max_components = len(data)

    # Extras.
    full_output = kwargs.pop("full_output", False)
    mp_queue = kwargs.pop("__mp_queue", False)

    # Defaults.
    kwargs.setdefault("covariance_type", "full")

    # Fit stuff.
    model = mixture.VBGMM(max_components, **kwargs)
    try:
        model.fit(data)
    except ValueError:
        logger.exception("Failed to infer cluster count by VBGMM:")
        if full_output:
            if mp_queue: mp_queue.put((-1, model))
            return (-1, model)
        if mp_queue: mp_queue.put(-1)
        return -1

    # Now actually find out how many components were used
    Y, components = model.predict(data), 0
    for i, (mean, covar) in enumerate(zip(model.means_, model._get_covars())):
        v, w = linalg.eigh(covar)
        u = w[0] / linalg.norm(w[0])
        if np.any(Y == i):
            components += 1

    if full_output:
        if mp_queue: mp_queue.put((components, model))
        return (components, model)
    if mp_queue: mp_queue.put(components)
    return components


def cluster_count_by_gmm(data, max_components=100, metric=None, **kwargs):
    """
    Determine the number of clusters by fitting Gaussian mixture models and 
    minimising some metric (e.g., AIC/BIC).

    :param data:
        The data array.

    :type data:
        :class:`numpy.array`

    :param max_components: [optional]
        The maximum number of components to consider. By default this will keep
        trying mixtures for 5 components past the minimum. This implicitly
        assumes there is a single minimum in the metric space. This value can
        be adjusted by providing `n_components_past_minimum`.

    :type max_components:
        int

    :param metric: [optional]
        The metric to use to determine how many components there are. Options
        are 'AIC' or 'BIC'. Default is both, in order of AIC, BIC.

    :type metric:
        str
    """
    logger.debug("Inferring cluster number by GMM {}".format(metric))

    if isinstance(metric, (str, unicode)):
        metric = metric.upper()

    # Remove unnecessary keywords.
    full_output = kwargs.pop("full_output", False)
    num_past_minimum = kwargs.pop("n_components_past_minimum", 5)
    kwargs.pop("n_components", None)

    # Defaults:
    kwargs.setdefault("covariance_type", "full")

    aics = []
    bics = []
    for i in range(1, max_components + 1):
        # Fit the data
        model = mixture.GMM(n_components=i, **kwargs)

        try:
            model.fit(data)
        except ValueError:
            logger.exception("Failed to infer cluster count by GMM/{}:".format(
                metric))
            if full_output:
                return [-1, -1, aics, bics]
            return [-1, -1]
        
        # Calculate statistics
        aics.append(model.aic(data))
        bics.append(model.bic(data))

        # Do we need to keep going?
        if metric == "AIC" \
            and (aics[-1] > min(aics) and i - np.argmin(aics) > num_past_minimum) \
        or metric == "BIC" \
            and (bics[-1] > min(bics) and i - np.argmin(bics) > num_past_minimum) \
        or metric == None \
            and (aics[-1] > min(aics) and i - np.argmin(aics) > num_past_minimum \
            and bics[-1] > min(bics) and i - np.argmin(bics) > num_past_minimum):
                break

    aics, bics = map(np.array, (aics, bics))
    if aics.size == max_components:
        logger.warn("Maximum number of components ({}) reached!".format(
            max_components))

    # Which one should we be returning on?
    num_clusters_aic = np.argmin(aics) + 1
    num_clusters_bic = np.argmin(bics) + 1

    if full_output:
        return [num_clusters_aic, num_clusters_bic, aics, bics]
    return [num_clusters_aic, num_clusters_bic]



__available_models = ("GMM", "GMM/AIC", "GMM/BIC", "DPGMM", "VBGMM")
def cluster_count(stars, data_columns, model, perturb_within_uncertainties=False,
    **kwargs):
    """
    Infer the number of star clusters in the given data.

    :param stars:
        The table of stars, containing abundance information.

    :type stars:
        :class:`astropy.table.Table`

    :param data_columns:
        Which columns from `stars` should be used to inform us about the number
        of star clusters in the data.

    :type data_columns:
        list of str

    :param model:
        The model to use to infer the number of star clusters. The available
        models are: {models}

    :type model:
        str

    :param perturb_within_uncertainties: [optional]
        If set to `True`, then the values within `data_columns` will be
        perturbed about the corresponding error values (assuming 'COLUMN' has
        an error column 'E_COLUMN'). Uncertainties are assumed to be normally
        distributed, and it is assumed the 'E_COLUMN' values represent the
        1-sigma uncertainty.

    :type perturb_within_uncertainties:
        bool
    """

    prefix = kwargs.pop("__mp_return_prefix", None)
    result = []
    if prefix is not None:
        result.append(prefix)
        logger.debug("Inferring on parallel job #{}".format(prefix))
    
    if not isinstance(stars, Table):
        raise TypeError("stars must be a :class:`astropy.table.Table` object")

    if not isinstance(data_columns, (tuple, list, np.array)):
        raise TypeError("data columns expected to be a tuple of str")

    # Check the columns exist.
    for c in data_columns:
        if c not in stars.dtype.names:
            raise ValueError("column '{}' is not in the stars table".format(c))
        if perturb_within_uncertainties \
        and "E_{}".format(c) not in stars.dtype.names:
            raise ValueError("cannot perturb within uncertainties for column "
                "'{0}' because it does not have an error column 'E_{0} exist "
                "in the star table".format(c))

    if model.lower() not in map(str.lower, __available_models):
        raise ValueError("model '{0}' is not recognised; available models are"
            " {1}".format(model, ", ".join(__available_models)))

    # Build data array
    data = np.zeros((len(stars), len(data_columns)), dtype=float)
    for i, c in enumerate(data_columns):
        if perturb_within_uncertainties:
            data[:, i] = np.random.normal(stars[c], stars["E_{}".format(c)])
        else:
            data[:, i] = stars[c][:]

    kwds = kwargs.copy()
    kwds.pop("metric", None) # Remove metric in case we have to specify it.

    # Do the fitting.
    if model == "GMM":
        _ = cluster_count_by_gmm(data, **kwds)
    elif model == "GMM/AIC":
        _ = cluster_count_by_gmm(data, metric="AIC", **kwds)
    elif model == "GMM/BIC":
        _ = cluster_count_by_gmm(data, metric="BIC", **kwds)
    elif model == "DPGMM":
        _ = cluster_count_by_dpgmm(data, **kwds)
    elif model == "VBGMM":
        _ = cluster_count_by_vbgmm(data, **kwds)
    else:
        raise WTFError()

    result += list(_)

    return result

# Update the docstring
cluster_count.__doc__ = cluster_count.__doc__.format(
    models=", ".join(__available_models))