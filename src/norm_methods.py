#!/usr/bin/env python3
import numpy as np
import scipy.stats
import sys

def combine_std_err(stderr_t, stderr_u):
    return math.sqrt(stderr_t**2 + stderr_u**2)

def np_combine_std_err(stderr_t, stderr_u):
    return math.sqrt(stderr_t**2 + stderr_u**2)

def calc_global_norm(x):
    """ boxplot outlier removal with upper 90% normalization """
    qs = np.quantile(x, [0.25, 0.75])
    iqr = np.diff(qs)[0]
    high_cutoff = qs[1] + (1.5 * iqr)
    
    # based on shapemapper approach, use either outlier or top 10% of
    # data, whichever removes the fewest points

    q = sum(x >= high_cutoff)
    ten_percent = sum(x > np.quantile(x, 0.90))
    
    if (q <= ten_percent):
      low_cutoff = np.quantile(x[x < high_cutoff], 0.90)
      norm_factor = np.mean(x[(x > low_cutoff) & (x < high_cutoff)])
      print("normalization factor = ", 
            norm_factor,
            " using boxplot",
            " with ", q, " removed",
            file = sys.stderr)
    else:
      idx = (x <= np.quantile(x, 0.90)) 
      # remove top 10%
      x = x[idx]
      # keep next top 10% for normalization
      new_idx = (x >= np.quantile(x, 0.90))
      norm_factor = np.mean(x[new_idx])
      print("normalization factor = ", 
            norm_factor,
            " using top 10% removal ",
            " with ", ten_percent, " removed",
            file = sys.stderr)
    return norm_factor

def windsorize(x, lim = 0.95):
    """
    Windsorize values
    """
    return scipy.stats.mstats.winsorize(x, limits = (None, 1.0 - lim)).data


def get_constraints(x, cutoff = 0.20, special_value = -999):
    """
    return a list formatted to indicate constrained bases
    for RNAfold
    """
    cutoff_value = max(x) * cutoff
    constrained_nts = np.logical_and(x >= cutoff_value,
                                     x != special_value)
    res = np.where(constrained_nts, "x", ".")
    
    return res.tolist()


