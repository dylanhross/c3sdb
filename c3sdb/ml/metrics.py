"""
    c3sdb/ml/metrics.py

    Dylan Ross (dylan.ross@pnnl.gov)

    module with functions to compute standard metrics for evaluating
    ML CCS prediction model
"""

from typing import Dict, List


import numpy as np
from numpy import typing as npt
from sklearn.metrics import r2_score, mean_squared_error
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec


def compute_metrics(y: npt.NDArray[np.float64], 
                    y_pred: npt.NDArray[np.float64]
                    ) -> Dict[str, float | List[float]] :
    """
    compute a standard array of metrics between a set of predicted and reference values
    
    - R-squared (R2)
    - root mean squared error (RMSE)
    - mean absolute error (MAE)
    - median absolute erre (MDAE)
    - mean relative error % (MRE)
    - median relative error % (MDRE)
    - cumulative error distribution at the <1, <3, <5, and <10% levels (CE135A)

    Parameters
    ----------
    y : ``numpy.ndarray(float)``
    y_pred : ``numpy.ndarray(float)``
        arrays of reference and predicted values

    Returns
    -------
    summary : ``dict(...)``
        dict with set of metrics
    """
    abs_y_err = np.abs(y_pred - y)
    r2 = r2_score(y, y_pred)
    mae = np.mean(abs_y_err)
    mdae = np.median(abs_y_err)
    mre = np.mean(100. * abs_y_err / y)
    mdre = np.median(100. * abs_y_err / y)
    rmse = np.sqrt(mean_squared_error(y, y_pred))
    y_err_percent = 100. * abs_y_err / y
    cum_err = np.cumsum(np.histogram(y_err_percent, [_ for _ in range(101)])[0])
    cum_err = 100. * cum_err / np.sum(cum_err)
    ce1, ce3, ce5, ceA = cum_err[0], cum_err[2], cum_err[4], cum_err[9]
    return {
        'R2': r2, 'MAE': mae, 'MDAE': mdae, 'MRE': mre, 'MDRE': mdre, 
        'RMSE': rmse, 'CE135A': [ce1, ce3, ce5, ceA]
    }


def compute_metrics_train_test(y_train: npt.NDArray[np.float64],
                               y_test: npt.NDArray[np.float64],
                               y_pred_train: npt.NDArray[np.float64],
                               y_pred_test: npt.NDArray[np.float64]
                               ) -> Dict[str, Dict[str, float | List[float]]] :
    """
    computes a standard set of performance metrics from separate predictions
    for a training and test dataset
    - training and test set R-squared (R2)
    - training and test set root mean squared error (RMSE)
    - training and test set mean absolute error (MAE)
    - training and test set median absolute error (MDAE)
    - training and test set mean relative error (MRE)
    - training and test set median relative error (MDRE)
    - cumulative error distribution at the <1, <3, <5, and <10% levels
        for training and test set (CE135A)
    
    Parameters
    ----------
    y_train : ``numpy.ndarray(float)``
    y_test : ``numpy.ndarray(float)``
    y_pred_train : ``numpy.ndarray(float)``
    y_pred_test : ``numpy.ndarray(float)``
        arrays of reference and predicted values for training and test datasets

    Returns
    -------
    summary : ``dict(...)``
        dict with set of metrics for training and test datasets
    """
    summary = {}
    for y, y_pred, lbl in [(y_train, y_pred_train, "train"), (y_test, y_pred_test, "test")]:
        # compute metrics
        summary[lbl] = compute_metrics(y, y_pred)
    return summary


def train_test_summary_figure(summary: Dict[str, Dict[str, float | List[float]]], 
                              fig_name: str, 
                              r2_range=[0.95, 1.]):
    """
    produces a summary figure displaying the results from `compute_metrics_train_test`
    
    Parameters
    ----------
    summary : ``dict(...)``
        summary dict returned from `compute_metrics_train_test`
    fig_name : ``str``
        file name to save the generated plot under 
    r2_range : ``list(float))``, default=[0.95, 1.]]
        lower and upper bounds of R-squared y axis 
"""
    fig = plt.figure(figsize=(5, 3))
    gs = GridSpec(1, 4, figure=fig, width_ratios=[1.2, 3, 2, 5])
    # R-squared
    ax1 = fig.add_subplot(gs[0])
    rsq_trn = summary['train']['R2']
    rsq_tst = summary['test']['R2']
    w1 = 0.15
    ax1.bar([1 - w1 / 2, 1 + w1 / 2], [rsq_trn, rsq_tst], color=['b', 'r'], width=w1)
    for d in ['top', 'right']:
        ax1.spines[d].set_visible(False)
    ax1.set_xticks([])
    ax1.set_ylabel(r'R$^2$')
    ax1.set_ylim(r2_range)
    ax1.set_xlim([0.75, 1.25])
    # MAE, MDAE and RMSE
    ax2 = fig.add_subplot(gs[1])
    mae_trn = summary['train']['MAE']
    mae_tst = summary['test']['MAE']
    mdae_trn = summary['train']['MDAE']
    mdae_tst = summary['test']['MDAE']
    mse_trn = summary['train']['RMSE']
    mse_tst = summary['test']['RMSE']
    ax2.bar([0.875, 1.125], [mae_trn, mae_tst], color=['b', 'r'], width=0.25)
    ax2.bar([1.875, 2.125], [mdae_trn, mdae_tst], color=['b', 'r'], width=0.25)
    ax2.bar([2.875, 3.125], [mse_trn, mse_tst], color=['b', 'r'], width=0.25)
    for d in ['top', 'right']:
        ax2.spines[d].set_visible(False)
    ax2.set_xticks([1, 2, 3])
    ax2.set_xticklabels(['MAE', 'MDAE', 'RMSE'], rotation='vertical')
    ax2.set_ylabel(r'CCS (Å$^2$)')
    # CE135A
    ax3 = fig.add_subplot(gs[3])
    x1 = [_ - 0.125 for _ in range(1, 5)]
    y1 = [100. * summary['train']['CE135A'][i] for i in range(4)]
    x2 = [_ + 0.125 for _ in range(1, 5)]
    y2 = [100. * summary['test']['CE135A'][i] for i in range(4)]
    ax3.bar(x1, y1, color='b', width=0.25)
    ax3.bar(x2, y2, color='r', width=0.25)
    for d in ['top', 'right']:
        ax3.spines[d].set_visible(False)
    ax3.set_xlabel('pred. error (%)')
    ax3.set_xticks([1, 2, 3, 4])
    ax3.set_xticklabels(['<1', '<3', '<5', '<10'])
    ax3.set_ylabel('proportion (%)')
    # between axes 2 and 3
    # MRE, MDRE
    ax23 = fig.add_subplot(gs[2])
    mre_trn = summary['train']['MRE']
    mre_tst = summary['test']['MRE']
    mdre_trn = summary['train']['MDRE']
    mdre_tst = summary['test']['MDRE']
    w23 = 0.25
    ax23.bar([1 - w23 / 2, 1 + w23 / 2], [mre_trn, mre_tst], color=['b', 'r'], width=w23)
    ax23.bar([2 - w23 / 2, 2 + w23 / 2], [mdre_trn, mdre_tst], color=['b', 'r'], width=w23)
    for d in ['top', 'right']:
        ax23.spines[d].set_visible(False)
    ax23.set_xticks([1, 2])
    ax23.set_xticklabels(['MRE', 'MDRE'], rotation='vertical')
    ax23.set_ylabel('%')
    plt.tight_layout()
    plt.savefig(fig_name, dpi=400, bbox_inches='tight')
