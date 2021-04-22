#!/usr/local/Cellar/python@3.9/3.9.1_6/bin/python3
from numpy import loadtxt, linspace, array, argwhere, count_nonzero, mean
from matplotlib import pyplot as plt
from matplotlib import rcParams


rcParams['font.size'] = 8


def top_n_tc(tcs, n):
    return array([mean(_[:n]) for _ in tcs])

# load TC data from file
print('loading fp2048_tc')
fp2048_tc = loadtxt('fp2048_tc.npy')
print('loading fp1024_tc')
fp1024_tc = loadtxt('fp1024_tc.npy')
print('loading fp0512_tc')
fp0512_tc = loadtxt('fp0512_tc.npy')
print('loading percent_error')
percent_error = loadtxt('percent_error.npy')


def fdr3(tcs):
    tc_threshold = linspace(0, 1, 33)
    fdr = []
    for tct in tc_threshold:
        i_above_threshold = argwhere(tcs >= tct)
        pe = percent_error[i_above_threshold]
        tp = count_nonzero(pe <= 3)
        fp = count_nonzero(pe > 3)
        fdr.append(100. * float(fp) / float(fp + tp))
    return tc_threshold, array(fdr)


def plot_tc_fdr3(tc, fp_size,):
    fig = plt.figure(figsize=(2, 2))
    ax = fig.add_subplot(111)
    for top_n, c in zip([5, 10, 50, 100], ['#0000FF', '#0066AA', '#00AA66', '#008800']):
        fdr3_ = fdr3(top_n_tc(tc, top_n))
        ax.plot(*fdr3_, '-', c=c, lw=1.5, label='N = {}'.format(top_n))
    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    ax.legend(frameon=False)
    ax.set_ylim([0, 15])
    ax.set_xlim([0, 1.01])
    ax.set_xlabel('TC threshold (mean(top N TCs))')
    ax.set_ylabel('FDR % (3% error)')
    plt.savefig('fp{:04d}_fdr3.png'.format(fp_size), dpi=400, bbox_inches='tight')
    plt.close()


plot_tc_fdr3(fp0512_tc, 512)
plot_tc_fdr3(fp1024_tc, 1024)
plot_tc_fdr3(fp2048_tc, 2048)
