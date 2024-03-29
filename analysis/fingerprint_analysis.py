#!/usr/local/Cellar/python@3.9/3.9.1_6/bin/python3
from rdkit import Chem, DataStructs
from matplotlib import pyplot as plt
from matplotlib import rcParams
from pickle import load
from numpy import array, mean, argwhere, savetxt

from c3sdb.ml.data import C3SD

rcParams['font.size'] = 8


# load the trained predictive model (to use its internal KMeans classifier)
with open('kmcm_svr_final_20210414.pickle', 'rb') as pf:
    kmc = load(pf)


# load the dataset
print('loading and preparing dataset')
data = C3SD('C3S.db', seed=2345)
data.assemble_features()
data.train_test_split('ccs')
data.center_and_scale()

# get the cluster assignments
X_clusters = array([kmc.kmeans_.predict(x.reshape(1, -1))[0] for x in data.SScaler_.transform(data.X_)])

# scale all of the X data and make predictions, compute error
print('scaling and predicting CCS')
X_scaled = data.SScaler_.transform(data.X_)
percent_error = array(abs(100. * (kmc.predict(X_scaled) - data.y_) / data.y_))


"""
# remove all of the errors below 5%
i_above_5 = argwhere(percent_error > 5)
X_scaled = X_scaled[i_above_5]
percent_error = percent_error[i_above_5]


# dump the >5% error compounds to file
print('dumping >5% error data to csv')
with open('gt5_percent_error.csv', 'w') as f:
    f.write('cmpd,adduct,src,clust,mz,ccs,ccs_err\n')
    for cmpd, adduct, src, clust, mz, ccs, ccs_err in zip(data.cmpd_[i_above_5],
                                                          data.adduct_[i_above_5],
                                                          data.src_[i_above_5],
                                                          X_clusters[i_above_5],
                                                          data.mz_[i_above_5],
                                                          data.ccs_[i_above_5],
                                                          percent_error):
        #print(cmpd, adduct, src, clust, mz, ccs, ccs_err)
        f.write('"{}","{}","{}",{:d},{:.4f},{:2f},{:2f}\n'.format(cmpd[0], adduct[0], src[0], clust[0], mz[0], ccs[0], ccs_err[0]))
"""

# make rdkit MOL objects for all of the compounds in the database
print('making MOL objects')
#print(data.smi_)
mols = []
for smi in data.smi_:
    mols.append(Chem.MolFromSmiles(smi))


# compute a couple different sized fingerprints for each of the structures
print('computing fingerprints (2048)')
fp2048 = [Chem.RDKFingerprint(mol, fpSize=2048) for mol in mols]
print('computing fingerprints (1024)')
fp1024 = [Chem.RDKFingerprint(mol, fpSize=1024) for mol in mols]
print('computing fingerprints (512)')
fp0512 = [Chem.RDKFingerprint(mol, fpSize=512) for mol in mols]


# compute tanimoto coefficients for all of the fingerprints
print('comparing all fingerprints')
fp2048_tc, fp1024_tc, fp0512_tc = [], [], []
for i in range(len(fp0512)):
    scores_2048, scores_1024, scores_0512 = [], [], []
    for j in range(len(fp0512)):
        if j != i:  # do not compare a fingerprint to itself
            print('{:6d} {:6d}\r'.format(i, j), end='')
            scores_2048.append(DataStructs.FingerprintSimilarity(fp2048[i], fp2048[j]))
            scores_1024.append(DataStructs.FingerprintSimilarity(fp1024[i], fp1024[j]))
            scores_0512.append(DataStructs.FingerprintSimilarity(fp0512[i], fp0512[j]))
    fp2048_tc.append(sorted(scores_2048, reverse=True))
    fp1024_tc.append(sorted(scores_1024, reverse=True))
    fp0512_tc.append(sorted(scores_0512, reverse=True))
print()


# dump TC data to file
print('dumping TC data to file')
savetxt('fp2048_tc.npy', array(fp2048_tc))
savetxt('fp1024_tc.npy', array(fp1024_tc))
savetxt('fp0512_tc.npy', array(fp0512_tc))
savetxt('percent_error.npy', percent_error)



def top_n_tc(tcs, n):
    """ returns an array of the mean of the top N tanimoto coefficients for each fingerprint """
    return array([mean(_[:n]) for _ in tcs])


# compute top N tcs for all of the fingerprint data
print('computing top N TCs')
# fp2048
fp2048_top05_tcs = top_n_tc(fp2048_tc, 5)
fp2048_top10_tcs = top_n_tc(fp2048_tc, 10)
fp2048_top50_tcs = top_n_tc(fp2048_tc, 50)
# fp1024
fp1024_top05_tcs = top_n_tc(fp1024_tc, 5)
fp1024_top10_tcs = top_n_tc(fp1024_tc, 10)
fp1024_top50_tcs = top_n_tc(fp1024_tc, 50)
# fp0512
fp0512_top05_tcs = top_n_tc(fp0512_tc, 5)
fp0512_top10_tcs = top_n_tc(fp0512_tc, 10)
fp0512_top50_tcs = top_n_tc(fp0512_tc, 50)

"""
def plot_tcs(tcs, fp_size, top_n):
    fig = plt.figure(figsize=(3.33, 3.33))
    ax = fig.add_subplot(111)
    ax.scatter(tcs, percent_error, s=0.5, alpha=0.2, c='k', edgecolors='none')
    for d in ['top', 'right']:
        ax.spines[d].set_visible(False)
    #ax.set_ylim([-0.05, 10])
    #ax.set_yticks([_ for _ in range(11)])
    ax.set_xlim([0, 1])
    ax.set_xlabel('mean(top {} TCs)'.format(top_n))
    ax.set_ylabel('abs. prediction error (%)')
    plt.savefig('fp{:04d}_top{:02d}_tc_vs_error_above_5.png'.format(fp_size, top_n), dpi=400, bbox_inches='tight')
    plt.close()

# plot error vs. TCs
print('plotting error vs. TCs')
# fp2048
plot_tcs(fp2048_top01_tcs, 2048, 1)
plot_tcs(fp2048_top05_tcs, 2048, 5)
plot_tcs(fp2048_top10_tcs, 2048, 10)
plot_tcs(fp2048_top50_tcs, 2048, 50)
# fp1024
plot_tcs(fp1024_top01_tcs, 1024, 1)
plot_tcs(fp1024_top05_tcs, 1024, 5)
plot_tcs(fp1024_top10_tcs, 1024, 10)
plot_tcs(fp1024_top50_tcs, 1024, 50)
# fp0512
plot_tcs(fp0512_top01_tcs, 512, 1)
plot_tcs(fp0512_top05_tcs, 512, 5)
plot_tcs(fp0512_top10_tcs, 512, 10)
plot_tcs(fp0512_top50_tcs, 512, 50)
"""