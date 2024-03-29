#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3
from numpy import array, arange, median, mean, abs
from numpy.linalg import norm
from pickle import load
from matplotlib import pyplot as plt
from matplotlib import rcParams
from sklearn.decomposition import PCA

from c3sdb.ml.data import C3SD

rcParams['font.size'] = 8


# load the trained predictive model (to use its internal KMeans classifier)
with open('kmcm_svr_final_20210414.pickle', 'rb') as pf:
    kmc = load(pf)


# load the dataset
data = C3SD('C3S.db', seed=2345)
data.assemble_features()
data.train_test_split('ccs')
data.center_and_scale()


# compute cluster assignments and distances
X_clusters = array([kmc.kmeans_.predict(x.reshape(1, -1))[0] for x in data.SScaler_.transform(data.X_)])
print('number of compounds in each cluster:', [len([l for l in kmc.kmeans_.labels_ if l == n]) for n in range(5)])


# compute cluster distance distributions
cluster_dists = [[], [], [], [], []]
cluster_dists_combined = []
for c, x in zip(X_clusters, data.SScaler_.transform(data.X_)):
    cluster_dists[c].append(norm(x - kmc.kmeans_.cluster_centers_[c]))
    cluster_dists_combined.append(norm(x - kmc.kmeans_.cluster_centers_[c]))
cluster_dists_combined = array(cluster_dists_combined)
fig = plt.figure(figsize=(3.33, 3.33))
ax = fig.add_subplot(111)
bins = arange(0, 20.001, 0.2)
y0 = 400
for cd, c in zip(cluster_dists, ['seagreen', 'darkkhaki', 'r', 'b', 'purple']):
    #plt.hist(cd, bins=bins, histtype='stepfilled', color=c, alpha=0.2)
    ax.hist(cd, bins=bins, histtype='step', color=c, lw=1.5)
    ax.text(10, y0, 'mean: {:.2f} median: {:.2f}'.format(mean(cd), median(cd)), fontweight='bold', fontsize=8, color=c)
    y0 -= 40
for d in ['top', 'right']:
    ax.spines[d].set_visible(False)
ax.set_xlabel('dist. to cluster center')
ax.set_ylabel('count')
plt.savefig('cluster_distance_distributions.png', dpi=400, bbox_inches='tight')
plt.close()


# plot cluster center distance against prediction error
X_scaled = data.SScaler_.transform(data.X_)
y_true = data.y_
percent_error = abs(100. * (kmc.predict(X_scaled) - y_true) / y_true)
fig = plt.figure(figsize=(3.33, 3.33))
ax = fig.add_subplot(111)
ax.scatter(cluster_dists_combined, percent_error, s=0.5, alpha=0.2, c='k', edgecolors='none')
for d in ['top', 'right']:
    ax.spines[d].set_visible(False)
ax.set_ylim([-0.05, 10])
ax.set_yticks([_ for _ in range(11)])
ax.set_xlim([0, 30])
ax.set_xlabel('dist. to cluster center')
ax.set_ylabel('abs. prediction error (%)')
plt.savefig('cluster_distance_vs_error.png', dpi=400, bbox_inches='tight')
plt.close()

# plot cluster center distance against prediction error (normalize distances by individual cluster medians)
dist_medians = array([[4.12, 8.11, 2.41, 4.47, 3.24][_] for _ in kmc.kmeans_.predict(X_scaled)])
fig = plt.figure(figsize=(3.33, 3.33))
ax = fig.add_subplot(111)
ax.scatter(cluster_dists_combined / dist_medians, percent_error, s=0.5, alpha=0.2, c='k', edgecolors='none')
for d in ['top', 'right']:
    ax.spines[d].set_visible(False)
ax.set_ylim([-0.05, 10])
ax.set_yticks([_ for _ in range(11)])
ax.set_xlim([0, 5])
ax.set_xlabel('dist. to cluster center (normalized)')
ax.set_ylabel('abs. prediction error (%)')
plt.savefig('cluster_distance_vs_error_normalized.png', dpi=400, bbox_inches='tight')
plt.close()


# plot distance from database center against prediction error
center = array([mean(_) for _ in X_scaled.T])
center_dists = array([norm(_ - center) for _ in X_scaled])
fig = plt.figure(figsize=(3.33, 3.33))
ax = fig.add_subplot(111)
ax.scatter(center_dists, percent_error, s=0.5, alpha=0.2, c='k', edgecolors='none')
for d in ['top', 'right']:
    ax.spines[d].set_visible(False)
ax.set_ylim([-0.05, 10])
ax.set_yticks([_ for _ in range(11)])
ax.set_xlim([0, 20])
ax.set_xlabel('dist. to center')
ax.set_ylabel('abs. prediction error (%)')
plt.savefig('center_distance_vs_error.png', dpi=400, bbox_inches='tight')
plt.close()


# plot distance to nearest neighbor against prediction error
# determine the nearest neighbors for each of the training data points
neighbor_dists = []
for i in range(data.N_):
    nds = []
    for j in range(data.N_):
        if j != i:
            nds.append(norm(X_scaled[i] - X_scaled[j]))
    neighbor_dists.append(min(nds))
fig = plt.figure(figsize=(3.33, 3.33))
ax = fig.add_subplot(111)
ax.scatter(neighbor_dists, percent_error, s=0.5, alpha=0.2, c='k', edgecolors='none')
for d in ['top', 'right']:
    ax.spines[d].set_visible(False)
ax.set_ylim([-0.05, 10])
ax.set_yticks([_ for _ in range(11)])
ax.set_xlim([0, 5])
ax.set_xlabel('dist. to nearest neighbor')
ax.set_ylabel('abs. prediction error (%)')
plt.savefig('neighbor_distance_vs_error.png', dpi=400, bbox_inches='tight')
plt.close()


# cluster distance matrix
dmat = [[0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0], [0, 0, 0, 0, 0]]
for i in range(5):
    for j in range(5):
        dmat[i][j] = norm(kmc.kmeans_.cluster_centers_[i] - kmc.kmeans_.cluster_centers_[j])
print('     c0    c1    c2    c3    c4')
print('c0 {:4.1f}  {:4.1f}  {:4.1f}  {:4.1f}  {:4.1f}'.format(dmat[0][0], dmat[0][1], dmat[0][2], dmat[0][3], dmat[0][4]))
print('c1       {:4.1f}  {:4.1f}  {:4.1f}  {:4.1f}'.format(dmat[1][1], dmat[1][2], dmat[1][3], dmat[1][4]))
print('c2             {:4.1f}  {:4.1f}  {:4.1f}'.format(dmat[2][2], dmat[2][3], dmat[2][4]))
print('c3                   {:4.1f}  {:4.1f}'.format(dmat[3][3], dmat[3][4]))
print('c4                         {:4.1f}'.format(dmat[4][4]))


# compute a 3-component PCA on the complete database using the full feature set
# use the StandardScaler to center and scale the data first
pca = PCA(n_components=0.95, svd_solver='full')
projections = pca.fit_transform(X_scaled).T
evr = pca.explained_variance_ratio_


# plot PC1 projection against prediction error
fig = plt.figure(figsize=(3.33, 3.33))
ax = fig.add_subplot(111)
ax.scatter(abs(projections[0]), percent_error, s=0.5, alpha=0.2, c='k', edgecolors='none')
for d in ['top', 'right']:
    ax.spines[d].set_visible(False)
ax.set_ylim([-0.05, 10])
ax.set_yticks([_ for _ in range(11)])
ax.set_xlim([0, 10])
ax.set_xlabel('abs. PC1 projection')
ax.set_ylabel('abs. prediction error (%)')
plt.savefig('PC1_vs_error.png', dpi=400, bbox_inches='tight')
plt.close()


# plot PC2 projection against prediction error
fig = plt.figure(figsize=(3.33, 3.33))
ax = fig.add_subplot(111)
ax.scatter(abs(projections[1]), percent_error, s=0.5, alpha=0.2, c='k', edgecolors='none')
for d in ['top', 'right']:
    ax.spines[d].set_visible(False)
ax.set_ylim([-0.05, 10])
ax.set_yticks([_ for _ in range(11)])
ax.set_xlim([0, 10])
ax.set_xlabel('abs. PC2 projection')
ax.set_ylabel('abs. prediction error (%)')
plt.savefig('PC2_vs_error.png', dpi=400, bbox_inches='tight')
plt.close()


# plot PC3 projection against prediction error
fig = plt.figure(figsize=(3.33, 3.33))
ax = fig.add_subplot(111)
ax.scatter(abs(projections[2]), percent_error, s=0.5, alpha=0.2, c='k', edgecolors='none')
for d in ['top', 'right']:
    ax.spines[d].set_visible(False)
ax.set_ylim([-0.05, 10])
ax.set_yticks([_ for _ in range(11)])
ax.set_xlim([0, 10])
ax.set_xlabel('abs. PC3 projection')
ax.set_ylabel('abs. prediction error (%)')
plt.savefig('PC3_vs_error.png', dpi=400, bbox_inches='tight')
plt.close()


cluster_center_pcs = [pca.transform([_])[0] for _ in kmc.kmeans_.cluster_centers_]
lab_to_c = ['seagreen', 'darkkhaki', 'r', 'b', 'purple']
c = [lab_to_c[_] for _ in kmc.kmeans_.predict(X_scaled)]
lab = [_ for _ in lab_to_c]


def pca_proj(pcA_n, pcB_n):
    pcA, pcB = projections[pcA_n - 1], projections[pcB_n - 1]
    pcA_evr, pcB_evr = evr[pcA_n - 1], evr[pcB_n - 1]
    fig = plt.figure(figsize=(3.33, 3.33))
    ax2 = fig.add_subplot(111)
    # PC A vs. PC B
    ax2.axvline(0, c='k', ls='--', lw=0.5)
    ax2.axhline(0, c='k', ls='--', lw=0.5)
    ax2.scatter(pcA, pcB, c=c, s=1, alpha=0.4, edgecolors='none')
    for d in ['top', 'left', 'bottom', 'right']:
        ax2.spines[d].set_visible(False)
    ax2.set_xlabel('PC{} ({:.1f} %)'.format(pcA_n, pcA_evr * 100))
    ax2.set_ylabel('PC{} ({:.1f} %)'.format(pcB_n, pcB_evr * 100))
    # add the cluster center projections onto the plots
    for i in range(5):
        cpA = cluster_center_pcs[i][pcA_n - 1]
        cpB = cluster_center_pcs[i][pcB_n - 1]
        ax2.scatter(cpA, cpB, c=lab_to_c[i], s=12, edgecolors='k', linewidths=0.7)
    plt.tight_layout()
    plt.savefig('KMC_cluster_PC{}{}_projections.png'.format(pcA_n, pcB_n), dpi=400, bbox_inches='tight')
    plt.close()


pca_proj(1, 2)
pca_proj(2, 3)
pca_proj(3, 4)


# plot m/z vs. CCS
fig = plt.figure(figsize=(3.33, 3.33))
ax = fig.add_subplot(111)
ccs_by_cluster = [[], [], [], [], []]
for l, cc in zip(kmc.kmeans_.labels_, data.ccs_):
    ccs_by_cluster[l].append(cc)
mean_cluster_ccs = [median(_) for _ in ccs_by_cluster]
print(mean_cluster_ccs)
mz_by_cluster = [[], [], [], [], []]
for l, m in zip(kmc.kmeans_.labels_, data.mz_):
    mz_by_cluster[l].append(m)
mean_cluster_mz = [median(_) for _ in mz_by_cluster]
print(mean_cluster_mz)
ax.scatter(data.mz_, data.ccs_, c=c, s=1, alpha=0.2, edgecolors='none')
for d in ['top', 'right']:
    ax.spines[d].set_visible(False)
ax.set_xlabel('m/z')
ax.set_ylabel(r'CCS (Å$^2$)')
ax.set_ylim([100, 380])
ax.set_xlim([50, 1200])
# add the cluster center projections onto the plots
for i in range(5):
    ax.scatter(mean_cluster_mz[i], mean_cluster_ccs[i],
               c=lab_to_c[i], s=16, edgecolors='k', linewidths=0.75)
plt.tight_layout()
plt.savefig('KMC_cluster_mz_v_ccs.png', dpi=400, bbox_inches='tight')
plt.close()
