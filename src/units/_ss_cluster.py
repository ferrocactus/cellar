# Constrained clustering
#from copkmeans.cop_kmeans import cop_kmeans
#from active_semi_clustering.semi_supervised.pairwise_constraints import PCKMeans

def new_hard_cluster(self, labels, indices=None, all_points=False):
    """
    Update the label of points whose indices are given in indices,
    or it all_points flag is set, then assume all labels are given.
    """
    if all_points == False:
        assert(indices is not None)
        self.labels[indices] = labels
    else:
        assert(len(labels) == len(labels))
        self.labels = labels

def new_soft_cluster(self, point_index, k=None):
    """
    Given a single point, find the cluster where that point belongs
    determined by using elbow heuristics and update the labels.
    """
    distances = np.linalg.norm(self.x_2d_emb[point_index] - self.x_2d_emb, axis=1)
    sorted_indices = np.argsort(distances)
    plt.plot(range(0, len(distances)), distances[sorted_indices])
    knee = KneeLocator(range(0, len(distances)), distances[sorted_indices],
                        curve='convex',
                        direction='increasing',
                        S=2)
    print('knee at', knee.knee)
    if k is not None:
        self.labels[sorted_indices[:k]] = self.n_clusters
        self.n_clusters += 1
        self.find_markers()
        self.convert_markers()