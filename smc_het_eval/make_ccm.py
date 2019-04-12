from SMCScoring import *
import numpy as np
import csv
import pdb
# np.set_printoptions(threshold='nan')

def closest_rand_reassign(in_clusters,p_err=0.1):
    clusters = np.copy(in_clusters)
    n_clusters = in_clusters.shape[1]
#    print(range(in_clusters.shape[0]))
 #   print(p_err)
    for j in range(in_clusters.shape[0]):
        if np.random.random() < p_err:
            cluster = np.argmax(clusters[j,:])
         #   print("cluster", cluster)
            if cluster == (n_clusters - 1):
                cluster = (n_clusters - 2)
            elif cluster == 0:
                cluster = 1
            else:
                if np.random.random() < 0.5:
                    cluster = (cluster + 1)
                else:
                    cluster = (cluster - 1)
            #cluster = cluster % n_clusters
           # print("cluster modulo", cluster)
            clusters[j,:] = 0
            clusters[j,cluster] = 1
  #  print("rand clusters")
  #  print(clusters)
    return( np.dot(clusters,clusters.T), clusters,)



def get_ccm_nvar(scenario, t_ccm=None, t_clusters=None, size_clusters=[100,100,100,100,100,100], prop_split=0.5, n_clusters=6, extra_prop=1.0/12.0, nssms=None):
#    print("in get_ccm", scenario)
    # pdb.set_trace()
    size_clusters = np.asarray(size_clusters)
 #   print(size_clusters)
    if t_clusters is None:
        t_clusters = np.zeros((np.sum(size_clusters),n_clusters))
    #    print(t_clusters.shape)
        for i in range(n_clusters):
       #     print(i)
       #     print(np.sum(size_clusters[0:i+1])-1)
            t_clusters[np.sum(size_clusters[0:i]):np.sum(size_clusters[0:i+1]),i] = 1 #assign each cluster
    if t_ccm is None and nssms is None:
        t_ccm = np.dot(t_clusters,t_clusters.T)
    #    print(t_clusters)
    #    print(t_ccm)
    if scenario == "Truth":
        clusters = np.copy(t_clusters)
        return (np.dot(clusters, clusters.T), clusters)
    elif "ParentIs" in scenario:
        clusters = np.copy(t_clusters)
        return (np.dot(clusters, clusters.T), clusters)
    elif scenario == 'OneCluster':
        if nssms is None:
        #    print(t_ccm.shape)
            out = np.ones(t_ccm.shape)
      #      print("made out")
            return (out, np.ones((t_ccm.shape[0],1)))
      #  return np.ones((nssms,nssms), dtype=np.int8)
    elif "NCluster" in scenario:
        if nssms is None:
            return (np.identity(t_ccm.shape[0]), np.identity(t_ccm.shape[0]))
        return (np.identity(nssms, dtype=np.int8),np.identity(nssms, dtype=np.int8))
    elif "SplitCluster" in scenario:
        clusters = np.zeros((np.sum(size_clusters),n_clusters+1))
        clusters[:,:-1] = np.copy(t_clusters)
        if scenario == "SplitClusterMidOneChild":
            clusters[np.sum(size_clusters[0:2])+int(size_clusters[2]/(1/(1-prop_split))):np.sum(size_clusters[0:3]),2] = 0
            clusters[np.sum(size_clusters[0:2])+int(size_clusters[2]/(1/(1-prop_split))):np.sum(size_clusters[0:3]),n_clusters] = 1
     #       print(clusters)
           # print(np.dot(clusters, clusters.T))
            return (np.dot(clusters,clusters.T),clusters)
        elif scenario == "SplitClusterMidMultiChild":
            clusters[np.sum(size_clusters[0:1])+int(size_clusters[1]/(1/(1-prop_split))):np.sum(size_clusters[0:2]),1] = 0
            clusters[np.sum(size_clusters[0:1])+int(size_clusters[1]/(1/(1-prop_split))):np.sum(size_clusters[0:2]),n_clusters] = 1
    #        print(clusters)
       #     print(np.dot(clusters, clusters.T))
            return (np.dot(clusters,clusters.T),clusters)
        elif "SplitClusterBot" in scenario:
            clusters[np.sum(size_clusters[0:(n_clusters-1)])+int(size_clusters[n_clusters-1]/(1/(1-prop_split))):np.sum(size_clusters[0:(n_clusters)]),(n_clusters-1)] = 0
            clusters[np.sum(size_clusters[0:n_clusters-1])+int(size_clusters[n_clusters-1]/(1/(1-prop_split))):np.sum(size_clusters[0:n_clusters]),n_clusters] = 1
    #        print(clusters)
    #        print(np.dot(clusters, clusters.T))
            return (np.dot(clusters,clusters.T),clusters)
    elif scenario == "MergeClusterBot":
        clusters = np.copy(t_clusters[:,:-1])
        clusters[np.sum(size_clusters[0:(n_clusters-1)]):np.sum(size_clusters[0:(n_clusters)]),n_clusters-2] = 1 #fix cluster 5 (originally cluster 6)
        clusters[np.sum(size_clusters[0:(n_clusters-2)]):np.sum(size_clusters[0:(n_clusters-1)]),n_clusters-2] = 0 #merge clusters 4 and 5 (from true phylogeny)
        clusters[np.sum(size_clusters[0:(n_clusters-2)]):np.sum(size_clusters[0:(n_clusters-1)]),n_clusters-3] = 1
    #    print("printing clusters:")
   #     print(clusters)
        return (np.dot(clusters,clusters.T),clusters)
    elif scenario == "MergeClusterMid&BotOneChild":
        clusters = np.copy(t_clusters[:,:-1])
        clusters[np.sum(size_clusters[0:5]):np.sum(size_clusters[0:6]), 2] = 1 #merge clusters 3 and 6
  #      print(clusters)
        return (np.dot(clusters,clusters.T),clusters)
    elif scenario == "MergeClusterMid&BotMultiChild":
        clusters = np.copy(t_clusters[:,:-1])
        clusters[np.sum(size_clusters[0:5]):np.sum(size_clusters[0:6]), 4] = 1 #fix cluster 5 (originally cluster 6)
        clusters[np.sum(size_clusters[0:4]):np.sum(size_clusters[0:5]), 1] = 1 #merge clusters 2 and 5 (from true phylogeny)
        clusters[np.sum(size_clusters[0:4]):np.sum(size_clusters[0:5]), 4] = 0
   #     print(clusters)
        return (np.dot(clusters,clusters.T),clusters)
    elif scenario == "MergeClusterTop&Mid":
        clusters = np.zeros((np.sum(size_clusters[0:6]),5))
        clusters[0:np.sum(size_clusters[0:2]), 0] = 1 #merged cluster from clusters 1 and 2 in true phylogeny
        clusters[np.sum(size_clusters[0:2]):np.sum(size_clusters[0:3]), 1] = 1
        clusters[np.sum(size_clusters[0:3]):np.sum(size_clusters[0:4]), 2] = 1
        clusters[np.sum(size_clusters[0:4]):np.sum(size_clusters[0:5]), 3] = 1
        clusters[np.sum(size_clusters[0:5]):np.sum(size_clusters[0:6]), 4] = 1
    #    print(clusters)
        return (np.dot(clusters,clusters.T),clusters)
    if 'Extra' in scenario:
        clusters = np.zeros((np.sum(size_clusters[0:n_clusters]),n_clusters+1))
        clusters[:,:-1] = np.copy(t_clusters)
    #    print(2.0/float(n_clusters))
       # rand_clusters = np.random.binomial(n=1, p=2.0/float(n_clusters), size=n_clusters)
       # print(rand_clusters)
        if "SmallExtra" in scenario:
            num_extra = int(extra_prop*np.sum(size_clusters[0:n_clusters]))

        elif "BigExtra" in scenario:
            num_extra = int(extra_prop*np.sum(size_clusters[0:n_clusters]))
        #    print("num extra big", num_extra, big_extra_prop*np.sum(size_clusters[0:n_clusters]))
     #   print("clusters_before",clusters)
        num_extra_clust = map(round,extra_prop*size_clusters)
      #  print "num_extra ", num_extra_clust
     #   print("num extra clust",num_extra_clust)
        for i in range(n_clusters):
      #      print(rand_clusters[i],i)
     #       if (rand_clusters[i]==1):

            if (i ==0):
                clusters[0:int(round(extra_prop*size_clusters[0])),i] = 0
                clusters[0:int(round(extra_prop*size_clusters[0])),n_clusters] = 1
            else:
                clusters[np.sum(size_clusters[0:i]):np.sum(size_clusters[0:i])+int(round(extra_prop*size_clusters[i])),i] = 0
                clusters[np.sum(size_clusters[0:i]):np.sum(size_clusters[0:i])+int(round(extra_prop*size_clusters[i])),n_clusters] = 1

            #     clusters[np.sum(size_clusters[0]):np.sum(size_clusters[0])+num_extra_clust,i] = 0
            #     clusters[np.sum(size_clusters[0]):np.sum(size_clusters[0])+num_extra_clust,n_clusters] = 1    
            # clusters[np.sum(size_clusters[0:i]):np.sum(size_clusters[0:i])+num_extra_clust,i] = 0
            # clusters[np.sum(size_clusters[0:i]):np.sum(size_clusters[0:i])+num_extra_clust,n_clusters] = 1
    #    print("clusters_after",clusters)
        return (np.dot(clusters,clusters.T),clusters)
    else:
        raise LookupError("Invalid scenario")
         
#if __name__ == '__main__':
#     test_truth_ccm, test_truth_clusters = get_ccm_nvar('Truth',size_clusters=[3,4,3],n_clusters=3)
#     test_one_ccm = get_ccm_nvar('OneCluster',size_clusters=[3,4,3],n_clusters=3)
#     print("truth",test_truth_ccm)
#     print("one",test_one_ccm)
#     sc = calculate2(test_truth_ccm,test_one_ccm,method='js_divergence')
#     print("sc",sc)

#, test_truth_ccm)
#
#get_ccm_nvar('SplitClusterMidOneChild',size_clusters=[3,4,4],n_clusters=3)
#get_ccm_nvar('SplitClusterMidMultiChild',size_clusters=[3,4,4],n_clusters=3)
#get_ccm_nvar('SplitClusterBot',size_clusters=[3,4,4],prop_split=0.5,n_clusters=3)
#get_ccm_nvar('MergeClusterBot',size_clusters=[2,2,3,5,3,2],n_clusters=6,prop_split=0.5)
#get_ccm_nvar('MergeClusterMid&BotOneChild',size_clusters=[2,2,3,5,3,2],n_clusters=6,prop_split=0.5)
#get_ccm_nvar('MergeClusterMid&BotMultiChild',size_clusters=[2,2,3,5,3,2],n_clusters=6,prop_split=0.5)
#test_clust,test_ccm = get_ccm_nvar('MergeClusterTop&Mid',size_clusters=[2,2,3,5,3,2],n_clusters=6,prop_split=0.5)
#    test_ccm,test_clust = get_ccm_nvar('SmallExtra',size_clusters=[6,6,6,6],n_clusters=4,prop_split=0.5,small_extra_prop=0.2)
#print(test_clust)
#test_rand, rand_ccm = closest_rand_reassign(test_clust,p_err=0.1)
#print(test_rand)