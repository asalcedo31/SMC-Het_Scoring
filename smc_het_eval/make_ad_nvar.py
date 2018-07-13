import numpy as np
from SMCScoring import *
from metric_behavior import *

np.set_printoptions(threshold='nan')
def run_case(scenario,nssm):
    truth = get_ad('Truth', size_clusters=nssm)
    new = get_ad_nvar(scenario,size_clusters=np.repeat(nssm,6))
    orig = get_ad(scenario, size_clusters=nssm)
    
    print("truth")
    print(truth)
    print("orig")
    print(orig)
    print('new')
    print(new)
    return new, orig

def test_same(scenario,nssm=2):
    print("Testing ", scenario)
    new, orig = run_case(scenario,nssm)
    if (np.array_equal(new,orig)):
        print('Pass')
        return('Pass')
    else: 
        print('Fail')
        return('Fail')
def get_ad_nvar(scenario,  t_ad=None, size_clusters=100, nssms=None, prop_split=0.5, lab=False, extra_prop = 4.0/12.0):
    '''Find the ancestry-descendant matrix for the given scenario

    Attributes:
    :param scenario: string representing the clustering scenario being evaluated
    :param t_ad: optional value for the true AD matrix, to avoid comupting it multiple times
    '''
    if t_ad is None and nssms is None:
        # t_ad =  np.zeros((np.sum(size_clusters),n_clusters))
        # for i in range(n_clusters):
        if (lab == True):

            labels = range(1,7)  
            print reduce(lambda x,y: np.concatenate(x,y), map(lambda x: np.repeat(labels[x],size_clusters[x]), range(6)))
            t_ad =  np.zeros((np.sum(size_clusters)+1,np.sum(size_clusters)+1))
            t_ad[np.sum(size_clusters),:np.sum(size_clusters)] = reducemap(lambda x: np.repeat(labels[x],size_clusters[x]), range(6))
            t_ad[0:np.sum(size_clusters),np.sum(size_clusters)] = map(lambda x: np.repeat(labels[x],size_clusters[x]), range(6))
        else:
            t_ad =  np.zeros((np.sum(size_clusters),np.sum(size_clusters)))
#        t_ad[size_clusters[0]:np.sum(size_clusters[0:2]),:] = 1
        t_ad[size_clusters[0]:np.sum(size_clusters[0:2]),np.sum(size_clusters[0:3]):np.sum(size_clusters[0:5])] = 1
        t_ad[0:size_clusters[0],size_clusters[0]:] = 1
        t_ad[np.sum(size_clusters[0:2]):np.sum(size_clusters[0:3]),np.sum(size_clusters[0:5]):np.sum(size_clusters[0:6])] = 1
       
       # t_ad[size_clusters:2*size_clusters, 3*size_clusters:5*size_clusters] = 2
       # t_ad[2*size_clusters:3*size_clusters, 5*size_clusters:6*size_clusters] = 3
    if scenario in ["Truth", "SplitClusterBotSame", "MergeClusterBot"]:
        return t_ad
    elif scenario == "SplitClusterBotDiff":
        ad = np.copy(t_ad)
    #    print(np.sum(size_clusters[0:5])+int(size_clusters[4]/(1/(1-prop_split))))
        ad[np.sum(size_clusters[0:5]):np.sum(size_clusters[0:5])+int(size_clusters[5]/(1/(1-prop_split))),np.sum(size_clusters[0:5])+int(size_clusters[5]/(1/(1-prop_split))):np.sum(size_clusters[0:6])] = 1.
        
        #ad[int(5*size_clusters):int(5.5*size_clusters),int(5.5*size_clusters):int(6*size_clusters)] = 1.
        return ad
    elif scenario == "SplitClusterMidOneChild":
        ad = np.copy(t_ad)
        ad[np.sum(size_clusters[0:2])+int(size_clusters[2]/(1/(1-prop_split))):np.sum(size_clusters[0:3]),np.sum(size_clusters[0:5]):np.sum(size_clusters[0:6])] = 0
        #ad[int(2.5*size_clusters):int(3*size_clusters),5*size_clusters:int(6*size_clusters)] = 0
        
        return ad
    elif scenario == "SplitClusterMidMultiChild":
        ad = np.copy(t_ad)
        ad[np.sum(size_clusters[0:1])+int(size_clusters[1]/(1/(1-prop_split))):np.sum(size_clusters[0:2]),np.sum(size_clusters[0:3]):np.sum(size_clusters[0:5])] = 0
      #  ad[int(1.5*size_clusters):int(2*size_clusters),3*size_clusters:int(5*size_clusters)] = 0
        
        return ad
    elif scenario == "MergeClusterMid&BotOneChild":
        ad = np.copy(t_ad)
        ad[np.sum(size_clusters[0:2]):np.sum(size_clusters[0:3]), np.sum(size_clusters[0:5]):np.sum(size_clusters[0:6])] = 0
        
      #  ad[int(2*size_clusters):int(3*size_clusters), 5*size_clusters:int(6*size_clusters)] = 0
        return ad
    elif scenario == "MergeClusterMid&BotMultiChild":
        ad = np.copy(t_ad)
        ad[size_clusters[0]:np.sum(size_clusters[0:2]), np.sum(size_clusters[0:4]):np.sum(size_clusters[0:5])] = 0
        ad[np.sum(size_clusters[0:4]):np.sum(size_clusters[0:5]), np.sum(size_clusters[0:3]):np.sum(size_clusters[0:4])] = 1
        # ad[int(size_clusters):int(2*size_clusters), int(4*size_clusters):int(5*size_clusters)] = 0
        # ad[int(4*size_clusters):int(5*size_clusters), int(3*size_clusters):int(4*size_clusters)] = 1
        
        return ad
    elif scenario == "MergeClusterTop&Mid":
        ad = np.zeros((np.sum(size_clusters),np.sum(size_clusters)))
        ad[0:np.sum(size_clusters[0:2]), np.sum(size_clusters[0:2]):] = 1
        ad[np.sum(size_clusters[0:2]):np.sum(size_clusters[0:3]), np.sum(size_clusters[0:5]):] = 1
             
        # ad = np.zeros((6*size_clusters,6*size_clusters))
        # ad[0:2*size_clusters, 2*size_clusters:] = 1
        # ad[2*size_clusters:3*size_clusters, 5*size_clusters:] = 1
        return ad
    elif scenario == "ParentIsSibling":
        ad = np.copy(t_ad)
        ad[np.sum(size_clusters[0:3]):np.sum(size_clusters[0:4]),np.sum(size_clusters[0:4]):np.sum(size_clusters[0:5])] = 1
            
        # ad[3*size_clusters:4*size_clusters,4*size_clusters:5*size_clusters] = 1
        return ad
    elif scenario == "ParentIsGrandparent":
        ad = np.copy(t_ad)
        ad[size_clusters[0]:np.sum(size_clusters[0:2]),np.sum(size_clusters[0:4]):np.sum(size_clusters[0:5])] = 0               
        #ad[size_clusters:2*size_clusters,4*size_clusters:5*size_clusters] = 0
        return ad
    elif scenario == "ParentIsAunt":
        ad = np.copy(t_ad)
        ad[size_clusters[0]:np.sum(size_clusters[0:2]),np.sum(size_clusters[0:4]):np.sum(size_clusters[0:5])] = 0
        ad[np.sum(size_clusters[0:2]):np.sum(size_clusters[0:3]),np.sum(size_clusters[0:4]):np.sum(size_clusters[0:5])] = 1    
        # ad[size_clusters:2*size_clusters,4*size_clusters:5*size_clusters] = 0
        # ad[2*size_clusters:3*size_clusters,4*size_clusters:5*size_clusters] = 1
        return ad
    elif scenario == "ParentIsCousin":
        ad = np.copy(t_ad)
        ad[size_clusters[0]:np.sum(size_clusters[0:2]),np.sum(size_clusters[0:4]):np.sum(size_clusters[0:5])] = 0
        ad[np.sum(size_clusters[0:2]):np.sum(size_clusters[0:3]),np.sum(size_clusters[0:4]):np.sum(size_clusters[0:5])] = 1
        ad[np.sum(size_clusters[0:5]):np.sum(size_clusters[0:6]),np.sum(size_clusters[0:4]):np.sum(size_clusters[0:5])] = 1
    
        # ad[size_clusters:2*size_clusters,4*size_clusters:5*size_clusters] = 0
        # ad[2*size_clusters:3*size_clusters,4*size_clusters:5*size_clusters] = 1
        # ad[5*size_clusters:6*size_clusters,4*size_clusters:5*size_clusters] = 1
        return ad
    elif scenario == "ParentIsSiblingWithChildren":
        ad = np.copy(t_ad)
        ad[np.sum(size_clusters[0:2]):np.sum(size_clusters[0:3]), range(size_clusters[0],np.sum(size_clusters[0:2]))+range(np.sum(size_clusters[0:3]),np.sum(size_clusters[0:5]))] = 1 #adjust cluster 3's ancestry
        # ad[2*size_clusters:3*size_clusters, range(size_clusters,2*size_clusters)+range(3*size_clusters,5*size_clusters)] = 1 #adjust cluster 3's ancestry
        return ad
    elif scenario == "ParentIsNieceWithChildren":
        ad = np.copy(t_ad)
        ad[np.sum(size_clusters[0:2]):np.sum(size_clusters[0:3]), range(size_clusters[0],np.sum(size_clusters[0:2]))+range(np.sum(size_clusters[0:3]),np.sum(size_clusters[0:5]))] = 1 #adjust cluster 3's ancestry
        ad[np.sum(size_clusters[0:5]):np.sum(size_clusters[0:6]), range(size_clusters[0],np.sum(size_clusters[0:2]))+range(np.sum(size_clusters[0:3]),np.sum(size_clusters[0:5]))] = 1 #adjust cluster 3's ancestry
        # ad[2*size_clusters:3*size_clusters, range(size_clusters,2*size_clusters)+range(3*size_clusters,5*size_clusters)] = 1 #adjust cluster 3's ancestry
        # ad[5*size_clusters:6*size_clusters, range(size_clusters,2*size_clusters)+range(3*size_clusters,5*size_clusters)] = 1 #adjust cluster 6's ancestry
        return ad
    elif scenario == "OneCluster":
        if nssms is None:
            return np.zeros(t_ad.shape)
        return np.zeros((nssms,nssms), dtype=np.int8)
    elif scenario == "NClusterOneLineage":
        # np.triu() returns a copy, this does the triu() in memory instead
#        if nssms is None:
#            return np.triu(np.ones(t_ad.shape))
        ad = np.ones(t_ad.shape, dtype=np.int8)
        for i in xrange(t_ad.shape[0]):
            for j in xrange(i + 1):
                ad[i, j] = 0
        return ad
    elif scenario == "NClusterTwoLineages":
        ad = np.triu(np.ones(t_ad.shape), k=1)
        ad[1:np.sum(size_clusters[0:3])+2,np.sum(size_clusters[0:3])+2:] = 0
        # ad[2:3*size_clusters+2,3*size_clusters+2:] = 0
        return ad
    elif scenario == "NClusterCorrectLineage":
        ad = np.triu(np.ones(t_ad.shape), k=1)
        ad[size_clusters[0]:np.sum(size_clusters[0:2]),range(np.sum(size_clusters[0:2]),np.sum(size_clusters[0:3]))+range(np.sum(size_clusters[0:5]),np.sum(size_clusters[0:6]))] = 0 # equivalent of cluster 2 from true AD matrix
        ad[np.sum(size_clusters[0:2]):np.sum(size_clusters[0:3]),np.sum(size_clusters[0:3]):np.sum(size_clusters[0:5])] = 0 # cluster 3 from true AD matrix
        ad[np.sum(size_clusters[0:3]):np.sum(size_clusters[0:4]),np.sum(size_clusters[0:4]):] = 0 # cluster 4 from true AD matrix
        ad[np.sum(size_clusters[0:4]):np.sum(size_clusters[0:5]),np.sum(size_clusters[0:5]):np.sum(size_clusters[0:6])] = 0 # cluster 5 from true AD matrix
   
        # ad[size_clusters:2*size_clusters,range(2*size_clusters,3*size_clusters)+range(5*size_clusters,6*size_clusters)] = 0 # equivalent of cluster 2 from true AD matrix
        # ad[2*size_clusters:3*size_clusters,3*size_clusters:5*size_clusters] = 0 # cluster 3 from true AD matrix
        # ad[3*size_clusters:4*size_clusters,4*size_clusters:] = 0 # cluster 4 from true AD matrix
        # ad[4*size_clusters:5*size_clusters,5*size_clusters:6*size_clusters] = 0 # cluster 5 from true AD matrix
        return ad
    if 'Extra' in scenario:
        ad = np.copy(t_ad)
        # print size_clusters,"\n"
        # map(lambda x: sys.stdout.write(str(-round(extra_prop*size_clusters[x]))+" "), range(6))
     #   num_extra = int(float(extra_prop)*np.sum(size_clusters)/size_clusters.shape[0])
      #  print "num_extra ", num_extra
        if scenario == "SmallExtraNewBot": #clusters 1, 3, and 6 are ancestors to the new extra cluster
            range_ones = [range(0,size_clusters[0]) ,range(np.sum(size_clusters[0:2]),np.sum(size_clusters[0:3])),range(np.sum(size_clusters[0:5]),np.sum(size_clusters[0:6])) ] # clusters 1, 3, and 6 will be ancestors to the new extra cluster
            if int(round(extra_prop*size_clusters[0])) >0 :
                #adjust cluster 1 for mutations lost from cluster 1
                range_ones[0] =  range(0+int(round(extra_prop*size_clusters[0])), np.sum(size_clusters[0]))
            if int(round(extra_prop*size_clusters[2])) >0 :
                #adjust cluster 3 for mutations lost from cluster 3
                range_ones[1] =  range(np.sum(size_clusters[0:2])+int(round(extra_prop*size_clusters[2])), np.sum(size_clusters[0:3]))
            if int(round(extra_prop*size_clusters[5])) >0 :
                #adjust cluster 6 for mutations lost from cluster 6
                range_ones[2] =  range(np.sum(size_clusters[0:5])+int(round(extra_prop*size_clusters[5])), np.sum(size_clusters[0:6]))

            range_ones = range_ones[0] + range_ones[1]+ range_ones[2]
         #   print "range\n",range_ones
            ad = np.copy(t_ad)
            if int(round(extra_prop*size_clusters[0])) >0 :
                ad[:,0:int(round(extra_prop*size_clusters[0]))] = 0
                ad[range_ones,0: int(round(extra_prop*size_clusters[0]))] = 1
                ad[0:int(round(extra_prop*size_clusters[0])),:] = 0
            for i in range(1,6):
                if int(round(extra_prop*size_clusters[i])) >0 :
                    ad[:,np.sum(size_clusters[0:i]):int(round(extra_prop*size_clusters[i]))] = 0
                    ad[range_ones,np.sum(size_clusters[0:i]): np.sum(size_clusters[0:i])+int(round(extra_prop*size_clusters[i]))] = 1
                    ad[np.sum(size_clusters[0:i]):np.sum(size_clusters[0:i])+int(round(extra_prop*size_clusters[i])),:] = 0
                # ad[:,size_clusters*i] = 0
                # ad[range(1,size_clusters)+range(size_clusters+1,2*size_clusters)+range(3*size_clusters+1,4*size_clusters),size_clusters*i] = 1
                # ad[size_clusters*i,:] = 0
         #   print(ad)
            return ad
        elif scenario == "SmallExtraCurBot":
            ad = np.copy(t_ad)
            range_ones = [range(0,size_clusters[0]) ,range(np.sum(size_clusters[0:2]),np.sum(size_clusters[0:3]))] # clusters 1 and 3 will be ancestors to the new extra cluster
            if int(round(extra_prop*size_clusters[0])) >0 :
                #adjust cluster 1 for mutations lost from cluster 1
                range_ones[0] =  range(0+int(round(extra_prop*size_clusters[0])), np.sum(size_clusters[0]))
            if int(round(extra_prop*size_clusters[2])) >0 :
                #adjust cluster 3 for mutations lost from cluster 3
                range_ones[1] =  range(np.sum(size_clusters[0:2])+int(round(extra_prop*size_clusters[2])), np.sum(size_clusters[0:3]))
            range_ones = range_ones[0] + range_ones[1]
            
            if int(round(extra_prop*size_clusters[0])) >0 :
                ad[:,0:int(round(extra_prop*size_clusters[0]))] = 0
                ad[range_ones,0:int(round(extra_prop*size_clusters[0]))] = 1
                ad[0:int(round(extra_prop*size_clusters[0])),:] = 0
            # print "range\n",range_ones
            for i in range(1,6):
                if int(round(extra_prop*size_clusters[i])) >0 :
                 #   print "start extra: ", np.sum(size_clusters[0:i]), "\n"
                #    print "end extra: ",  np.sum(size_clusters[0:i])+int(round(extra_prop*size_clusters[i])), "\n"
                    ad[:,np.sum(size_clusters[0:i]): np.sum(size_clusters[0:i])+int(round(extra_prop*size_clusters[i]))] = 0
                  #  print(range(int(round(extra_prop*size_clusters[1])),size_clusters[0])+range(np.sum(size_clusters[0:2])+int(round(extra_prop*size_clusters[1])),np.sum(size_clusters[0:3])))
                  #  range_ones= range_ones+
                    ad[range_ones,np.sum(size_clusters[0:i]):np.sum(size_clusters[0:i])+int(round(extra_prop*size_clusters[i])) ] = 1
             #       ad[range(int(round(extra_prop*size_clusters[1])),size_clusters[0])+range(np.sum(size_clusters[0:2])+int(round(extra_prop*size_clusters[1])),np.sum(size_clusters[0:3])),np.sum(size_clusters[0:i])] = 1
               
                    ad[np.sum(size_clusters[0:i]): np.sum(size_clusters[0:i])+int(round(extra_prop*size_clusters[i])),:] = 0
            
            # for i in range(0,6):
            #     ad[:,size_clusters*i] = 0
            #     ad[range(1,size_clusters)+range(2*size_clusters+1,3*size_clusters),size_clusters*i] = 1
            #     ad[size_clusters*i,:] = 0
            return ad
        elif scenario == "SmallExtraMid":
            ad = np.copy(t_ad)
            ad[:,0:int(round(extra_prop*size_clusters[0]))] = 0
            ad[range(int(round(extra_prop*size_clusters[0])),size_clusters[0]),0:int(round(extra_prop*size_clusters[0]))] = 1
            ad[0:int(round(extra_prop*size_clusters[0])),:] = 0
            for i in range(1,6):
                ad[:,np.sum(size_clusters[0:i]):np.sum(size_clusters[0:i])+int(round(extra_prop*size_clusters[i]))] = 0
                ad[range(int(round(extra_prop*size_clusters[0])),size_clusters[0]),np.sum(size_clusters[0:i]): np.sum(size_clusters[0:i])+int(round(extra_prop*size_clusters[i]))] = 1
                ad[np.sum(size_clusters[0:i]):np.sum(size_clusters[0:i])+int(round(extra_prop*size_clusters[i])),:] = 0
            # for i in range(0,6):
            #     ad[:,size_clusters*i] = 0
            #     ad[range(1,size_clusters),size_clusters*i] = 1
            #     ad[size_clusters*i,:] = 0
            return ad
        elif scenario == 'SmallExtraTop':
            ad = np.copy(t_ad)
            range_ones=[]
            
            if int(round(extra_prop*size_clusters[0])) >0 :
                range_ones = range_ones + range(int(round(extra_prop*size_clusters[0])),size_clusters[0])
            else:
                range_ones = range(0,size_clusters[0])
            for i in range(1,6):
                if int(round(extra_prop*size_clusters[i])) >0 :
          #          print "true", i
                    range_ones = range_ones + range(np.sum(size_clusters[0:i])+int(round(extra_prop*size_clusters[i])),np.sum(size_clusters[0:i+1]))
                else:
                    range_ones = range_ones + range(np.sum(size_clusters[0:i]),np.sum(size_clusters[0:i+1]))
             #   print range_ones
            if int(round(extra_prop*size_clusters[0])) >0 :
                ad[0:int(round(extra_prop*size_clusters[0])),range_ones] = 1

                ad[:,0:int(round(extra_prop*size_clusters[0]))] = 0
               # range(int(round(extra_prop*size_clusters[0])),size_clusters[0])+
               # range(size_clusters[0]+int(round(extra_prop*size_clusters[0])),np.sum(size_clusters[0:2]))+
               # range(np.sum(size_clusters[0:2])+int(round(extra_prop*size_clusters[1])),np.sum(size_clusters[0:3]))+
               # range(np.sum(size_clusters[0:3])+int(round(extra_prop*size_clusters[2])),np.sum(size_clusters[0:4]))+range(np.sum(size_clusters[0:4])+int(round(extra_prop*size_clusters[3])),np.sum(size_clusters[0:5]))+range(np.sum(size_clusters[0:5])+int(round(extra_prop*size_clusters[4])),np.sum(size_clusters[0:6]))] = 1
            # print "range", range_ones
            for i in range(1,6):
                if int(round(extra_prop*size_clusters[i])) >0 :
                    ad[:,np.sum(size_clusters[0:i]):np.sum(size_clusters[0:i]) +int(round(extra_prop*size_clusters[i]))] = 0
                    ad[np.sum(size_clusters[0:i]):np.sum(size_clusters[0:i])+int(round(extra_prop*size_clusters[i])), range_ones] = 1
                    
                   # range(int(round(extra_prop*size_clusters[0])),size_clusters[0])+
                   # range(size_clusters[0]+int(round(extra_prop*size_clusters[0])),np.sum(size_clusters[0:2]))+
                   # range(np.sum(size_clusters[0:2])+int(round(extra_prop*size_clusters[1])),np.sum(size_clusters[0:3]))+
                   # range(np.sum(size_clusters[0:3])+int(round(extra_prop*size_clusters[2])),np.sum(size_clusters[0:4]))+range(np.sum(size_clusters[0:4])+int(round(extra_prop*size_clusters[3])),np.sum(size_clusters[0:5]))+range(np.sum(size_clusters[0:5])+int(round(extra_prop*size_clusters[4])),np.sum(size_clusters[0:6]))] = 1
            # for i in range(0,6):
            #     ad[:,size_clusters*i] = 0
            #     ad[size_clusters*i,
            #        range(1,size_clusters)+
            #        range(size_clusters+1,2*size_clusters)+
            #        range(2*size_clusters+1,3*size_clusters)+
            #        range(3*size_clusters+1,4*size_clusters)+range(4*size_clusters+1,5*size_clusters)+range(5*size_clusters+1,6*size_clusters)] = 1
            return ad
        
    else:
        raise LookupError("Invalid scenario")
# print("Truth")
# t=get_ad_nvar('Truth',size_clusters=[3,2,3,3,2,2])
# t=get_ad_nvar('Truth',size_clusters=[3,2,0,3])
# print t
# test_same('SplitClusterBotDiff',nssm=3)
# test_same('SplitClusterBotSame',nssm=3)
# test_same('SplitClusterMidMultiChild',nssm=3)
# test_same('SplitClusterMidOneChild',nssm=3)

# # get_ad_nvar('SmallExtraNewBot',size_clusters=[2,2,2,2,2,2])
# # print("MergeClusterTop&Mid")
# print("Truth")
# truth = get_ad_nvar('Truth',size_clusters=[2,2,4,2,2,4],lab=False)
# #test = get_ad_nvar('SmallExtraTop',size_clusters=[4,3,2,2,3,2], extra_prop = float(6.0/12.0))
# test = get_ad_nvar('MergeClusterMid&BotMultiChild',size_clusters=[2,4,2,2,2,4], extra_prop = float(2.0/12.0))
# print("test")
# print(test)
# print("truth")
# print(truth)
# # # #get_ad_nvar('SplitClusterMidMultiChild',size_clusters=[4,4,4,4,4,4])
# orig = get_ad('SmallExtraNewBot', size_clusters=2)
# print("orig")
# print(orig)

# orig = get_ad('SmallExtraCurBot', size_clusters=2)
# print("orig")

# print(orig)
# get_ccm_nvar('SplitClusterMidOneChild',size_clusters=[3,4,4],n_clusters=3)