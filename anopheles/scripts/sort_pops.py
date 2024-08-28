import pandas as pd,  numpy as np

"""
split Anopheles gambiae metadata by Admixture results
"""

md = pd.read_csv('data/ag1000g_v3_gambiae.txt', sep='\t')
ax = pd.read_csv('data/gamb.k3.txt', sep=' ')

samples = ax.iloc[:,0]
admix = ax.iloc[:,2:] # first two columns are sampleID and country

def PopID(samples, admix, minimum_id=0.9):
    """
    assign samples to admixture populations
    """
    # mask mixed individuals
    unmix = np.max(admix, axis=1) >= minimum_id
    admix = admix[unmix]
    samples = samples[unmix]
    # assign to populations
    pops = np.argmax(admix, axis=1)
    return samples, pops

def PopDF(popID, samples, pops, metadata):
    """
    separate metadata for given admixture population
    """
    popsamples = samples[pops==popID]
    popmetadata = metadata[metadata['sampleID'].isin(popsamples)]
    return popmetadata
    
samples, pops = PopID(samples, admix)
for popID in np.unique(pops):
    popdf = PopDF(popID, samples, pops, md)
    popdf.to_csv(f'data/admixture_k{popID}_metadata.txt', sep='\t')
