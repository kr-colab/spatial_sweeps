import pandas as pd, numpy as np
from matplotlib import pyplot as plt
import tskit, pyslim, msprime
from scipy.spatial import distance, ConvexHull
from scipy.spatial.distance import cdist
from shapely.geometry import Polygon, polygon
import os, subprocess
import itertools as it
import allel

def process_tree(treepath, logpath):
    # recapitate
    tree = tskit.load(treepath)
    ancestral_Ne=len(tree.tables.individuals)
    rtree = pyslim.recapitate(tree, 
                          recombination_rate=1e-8,
                          ancestral_Ne=ancestral_Ne,
                          random_seed=666)

    # simplify
    alive = pyslim.individuals_alive_at(rtree, 1)
    nodes = []
    for i in alive:
        nodes.extend(rtree.individual(i).nodes)
    stree = rtree.simplify(nodes)

    # add neutral variation
    next_id = pyslim.next_slim_mutation_id(stree)
    tree = msprime.sim_mutations(stree,
                             rate=1e-8,
                             model=msprime.SLiMMutationModel(type=0, next_id=next_id),
                             keep=True)
    return tree

def polygon_area(coords):
    """
    calculate area of polygon on SLiM landscape
    """
    poly = Polygon(coords)
    poly = polygon.orient(poly)
    return poly.area

def spatial_spread(locs):
    nlocs = len(locs)
    if nlocs == 2:
        distance = np.sqrt( (locs[0,0] - locs[1,0])**2 + (locs[0,1] - locs[1,1])**2 )
        area = 0
    elif nlocs > 2:
        hull = ConvexHull(locs)
        coords = list(zip(locs[hull.vertices,0],locs[hull.vertices,1]))
        area = polygon_area(coords)
        hullpoints = locs[hull.vertices,:]
        distance = np.max(cdist(hullpoints, hullpoints, metric='euclidean'))
    else:
        distance = 0
        area = 0
    return distance, area

def tree_frequency_area(tree):
    ind_count = len(tree.tables.individuals)
    ind_locations = tree.tables.individuals.location.reshape(ind_count, 3)[:,0:2]
    node_locations = np.repeat(ind_locations, 2, axis=0)

    # record frequencies and distances
    frequencies=np.repeat(np.nan, 100000000)
    areas = np.repeat(np.nan, 100000000)
    distances = np.repeat(np.nan, 100000000)
    positions = np.repeat(np.nan, 100000000)
    scoeffs = np.repeat(np.nan, 100000000)

    I = 0
    for var in tree.variants():
        allele = var.alleles[-1]
        frequency = var.frequencies()[allele]
        carrier_locations = node_locations[var.genotypes == 1]
        carrier_locations = np.unique(carrier_locations, axis=0)
        distance, area = spatial_spread(carrier_locations)
        selection_coeff = var.site.mutations[0].metadata['mutation_list'][0]['selection_coeff']
        position = var.site.position
        frequencies[I] = frequency
        areas[I] = area
        distances[I] = distance
        positions[I] = position
        scoeffs[I] = selection_coeff
        I += 1

    mask = ~np.isnan(frequencies)
    frequencies = frequencies[mask]
    areas = areas[mask]
    distances = distances[mask]
    positions = positions[mask]
    scoeffs = scoeffs[mask]

    df = pd.DataFrame({'freq':frequencies,
                  'area':areas,
                  'distance':distances,
                  'position':positions,
                  'scoef':scoeffs})
    return df


#sel_allele = scoeffs[scoeffs > 0][0]
#sel_freq = frequencies[scoeffs > 0][0]
#sel_area = areas[scoeffs > 0][0]

#print(f'{simID},{simTK},{sel_allele},{sel_freq},{sel_area}')
