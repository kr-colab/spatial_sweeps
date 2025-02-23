import pandas as pd, numpy as np
import allel, math, argparse
from scipy.spatial import ConvexHull, distance
from scipy.spatial.distance import cdist 
from pyproj import Geod
from shapely.geometry import Polygon, polygon
from geopy.distance import distance as geodist

# spatial genome scan



parser=argparse.ArgumentParser()
parser.add_argument('--vcf', help='path to VCF')
parser.add_argument('--annotated', default=False, action='store_true', help='flag for VCFs annotated by SNPEff')
parser.add_argument('--sample_data',  help='path to sample data TSV')
parser.add_argument('--min_locs', default=3, type=int, help='minimum number of locations required to record a given allele. \
                                                    area recording works as follows for N given locations: \
                                                    if N>=3: calculate area of convex hull over sampling locations in square km \
                                                    if N==2: calculate area of a transect between two sampling locations. \
                                                             transect width is determined by --transect (default = 1km). \
                                                    if N==1: return default area covered by a given sampling site')
parser.add_argument('--transect', default=1, type=float, help='default transect width between two sampling sites')
parser.add_argument('--sample_area', default=1, type=float, help='default area (in km2) covered by a given sampling site')
parser.add_argument('--keep_locs', default=False, action='store_true', help='save locations of allele carriers? \
                                                                            if True, outputs seperate tsv appended with _allele_locations.txt \
                                                                            with the following format\
                                                                             ')
parser.add_argument('--out', help='output path')
args = parser.parse_args()
Min_locs=args.min_locs
Transect=args.transect
Sample_area=args.sample_area


def polygon_area(coords):
    """
    calculate area over a set of coordinates in square KM
    input: list of [[lat, long]] coords of polygon axes
    """
    poly = Polygon(coords)
    poly = polygon.orient(poly)
    geod = Geod(ellps="WGS84")
    poly_area, poly_perimiter = geod.geometry_area_perimeter(poly)
    return poly_area

def polygon_distance(coords):
    """
    calculate distance between two coords
    input: list of [[lat, long]] coords of polygon axes
    """

    geod = Geod(ellps="WGS84")
    az1, az2, dist = geod.inv(lons1=coords[0][1],
                              lats1=coords[0][0],
                              lons2=coords[1][1],
                              lats2=coords[1][0])
    return dist

def load_data(vcf, sample_data, annotated):
    # read metadata
    metadata = pd.read_csv(sample_data, sep='\t')
    if annotated:
        # read vcf data
        callset = allel.read_vcf(vcf, 
                                 fields=['samples',
                                         'calldata/GT',
                                         'variants/FILTER_PASS',
                                         'variants/POS',
                                         'ANN'],
                                 numbers={'ANN': 3},
                                 transformers=allel.ANNTransformer())
        # read metadata
        samples = callset['samples'] # sample IDs
        # deal with vcf
        gt = allel.GenotypeArray(callset['calldata/GT']) # genotype array
        ANN_Annotation = callset['variants/ANN_Annotation'] # SNP annotation
        ANN_Annotation_Impact = callset['variants/ANN_Annotation_Impact'] # SNP impact
        position = callset['variants/POS'] # SNP position
        
        # cut down genotypes to only included samples
        samples_list = list(samples)
        metadata = metadata[[i in samples_list for i in metadata['sampleID']]]
        samples_callset_index = [samples_list.index(s) for s in metadata['sampleID']]
        samples = np.array([samples[s] for s in samples_callset_index])
        metadata['callset_index'] = samples_callset_index
        indexes = metadata.callset_index.values.sort()    
        gt = gt.take(samples_callset_index, axis=1)
        return gt, metadata, ANN_Annotation, ANN_Annotation_Impact, position
    else:
        # read vcf data
        callset = allel.read_vcf(vcf, 
                                 fields=['samples',
                                         'calldata/GT',
                                         'variants/FILTER_PASS',
                                         'variants/POS'])
        samples = callset['samples'] # sample IDs
        # deal with vcf
        gt = allel.GenotypeArray(callset['calldata/GT']) # genotype array
        position = callset['variants/POS'] # SNP position
        
        # cut down genotypes to only included samples
        samples_list = list(samples)
        metadata = metadata[[i in samples_list for i in metadata['sampleID']]]
        samples_callset_index = [samples_list.index(s) for s in metadata['sampleID']]
        samples = np.array([samples[s] for s in samples_callset_index])
        metadata['callset_index'] = samples_callset_index
        indexes = metadata.callset_index.values.sort()    
        gt = gt.take(samples_callset_index, axis=1)
        return gt, metadata, position
    


def filter_alleles(gt, position, annotated, ANN_Annotation=None, ANN_Annotation_Impact=None):
    ac = gt.count_alleles()
    allele_filter = np.sum(ac[:,1:4], axis=1) > 1
    gt = gt[allele_filter]
    position = position[allele_filter]
    ac = ac[allele_filter]
    if annotated:
        ANN_Annotation = ANN_Annotation[allele_filter]
        ANN_Annotation_Impact = ANN_Annotation_Impact[allele_filter]
        return gt, ANN_Annotation, ANN_Annotation_Impact, position, ac
    else:
        return gt, position, ac

def at_position(position_index, position, ac, annotated, ANN_Annotation=None, ANN_Impact=None):
    # pull out genotype information at genomic position
    position_position = position[position_index]
    position_count = ac[position_index]
    if annotated:
        position_annotation = ANN_Annotation[position_index]
        position_impact = ANN_Annotation_Impact[position_index]
        return position_annotation, position_impact, position_position, position_count
    else:
        return position_position, position_count

def spatial_spread(locs, min_locs, transect, sample_area):
    """
    calculate the spatial spread of individuals carrying given allele in km2
    """
    nlocs = len(locs)
    if nlocs < min_locs:
        return np.nan, np.nan, np.nan
    elif nlocs == 1:
        return 1, sample_area, 0
    elif nlocs == 2:
        coords = list(zip(locs[:,0],locs[:,1]))
        area = polygon_distance(coords) / 1000
        
        # DEBUG!!
        #maxdist = max(cdist(coords, coords, lambda u, v: geodist(u, v).kilometers))
        maxdist=0
        
        return nlocs, area, maxdist
    else: 
        # make convex hull over locations
        hull = ConvexHull(locs)
        coords = list(zip(locs[hull.vertices,0],locs[hull.vertices,1]))
        area = polygon_area(coords) / 1000000
        
        # DEBUG!!
        #maxdist = max(cdist(coords, coords, lambda u, v: geodist(u, v).kilometers))
        maxdist=0


        return nlocs, area, maxdist


def alt_allele(alt, position_index, position_count, metadata, annotated,  position_annotation=None):
    """
    for alternate allele (either [1, 2, 3]),
    get frequency, annotation info, and spatial spread of carriers
    """
    # for alternate alle, get frequency and annotation
    if annotated:
        allele_frequency = position_count[alt] / (len(metadata) * 2)
        allele_impact = position_impact[alt-1]
        allele_annotation = position_annotation[alt-1]
    else:
        allele_frequency = position_count[alt] / (len(metadata) * 2)
    # individuals in GT array with that alternate allele
    mask = np.sum(gt[index] == alt, axis=1) > 1
    carriers = metadata[mask]
    locs = np.unique(np.array(carriers[['x','y']]), axis=0) # unique locations of allele carriers
    locs_x = locs[:,0]
    locs_y = locs[:,1]
    nlocs, area, maxdist = spatial_spread(locs, Min_locs, Transect, Sample_area) 
    if annotated:
        return allele_frequency, allele_annotation, allele_impact, area, nlocs, maxdist
    else:
        return allele_frequency, area, nlocs, maxdist
 

if args.annotated:
    gt, metadata, ANN_Annotation, ANN_Annotation_Impact, position = load_data(args.vcf, 
                                                                                args.sample_data,
                                                                                args.annotated)
    gt, ANN_Annotation, ANN_Annotation_Impact, position, ac = filter_alleles(gt, position, args.annotated,
                                                                                ANN_Annotation, 
                                                                                ANN_Annotation_Impact, 
                                                                                )

    SNP_position = np.full(len(position)*3, np.nan)
    SNP_alternate = np.full(len(position)*3, np.nan)
    SNP_frequency = np.full(len(position)*3, np.nan)
    SNP_annotation = np.full(len(position)*3, np.nan, dtype=object)
    SNP_impact = np.full(len(position)*3, np.nan, dtype=object)
    SNP_area = np.full(len(position)*3, np.nan)
    SNP_locs = np.full(len(position)*3, np.nan)
    SNP_maxdist = np.full(len(position)*3, np.nan)

    array_place = 0
    for index in range(len(position)):
        position_annotation, position_impact, position_position, position_count = at_position(index, position, ac, args.annotated, ANN_Annotation, ANN_Annotation_Impact)
        for alt_index in [1, 2, 3]:
            allele_frequency, allele_annotation, allele_impact, area, nlocs, maxdist = alt_allele(alt_index, index, position_count, metadata, args.annotated, position_annotation)
            if allele_frequency > 0:
                SNP_position[array_place] = position_position
                SNP_alternate[array_place] = alt_index
                SNP_frequency[array_place] = allele_frequency
                SNP_annotation[array_place] = allele_annotation
                SNP_impact[array_place] = allele_impact
                SNP_area[array_place] = area
                SNP_locs[array_place] = nlocs
                SNP_maxdist[array_place] = maxdist
                
                array_place += 1


    SNPmask = np.isnan(SNP_area)
    SNP_annotation = SNP_annotation[ ~ SNPmask ]
    SNP_impact = SNP_impact[ ~ SNPmask ]
    SNP_position = SNP_position[ ~ SNPmask ]
    SNP_alternate = SNP_alternate[ ~ SNPmask ]
    SNP_frequency = SNP_frequency[ ~ SNPmask ]
    SNP_area = SNP_area[ ~ SNPmask ]
    SNP_locs = SNP_locs[ ~ SNPmask ]
    SNP_maxdist = SNP_maxdist[ ~ SNPmask ]

    df = pd.DataFrame({'position':SNP_position,
                        'alternate':SNP_alternate,
                        'frequency':SNP_frequency,
                        'annotation':SNP_annotation,
                        'impact':SNP_impact,
                        'area':SNP_area})
    df.to_csv(args.out+'_spatial_genome_scan.txt', sep='\t') 

else:
    gt, metadata, position = load_data(args.vcf, 
                                                                                args.sample_data,
                                                                                args.annotated)
    gt, position, ac = filter_alleles(gt, position, args.annotated)

    SNP_position = np.full(len(position)*3, np.nan)
    SNP_alternate = np.full(len(position)*3, np.nan)
    SNP_frequency = np.full(len(position)*3, np.nan)
    SNP_area = np.full(len(position)*3, np.nan)
    SNP_locs = np.full(len(position)*3, np.nan)
    SNP_maxdist = np.full(len(position)*3, np.nan)

    array_place = 0
    for index in range(len(position)):
        position_position, position_count = at_position(index, position, ac, args.annotated)
        for alt_index in [1, 2, 3]:
            if alt_index >= len(position_count):
                pass
            else:
                allele_frequency, area, nlocs, maxdist = alt_allele(alt_index, index, position_count, metadata, args.annotated)
                if allele_frequency > 0:
                    SNP_position[array_place] = position_position
                    SNP_alternate[array_place] = alt_index
                    SNP_frequency[array_place] = allele_frequency
                    SNP_area[array_place] = area
                    SNP_locs[array_place] = nlocs
                    SNP_maxdist[array_place] = maxdist
                    
                    array_place += 1

    SNPmask = np.isnan(SNP_area)
    SNP_position = SNP_position[ ~ SNPmask ]
    SNP_alternate = SNP_alternate[ ~ SNPmask ]
    SNP_frequency = SNP_frequency[ ~ SNPmask ]
    SNP_area = SNP_area[ ~ SNPmask ]
    SNP_locs = SNP_locs[ ~ SNPmask ]
    SNP_maxdist = SNP_maxdist[ ~ SNPmask ]
    
    df = pd.DataFrame({'position':SNP_position,
                        'alternate':SNP_alternate,
                        'frequency':SNP_frequency,
                        'area':SNP_area})
    df.to_csv(args.out+'_spatial_genome_scan.txt', sep='\t') 
 