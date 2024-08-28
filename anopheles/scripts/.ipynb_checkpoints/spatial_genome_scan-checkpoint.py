import pandas as pd, numpy as np
import allel, math, argparse
from scipy.spatial import ConvexHull, distance
from pyproj import Geod
from shapely.geometry import Polygon, polygon

# spatial genome scan

parser=argparse.ArgumentParser()
parser.add_argument('--vcf', help='path to annotated VCF')
parser.add_argument('--sample_data',  help='path to sample data TSV')
parser.add_argument('--out', help='output path')
args = parser.parse_args()

def polygon_area(coords):
    """
    calculate area in square KM
    input: list of [[lat, long]] coords of polygon axes
    """
    poly = Polygon(coords)
    poly = polygon.orient(poly)
    geod = Geod(ellps="WGS84")
    poly_area, poly_perimiter = geod.geometry_area_perimeter(poly)
    return poly_area

def load_data(vcf, sample_data):
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
    metadata = pd.read_csv(sample_data, sep='\t')
    samples = callset['samples'] # sample IDs
    # deal with vcf
    gt = allel.GenotypeArray(callset['calldata/GT']) # genotype array
    ANN_Annotation = callset['variants/ANN_Annotation'] # SNP annotation
    ANN_Annotation_Impact = callset['variants/ANN_Annotation_Impact'] # SNP impact
    position = callset['variants/POS'] # SNP position
    # reorder metadata to match VCF
    metadata = metadata.set_index('sampleID')
    metadata = metadata.loc[samples]
    return gt, metadata, ANN_Annotation, ANN_Annotation_Impact, position

def filter_alleles(gt, ANN_Annotation, ANN_Annotation_Impact, position):
    ac = gt.count_alleles()
    allele_filter = np.sum(ac[:,1:4], axis=1) > 1
    gt = gt[allele_filter]
    position = position[allele_filter]
    ANN_Annotation = ANN_Annotation[allele_filter]
    ANN_Annotation_Impact = ANN_Annotation_Impact[allele_filter]
    ac = ac[allele_filter]
    return gt, ANN_Annotation, ANN_Annotation_Impact, position, ac

def at_position(position_index):
    # pull out genotype information at genomic position
    position_annotation = ANN_Annotation[position_index]
    position_impact = ANN_Annotation_Impact[position_index]
    position_position = position[position_index]
    position_count = ac[position_index]
    return position_annotation, position_impact, position_position, position_count

def alt_allele(alt, position_index):
    """
    for alternate allele (either [1, 2, 3]),
    get frequency, annotation info, and spatial spread of carriers
    """
    # for alternate alle, get frequency and annotation
    allele_frequency = position_count[alt] / sum(position_count)
    allele_annotation = position_annotation[alt-1]
    allele_impact = position_impact[alt-1]
    # individuals in GT array with that alternate allele
    mask = np.sum(gt[position_index] == alt, axis=1) > 0
    carriers = metadata[mask]
    locs = np.unique(np.array(carriers[['x','y']]), axis=0) # unique locations of allele carriers
    locs_x = locs[:,0]
    locs_y = locs[:,1]
    area = spatial_spread(locs)
    return allele_frequency, allele_annotation, allele_impact, area
 
def spatial_spread(locs):
    """
    calculate the spatial spread of individuals carrying given allele in km2
    """
    if len(locs) > 2: # if possible,
        # make convex hull over locations
        hull = ConvexHull(locs)
        coords = list(zip(locs[hull.vertices,0],locs[hull.vertices,1]))
        area = polygon_area(coords) / 1000000
    else:
        area = 0
    return area

gt, metadata, ANN_Annotation, ANN_Annotation_Impact, position = load_data(args.vcf, 
                                                                            args.sample_data)
gt, ANN_Annotation, ANN_Annotation_Impact, position, ac = filter_alleles(gt, 
                                                                            ANN_Annotation, 
                                                                            ANN_Annotation_Impact, 
                                                                            position)

SNP_position = np.full(len(position)*3, np.nan)
SNP_alternate = np.full(len(position)*3, np.nan)
SNP_frequency = np.full(len(position)*3, np.nan)
SNP_annotation = np.full(len(position)*3, np.nan, dtype=object)
SNP_impact = np.full(len(position)*3, np.nan, dtype=object)
SNP_area = np.full(len(position)*3, np.nan)

array_place = 0
for index in range(len(position)):
    position_annotation, position_impact, position_position, position_count = at_position(index)
    for alt_index in [1, 2, 3]:
        allele_frequency, allele_annotation, allele_impact, area = alt_allele(alt_index, index)
        if allele_frequency > 0:
            SNP_position[array_place] = position_position
            SNP_alternate[array_place] = alt_index
            SNP_frequency[array_place] = allele_frequency
            SNP_annotation[array_place] = allele_annotation
            SNP_impact[array_place] = allele_impact
            SNP_area[array_place] = area
            
            array_place += 1


SNP_annotation = SNP_annotation[ ~ np.isnan(SNP_position) ]
SNP_impact = SNP_impact[ ~ np.isnan(SNP_position) ]
SNP_position = SNP_position[ ~ np.isnan(SNP_position) ]
SNP_alternate = SNP_alternate[ ~ np.isnan(SNP_alternate) ]
SNP_frequency = SNP_frequency[ ~ np.isnan(SNP_frequency) ]
SNP_area = SNP_area[ ~ np.isnan(SNP_area) ]

df = pd.DataFrame({'position':SNP_position,
                    'alternate':SNP_alternate,
                    'frequency':SNP_frequency,
                    'annotation':SNP_annotation,
                    'impact':SNP_impact,
                    'area':SNP_area})
df.to_csv(args.out+'_spatial_genome_scan.txt', sep='\t') 
