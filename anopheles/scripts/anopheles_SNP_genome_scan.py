import pandas as pd, numpy as np
import allel, argparse, re
from scipy.spatial import ConvexHull, distance
import malariagen_data

"""
analyze the frequency, spatial spread, and effect of SNPs in an Anopheles gene
input: 
    --transcript name of gene or transcript to analyze
    --chromosome name of chromosome being analyzed (2L, 2R, 3L, 3R, X)
    --output path to output csv
output:
    csv with columns 'mean_af': average allele frequency across sampling locations, 
                     'effect': imputed SNP effect, 
                     'impact': imputed SNP impact, 
                     'area': spatial spread of allele carriers, 
                     'position': genomic position
"""

parser=argparse.ArgumentParser()
parser.add_argument('--transcript', default=None, help='name of transcript to analyze')
parser.add_argument('--chromosome', default=None, help='chromosome being analyzed')
parser.add_argument('--out', help='path to output csv')
args = parser.parse_args()

# load Anopheles genome data
ag3 = malariagen_data.Ag3(
    "simplecache::gs://vo_agam_release",
    simplecache=dict(cache_storage="gcs_cache"),
)


# anopheles gambiae samples
df_samples = ag3.sample_metadata()
gambiae_samples = df_samples[df_samples.taxon == 'gambiae'].reset_index(drop=True)

def freq_area(transcript):
    # get snp calls for gene
    ds_snps = ag3.snp_calls(transcript)
    ds_snps = ds_snps.set_index(variants='variant_position', samples='sample_id')
    ds_snps = ds_snps.sel(samples = gambiae_samples.sample_id.to_numpy()) # only gambiae samples
    
    genotype = allel.GenotypeArray(ds_snps.call_genotype)
    genotype.index = ds_snps.variants
    
    snp_allele_freq = ag3.snp_allele_frequencies(transcript, sample_query="taxon == 'gambiae'", 
                           cohorts = 'admin2_year',
                           drop_invariant=False, effects=True)
    
    # mean SNP allele frequency across sampling locs
    cols = snp_allele_freq.columns
    reg = re.compile(r'^frq_')
    cols = list(filter(reg.search, cols))
    snp_allele_freq['mean_af'] = snp_allele_freq[cols].mean(axis=1)

    # new dataframe
    df = snp_allele_freq[['mean_af', 'effect', 'impact']]
    df['area'] = np.nan
    df['position'] = df.index.get_level_values('position')
    df = df.sort_index()
    positions = np.unique(df.index.get_level_values('position'))
    
    # loop thru snps
    for i in range(len(positions)):
        gt = genotype[i]
        position = positions[i]
        alleles = df.loc[args.chromosome, position]

        # get reference and alts for easy indexing
        ref_alleles = alleles.index.get_level_values('ref_allele')
        alt_alleles = alleles.index.get_level_values('alt_allele')

        # don't consider non-variant sites
        if sum(alleles.mean_af) <= 0:
            pass
        else:
            for p in range(len(alleles)):
                locs = None
                area = np.nan
                alt = p + 1 # alternate allele
                carriers = np.where(gt == alt)[0] # allele carriers

                if len(carriers) > 0:
                    # locations of allele carriers
                    locs = np.array(gambiae_samples.iloc[carriers][['longitude', 'latitude']])
                    locs = np.unique(locs, axis=0)
                    nlocs = len(locs)
                    if nlocs > 2:
                        area = ConvexHull(locs).area
                    elif nlocs == 2:
                        area = distance.euclidean(locs[0, :], locs[1, :])
                df.loc[(args.chromosome, position, ref_alleles[p], alt_alleles[p]), 'area'] = area
                
    return df.to_numpy()

vfreq_area = np.vectorize(freq_area)

if not args.transcript:
    # run on all genes in the chromosome!
    df_geneset = ag3.geneset()
    chrom_genes = df_geneset[(df_geneset.contig == args.chromosome) & (df_geneset.type == 'mRNA')].ID.to_numpy()
    df = vfreq_area(chrom_genes)
    df = np.concatenate(df)
    df = pd.DataFrame(df, columns = ['mean_af', 'effect', 'impact', 'area', 'position'])
    df.to_csv(args.out, sep='\t')
    
else:
    df = freq_area(args.transcript)
    df = pd.DataFrame(df, columns = ['mean_af', 'effect', 'impact', 'area', 'position'])
    df.to_csv(args.out, sep='\t')
    
   


"""
# snp calls for transcript
ds_snps = ag3.snp_calls(args.transcript)
ds_snps = ds_snps.set_index(variants="variant_position", samples="sample_id")
ds_snps = ds_snps.sel(samples=gambiae_samples.sample_id.to_numpy())

genotype = allel.GenotypeArray(ds_snps.call_genotype)
genotype.index = ds_snps.variants

snp_allele_freq = ag3.snp_allele_frequencies(args.transcript, sample_query="taxon == 'gambiae'", 
                           cohorts = 'admin2_year', ### NOT SURE WHAT THIS IS DOING BUT IT'S REQUIRED FOR THE FUNCTION
                           drop_invariant=False, effects=True)

# mean SNP allele frequency across sampling locs
cols = snp_allele_freq.columns
reg = re.compile(r'^frq_')
cols = list(filter(reg.search, cols))
snp_allele_freq['mean_af'] = snp_allele_freq[cols].mean(axis=1)

# new dataframe
df = snp_allele_freq[['mean_af', 'effect', 'impact']]
df['area'] = np.nan
df = df.sort_index()
positions = np.unique(df.index.get_level_values('position'))

# loop thru snps
for i in range(len(positions)):
    gt = genotype[i]
    position = positions[i]
    alleles = df.loc['2L', position]

    # get reference and alts for easy indexing
    ref_alleles = alleles.index.get_level_values('ref_allele')
    alt_alleles = alleles.index.get_level_values('alt_allele')

    # don't consider non-variant sites
    if sum(alleles.mean_af) <= 0:
        pass
    else:
        for p in range(len(alleles)):
            locs = None
            area = np.nan
            alt = p + 1 # alternate allele
            carriers = np.where(gt == alt)[0] # allele carriers

            if len(carriers) > 0:
                # locations of allele carriers
                locs = np.array(gambiae_samples.iloc[carriers][['longitude', 'latitude']])
                locs = np.unique(locs, axis=0)
                nlocs = len(locs)
                if nlocs > 2:
                    area = ConvexHull(locs).area
                elif nlocs == 2:
                    area = distance.euclidean(locs[0, :], locs[1, :])
            df.loc[(args.chromosome, position, ref_alleles[p], alt_alleles[p]), 'area'] = area
            
df.to_csv(args.out, sep='\t')
"""
