import seaborn as sns
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
import multiprocessing


def compute_centers(closeness_df):
    closeness_df['centered_locus'] = round(
        (closeness_df['end'] - closeness_df['start']) / 2
        + closeness_df['start']).astype(int)

    closeness_df['close_centered_locus'] = round(
        (closeness_df['close_end'] - closeness_df['close_start']) / 2
        + closeness_df['close_start']).astype(int)

    closeness_df['centered_distance'] = closeness_df['close_centered_locus'] - closeness_df['centered_locus']
    return closeness_df


def compute_closeness_df(sorted_bed, sorted_bed_with):
    closeness_columns = ['chrom', 'start', 'end', 'name', 'score', 'strand',
                         'close_chrom', 'close_start', 'close_end',
                         'close_name', 'close_score', 'close_strand', 'distance']
    closeness = sorted_bed.closest(sorted_bed_with, D='a')
    closeness_df = closeness.to_dataframe()
    closeness_df.columns = closeness_columns
    closeness_df = compute_centers(closeness_df)
    return closeness_df


def build_closeness_analysis(bed, bed_with, assembly, conserve_chrom=False):
    closeness_columns = ['chrom', 'start', 'end', 'name', 'score', 'strand',
                         'close_chrom', 'close_start', 'close_end',
                         'close_name', 'close_score', 'close_strand', 'distance']

    shuffled_bed = bed.shuffle(genome=assembly, chrom=conserve_chrom).sort()
    shuffled_bed_with = bed_with.shuffle(genome=assembly, chrom=conserve_chrom).sort()

    sorted_bed = bed.sort()
    sorted_bed_with = bed_with.sort()

    # real on real
    closeness = sorted_bed.closest(sorted_bed_with, D='a')
    closeness_df = closeness.to_dataframe()
    closeness_df.columns = closeness_columns
    closeness_df = compute_centers(closeness_df)

    # genomic background1: real on shuffled
    shuffled_closeness = sorted_bed.closest(shuffled_bed_with, D='a')
    shuffled_closeness_df = shuffled_closeness.to_dataframe()
    shuffled_closeness_df.columns = closeness_columns
    shuffled_closeness_df = compute_centers(shuffled_closeness_df)

    # genomic background2: shuffled on shuffled
    totally_shuffled_closeness = shuffled_bed.closest(shuffled_bed_with, D='a')
    totally_shuffled_closeness_df = totally_shuffled_closeness.to_dataframe()
    totally_shuffled_closeness_df.columns = closeness_columns
    totally_shuffled_closeness_df = compute_centers(totally_shuffled_closeness_df)

    sns.set(style="whitegrid")
    sns.set_context("poster", font_scale=0.8, rc={"lines.linewidth": 1})
    g1 = sns.distplot(closeness_df[['centered_distance']], kde=False, bins=1000)
    g2 = sns.distplot(shuffled_closeness_df[['centered_distance']], kde=False, bins=1000)
    g3 = sns.distplot(totally_shuffled_closeness_df[['centered_distance']], kde=False, bins=10000)

    g3.set_xlim(-10000, 10000)


def compute_density(bed, bed_with, assembly, window, conserve_chrom=False):
    sorted_bed = bed.sort()
    sorted_bed_with = bed_with.sort()

    closeness_df = compute_closeness_df(sorted_bed, sorted_bed_with)

    density_hist = np.histogram(closeness_df['centered_distance'], bins='fd', density=True)
    density_df = pd.DataFrame(density_hist[0])
    density_df.columns = ['density']
    density_df['bin_edges'] = density_hist[1][1:].astype(float)
    density_df['abs_bin_edges'] = abs(density_df['bin_edges'])

    half_window = window / 2

    return sum(density_df.query("abs_bin_edges <= @half_window")['density'])


def compute_shuffled_density(bed, bed_with, assembly, window, conserve_chrom=False):
    sorted_bed = bed.sort()
    shuffled_bed = bed_with.shuffle(genome=assembly, chrom=conserve_chrom).sort()

    closeness_df = compute_closeness_df(sorted_bed, shuffled_bed)

    density_hist = np.histogram(closeness_df['centered_distance'], bins='fd', density=True)
    density_df = pd.DataFrame(density_hist[0])
    density_df.columns = ['density']
    density_df['bin_edges'] = density_hist[1][1:].astype(float)
    density_df['abs_bin_edges'] = abs(density_df['bin_edges'])

    half_window = window / 2

    return sum(density_df.query("abs_bin_edges <= @half_window")['density'])


def compute_density_enrichment(bed, bed_with, assembly, window, samples_num):
    print("Computing densities within window of", window, "centered on enhancers centers")

    num_cores = multiprocessing.cpu_count()

    densities = Parallel(n_jobs=num_cores)(delayed(compute_shuffled_density)
                                           (bed, bed_with, assembly, window, False)
                                           for i in range(0, samples_num))

    real_density = compute_density(bed, bed_with, assembly, window, False)

    print(real_density)

    shuffled_mean_density = np.mean(densities)
    shuffled_std = np.std(densities)

    print(shuffled_mean_density, shuffled_std)

    z_shuffled = (real_density - shuffled_mean_density) / shuffled_std

    return z_shuffled
