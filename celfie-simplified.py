#!/usr/bin/env python
import argparse
import os
import sys
import subprocess
from io import BytesIO

assert sys.version_info >= (
    3,
    7,
    9,
), "This script requires Python 3.7.9 or greater for (subprocess.run features)."

import tqdm
import numpy as np
import pandas as pd

np.seterr(divide="ignore", invalid="ignore")


def add_pseudocounts(value, array, meth, meth_depths):
    """finds values of gamma where logll cannot be computed, adds pseudo-counts to make
    computation possible

    value: checks for a value that will prevent computation; either 0 or 1
    array: gamma array to check for inproper value
    meth: np array of methylation counts
    meth_depths: np array of total number of reads (meth counts + unmethylated counts)
    """

    axis0, axis1 = np.where(
        array == value  # find indices where value isn't able to be computed
    )

    # TODO: Fix mutation of inputs
    meth[axis0, axis1] += 1  # add one read to methylated counts
    meth_depths[axis0, axis1] += 2  # adds two reads to total counts


def expectation(gamma, alpha):
    """calculates the components needed for loglikelihood for each iteration of gamma and alpha

    gamma: np matrix of the estimated 'true' methylation proportions
    alpha: np matrix of estimated mixing proportions
    """

    alpha = alpha.T[:, np.newaxis, :]
    gamma = gamma[..., np.newaxis]

    p0 = (1.0 - gamma) * alpha
    p1 = gamma * alpha

    p0 /= np.nansum(p0, axis=0)[np.newaxis, ...]
    p1 /= np.nansum(p1, axis=0)[np.newaxis, ...]

    return p0, p1


def compute_log_likelihood(p0, p1, x_depths, x, y_depths, y, gamma, alpha):
    """calculates the log likelihood P(X, Z, Y | alpha, gamma)

    p0: probability that read is methylated
    p1: probability read is unmethylated
    x_depths: input read depths
    x: input methylated reads
    y_depths: reference matrix read depths
    y: reference methylated counts
    gamma: estimated true methylation proportions
    alpha: estimated mixing proportions
    """

    # Reshape arrays for faster computation
    alpha = alpha.T[:, np.newaxis, :]
    gamma = gamma[..., np.newaxis]

    y = y[..., np.newaxis]
    y_depths = y_depths[..., np.newaxis]

    x = x.T[np.newaxis, ...]
    x_depths = x_depths.T[np.newaxis, ...]

    log_likelihood = 0
    log_likelihood += np.sum((y + p1 * x) * np.log(gamma))
    log_likelihood += np.sum((y_depths - y + p0 * (x_depths - x)) * np.log(1.0 - gamma))
    log_likelihood += np.sum((p1 * x + (x_depths - x) * p0) * np.log(alpha))

    return log_likelihood


def maximization(p0, p1, x, x_depths, y, y_depths):
    """maximizes log-likelihood, calculated in the expectation step
    calculates new alpha and gamma given these new parameters

    p0: probability that read is methylated
    p1: probability read is unmethylated
    x_depths: input read depths
    x: input methylated reads
    y_depths: reference matrix read depths
    y: reference methylated counts
    """

    individuals = p0.shape[2]

    # initialize vector
    ones_vector = np.ones(shape=(y.shape[0]))
    new_alpha = np.zeros((x.shape[0], y.shape[0]))

    # in case of overflow or error, transform nans to 0 and inf to large float
    p0 = np.nan_to_num(p0)
    p1 = np.nan_to_num(p1)
    x = np.nan_to_num(x)
    x_depths = np.nan_to_num(x_depths)

    # break up calculation into two terms
    term0 = 0
    term1 = 0

    for n in range(individuals):

        new_alpha[n, :] = np.dot(p1[:, :, n], x[n, :]) + np.matmul(
            p0[:, :, n], (x_depths[n, :] - x[n, :])
        )

        term1 += p1[:, :, n] * (np.outer(ones_vector, x[n, :]))
        term0 += p0[:, :, n] * (np.outer(ones_vector, x_depths[n, :] - x[n, :]))

    gamma = (term1 + y) / (term0 + term1 + y_depths)  # calculate new gamma

    # check if gamma goes out of bounds, if so add psuedocounts to misbehaving y values.
    if (0 in gamma) or (1 in gamma):
        add_pseudocounts(1, gamma, y, y_depths)
        add_pseudocounts(0, gamma, y, y_depths)
        gamma = (term1 + y) / (term0 + term1 + y_depths)  # recalculate gamma

    # return alpha to be normalized to sum to 1
    normalized_new_alpha = new_alpha / np.sum(new_alpha, axis=1)[:, np.newaxis]
    return normalized_new_alpha, gamma


def expectation_maximization(
    x, x_depths, y, y_depths, num_iterations, convergence_criteria
):
    """take in the input cfdna matrices and the reference data and
    runs the EM for the specified number of iterations, or stops once the
    convergence_criteria is reached

    x: methylated cfDNA read counts
    x_depths: depth of cfDNA
    y: methylated reference counts
    y_depths: depth of cfDNA
    convergence_criteria: difference between alpha + gamma before stopping

    """

    # randomly intialize alpha for each iteration
    alpha = np.random.uniform(size=(x.shape[0], y.shape[0]))
    alpha /= np.sum(alpha, axis=1)[:, np.newaxis]  # make alpha sum to 1

    # begin by checking for instances where there are no counts for y or y_depths
    add_pseudocounts(1, np.nan_to_num(y / y_depths), y, y_depths)
    add_pseudocounts(0, np.nan_to_num(y / y_depths), y, y_depths)

    # intialize gamma to reference values
    gamma = y / y_depths

    # perform EM for a given number of iterations
    for _ in range(num_iterations):

        p0, p1 = expectation(gamma, alpha)
        a, g = maximization(p0, p1, x, x_depths, y, y_depths)

        # check convergence of alpha and gamma
        alpha_diff = np.mean(abs(a - alpha)) / np.mean(abs(alpha))
        gamma_diff = np.mean(abs(g - gamma)) / np.mean(abs(gamma))

        if (
            alpha_diff + gamma_diff < convergence_criteria
        ):  # if convergence criteria, break
            break
        # set current evaluation of alpha and gamma
        alpha = a
        gamma = g

    # print ll for random restarts
    log_likelihood = compute_log_likelihood(
        p0, p1, x_depths, x, y_depths, y, gamma, alpha
    )

    return alpha, gamma, log_likelihood


def define_arrays(input_bed, tim_matrix_bed, num_unk):
    """
    takes input data matrix- cfDNA and reference, and creates the arrays to run in EM. Adds
    specified number of unknowns to estimate

    sample: pandas dataframe of data (samples and reference). Assumes there is 3 columns (chrom, start, end)
    before the samples and before the reference
    num_samples: number of samples to deconvolve
    num_unk: number of unknowns to estimate
    """

    test = input_bed.iloc[:, 3:].values.T
    train = tim_matrix_bed.iloc[:, 3:].values.T

    x_df = test[::2, :]
    x_depths = test[1::2, :]

    y_df = train[::2, :]
    y_depths = train[1::2, :]

    # add N unknown components
    unknown = np.zeros((num_unk, y_depths.shape[1]))
    y_depths_unknown = np.append(y_depths, unknown, axis=0)
    y_unknown = np.append(y_df, unknown, axis=0)

    return (
        np.nan_to_num(x_df),
        np.nan_to_num(x_depths),
        np.nan_to_num(y_unknown),
        np.nan_to_num(y_depths_unknown),
    )


def write_output(output_file, output_matrix, header, index):
    """
    write estimated methylation proportions and tissue proportions as txt file

    output_file: outputfile name
    output_matrix: celfie estimate
    header: tissue names
    index: either number of cpgs or number of samples, depending on type of output
    written
    """

    output = pd.DataFrame(output_matrix)
    output.columns = header
    output.insert(
        0, "", index
    )  # insert either the sample names or cpg numbers as first col

    output.to_csv(output_file, sep="\t", index=False)


def validate_and_return_header_names(name_list):
    """
    Ensure an input list of sample / TIM matrix names is valid.

    Returns:
    sample_names: list of sample or TIM tissue names, e.g.:
       ['sample1', 'sample2', 'sample3', 'sample4', ...]
       or
       ['erythrocyte', 'lymphocyte', 'monocyte', ...]
    """
    assert name_list[0] == "chrom"

    # We expect an even # of columns (paired methylation & read data) after the chrom/start/end columns
    assert (len(name_list) - 3) % 2 == 0

    meth_names = [e[:-5] for e in name_list[3::2]]
    depth_names = [e[:-6] for e in name_list[4::2]]
    assert meth_names == depth_names

    return meth_names


def main(parsedargs):
    """
    Main function for running the EM algorithm.
    """
    os.makedirs(parsedargs.output_directory, exist_ok=True)
    print("Writing to: " + parsedargs.output_directory + "/")

    ## Loosely validate the TIM matrix
    print(f"Loading TIM matrix: {parsedargs.tim_matrix_bed}")
    tim_matrix_df = pd.read_csv(parsedargs.tim_matrix_bed, delim_whitespace=True, header=0)

    # .bed with tab or space after #
    if tim_matrix_df.columns[0] == "#":
        tim_cols_shifted = tim_matrix_df.columns[1:]
        tim_matrix_df = tim_matrix_df[tim_matrix_df.columns[:-1]]
        tim_matrix_df.columns = tim_cols_shifted
    tim_entry_names = validate_and_return_header_names(tim_matrix_df.columns)
    print(f"\tNumber of tissues in TIM matrix: {len(tim_entry_names)}")

    ## Parse input beds
    print(f"Loading data: {parsedargs.input_bed}")

    # Load the input bed, with our without a header
    with open(parsedargs.input_bed, "r", encoding="utf-8") as input_bed_fh:
        input_sample_first_line = input_bed_fh.readline()
        input_sample_has_header = input_sample_first_line.startswith("#")
        input_sample_second_line = input_bed_fh.readline()
        input_bed_number_of_sample_columns = len(
            input_sample_second_line.split("\t")[3:]
        )
    print(f"\tNumber of samples: {int(input_bed_number_of_sample_columns/2)}")

    USE_HEADER = 0
    if not input_sample_has_header:
        print(
            "\tNote: Input sample .bed file does not have a header. Samples will be labeled 'sample1', 'sample2', etc."
        )
        USE_HEADER = None

    # Validate the second line of input bed
    i_chr, i_start, i_end = input_sample_second_line.split("\t")[:3]
    i_size = int(i_end) - int(i_start)
    assert i_chr.startswith("chr")

    if parsedargs.skip_validation:
        print("WARNING: Skipping validation of input BED file loci sizes.")
        print(
            "Results may be odd if run on .bed files without individual CpG methylation count & coverage data."
        )
    else:
        # Validate the second line of input bed
        # We expect WGBS (or similar) data, with entries of size 1 or 2 basepairs.
        assert i_size in (1, 2)

    # Run bedtools map to sum features that overlap with the TIM matrix.
    # This is equivalent to:
    #  bedtools map -a tim_matrix.bed -b sample.bed -c 4,5 -null 0
    print("Running bedtools map...")
    print(
        "\tThis computes the sum of the features (both # methylated reads and # total reads) in the sample that overlap with the TIM matrix."
    )

    # Unfortuantely, we can't use pybedtools here, because it tries to be too smart and will mis-interpret
    # long pseudo-.bed files as .sam files, which results in nebulous errors (see: https://github.com/daler/pybedtools/issues/363)
    # Instead, we spawn a subprocess to run bedtools map:

    # We only want the first three columns (chrom, start, end)
    cut_command = f"cut -f1-3 {parsedargs.tim_matrix_bed}".split()
    cut_job = subprocess.run(cut_command, check=True, stdout=subprocess.PIPE)

    # Sum of all sample columns (4,5,6,...) for `bedtools map`
    COLUMNS_TO_SUM = str(
        list(range(4, 4 + input_bed_number_of_sample_columns))
    ).replace(" ", "")[1:-1]
    bedtools_command = (
        f"bedtools map -a stdin -b {parsedargs.input_bed} -c {COLUMNS_TO_SUM} -null 0".split()
    )
    bedtools_job = subprocess.run(
        bedtools_command, check=True, input=cut_job.stdout, capture_output=True
    )

    # Load the bedtools output (a .bed) as a pandas dataframe
    mapped_bed_df = pd.read_csv(
        BytesIO(bedtools_job.stdout), delim_whitespace=True, header=USE_HEADER
    )

    if USE_HEADER:
        sample_names = validate_and_return_header_names(mapped_bed_df.columns)
    else:
        sample_names = [
            "sample" + str(e) for e in range(1, input_bed_number_of_sample_columns)
        ]

    # Same number of rows in both (each row is one TIM matrix entry)
    assert mapped_bed_df.shape[0] == tim_matrix_df.shape[0]
    assert mapped_bed_df.shape[1] == len(sample_names * 2) + 3

    # make input arrays and add the specified number of unknowns
    x, x_depths, y, y_depths = define_arrays(
        mapped_bed_df, tim_matrix_df, parsedargs.unknowns
    )

    print("Starting computation...")

    # Run EM with the specified iterations and convergence criteria
    random_restarts = []

    for _ in tqdm.trange(parsedargs.random_restarts):
        alpha, gamma, ll = expectation_maximization(
            x, x_depths, y, y_depths, parsedargs.max_iterations, parsedargs.convergence
        )
        random_restarts.append((ll, alpha, gamma))

    # pick best random restart per replicate
    _, alpha_max, gamma_max = max(random_restarts)

    # get header for output files
    # N: samples here was: [nonpreg1, nonpreg2...]
    # tissues was: [dentricit, epithel...]

    # Save our results.
    write_output(
        f"{parsedargs.output_directory}/tissue_proportions.txt",
        alpha_max,
        tim_entry_names,
        sample_names,
    )
    write_output(
        f"{parsedargs.output_directory}/methylation_proportions.txt",
        gamma_max.T,
        tim_entry_names,
        list(range(len(gamma_max[1]))),
    )

    print("Done!")
    print(f"\tResult saved to: {os.getcwd()}/{parsedargs.output_directory}")


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="CelFiE - Cell-free DNA decomposition. CelFie estimates the cell type of origin proportions of a cell-free DNA sample."
    )
    parser.add_argument(
        "--input-bed", required=True, help="Your unknown sample(s) .bed file."
    )
    parser.add_argument(
        "--tim-matrix-bed",
        required=True,
        help="Your pre-trained tissue informative marker (TIM) matrix .bed.",
    )
    parser.add_argument(
        "--output-directory",
        required=True,
        help="Output directory. Any existing output files will be overwritten.",
    )
    parser.add_argument(
        "--skip-validation",
        default=False,
        action="store_true",
        help="Don't validate the input BED; this will run `bedtools map` regardless of the input BED contents. Use with caution.",
    )
    parser.add_argument(
        "--max-iterations",
        default=1000,
        type=int,
        help="How long the EM should iterate before stopping, unless convergence criteria is met. Default: 1000.",
    )
    parser.add_argument(
        "--unknowns",
        default=0,
        type=int,
        help="Number of unknown categories to be estimated along with the reference data. Default: 0.",
    )
    parser.add_argument(
        "--convergence",
        default=0.0001,
        type=float,
        help="Convergence criteria for EM. Default: 0.0001.",
    )
    parser.add_argument(
        "--random-restarts",
        default=10,
        type=int,
        help="Perform several random restarts and select the one with the highest log-likelihood. Default: 10.",
    )

    main(parser.parse_args())
