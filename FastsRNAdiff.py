import argparse
import sys
import time
import math
import numpy as np
import scipy.stats as st
from pathlib import Path
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
from collections import defaultdict
from typing import Dict, List, Tuple, Any
from Visualization import Visualization

# Set the working directory for the software
def set_base_output(p: str = "./AnalysisRst"):
    global BASE_OUTPUT
    BASE_OUTPUT = Path(p).resolve()
    BASE_OUTPUT.mkdir(parents=True, exist_ok=True)

# Set output path
def out_path(*parts) -> Path:
    p = BASE_OUTPUT.joinpath(*parts)
    p.parent.mkdir(parents=True, exist_ok=True)
    return p

# If users provide Dicer enzyme information for small RNA clusters, further filtering
# can be performed to obtain more mature small RNA clusters.
def filterDicer(input_files: List[str],
                mut_type: str,
                outdir: Path) -> list[Any]:
    """
        filterDicer :
            Filter small RNA clusters based on whether DicerCall is present in the 20-24 nt range.
        Parameters:
            input_files: List[str]
                List of the file to be analyzed passed in by the user.
            mut_type: str
                Mutant type
            outdir: Path
                output dictionary
        Returns:
            List of the filtered small RNA cluster result files.
    """
    # Create the output directory;
    # if it does not exist, it will be created automatically.
    output_dir = outdir / "filterDicer"
    output_dir.mkdir(parents=True, exist_ok=True)
    filtered_files = []

    # Iterate through the list of input files and process each one individually.
    for cnt, file in enumerate(input_files, 1):
        try:
            src = Path(file)
            output_name = f"filterDicer_{mut_type}_{cnt}.txt"
            filtered_file = output_dir / output_name

            with open(src, "r") as fin, \
                    open(filtered_file, "w") as fon:
                # Read the header line and extract the indices of the required columns.
                header = fin.readline().strip().split()
                # There are three required attribute columns: #Locus    Reads(Total)    Rep-total
                # There is an optional DicerCall feature. In ShortStack analysis results, the DicerCall column is included.
                # If the DicerCall column is missing, the software will still function normally.
                indices = {
                    'locus': header.index('#Locus'),
                    'reads': header.index('Reads') if 'Reads' in header else header.index('Total') \
                        if 'Total' in header else -1,
                    'rep_total': header.index('Rep-total') if 'Rep-total' in header else header.index('Rep-Total') \
                        if 'Rep-Total' in header else -1,
                    'dicer': header.index('DicerCall') if 'DicerCall' in header else -1 #Non-essential attributes
                }
                # Write the header line to the output file.
                fon.write("#Locus\tReads\tRep-total\tDicerCall\n")

                # Declaration buffer
                buffer = []
                for line in fin:
                    parts = line.strip().split()
                    # Filter cluster
                    if not (parts[indices['dicer']] == 'N' and float(parts[indices['rep_total']]) >= 1):
                        line_data = [
                            parts[indices['locus']],
                            parts[indices['reads']],
                            parts[indices['rep_total']],
                            parts[indices['dicer']]
                        ]
                        # Temporarily store results in a buffer.
                        # When the buffer is full, write all data to the file at once to reduce write overhead.
                        buffer.append('\t'.join(line_data) + '\n')
                        if len(buffer) >= 1000:
                            fon.writelines(buffer)
                            buffer.clear()
                if buffer:
                    fon.writelines(buffer)

            filtered_files.append(str(filtered_file))
        except Exception as e:
            print(f"[Error] Precondition Failed : {file} -> {e}")
            continue

    return filtered_files

# Intersection check
def intersection(input_files: Dict[str,List[str]],
                 outdir: Path) -> str:
    """
        intersection:
            Perform cross-validation to extract common small RNA clusters.
            This function performs cross-validation on filtered small RNA
        analysis files to extract common small RNA clusters, ensuring
        stability and convenience for subsequent software analysis and
        visualization.
        Parameters:
            input_files: List[str]
                List of filtered small RNA analysis files
            outdir: Path
                output dictionary
        Returns:
            str: Path to the intersection results file
    """
    # Create dictionary
    output_file = outdir / "intersection.txt"
    # Check whether the input file list is empty
    if not input_files:
        raise ValueError("intersection: The input file list is empty.")

    filter_files = [item for file_dict in input_files.values()
                   for item in file_dict]

    try:
        # Read the first file to obtain the initial set of common sites.
        with open(filter_files[0], "r") as fin:
            next(fin)
            common_loci = {line.strip().split()[0] for line in fin}
    except Exception as e:
        print(f"[Error] Read Files Failed: {filter_files[0]} -> {e}")
        sys.exit(1)
    # Process the remaining files sequentially and update the set of common sites.
    for file in filter_files[1:]:
        try:
            with open(file, "r") as fin:
                next(fin)
                current_loci = {line.strip().split()[0] for line in fin}
                common_loci.intersection_update(current_loci)
                # If the public site is empty, terminate the loop early.
                if not common_loci:
                    break
        except Exception as e:
            print(f"[Error] Read Files Failed: {file} -> {e}")
            sys.exit(1)
    # Write the results to the output file.
    with output_file.open("w") as fos:
        fos.write("#Locus\n")
        fos.writelines(f"{locus}\n" for locus in sorted(common_loci))

    return str(output_file)

def pickUp_readCount(input_files: List[str],
                     base_file: str,
                     outdir: Path) -> List[str]:
    """
        pickUp_readCount:
            Pick up read counts from input filesFilter non-public
        small RNA cluster results, remove redundant columns, and
        extract key statistical test indicators(Locus,repTotal,
        subtract,total).
        Parameters:
            input_files: List[str]
                User input file with filtered small RNA cluster results。
            base_file:
                intersection file by Intersection function.
            outdir: Path
                output dictionary
        Returns:
            The result file after secondary filtering
    """
    with open(base_file, "r") as fin:
        next(fin)
        locus_set = {line.strip().split()[0] for line in fin}

    output_dir = outdir / "readCounts"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_files = []

    for input_file in input_files:
        try:
            path = Path(input_file)
            output_file = output_dir / f"readCounts_{path.name[12:]}"
            with open(input_file, "r") as fin, output_file.open("w") as fos:
                fos.write("#Locus\trepTotal\tsubtract\ttotal\n")
                buffer = []
                next(fin)
                for line in fin:
                    parts = line.strip().split()
                    locus = parts[0]
                    if locus in locus_set:
                        try:
                            # repTotal: Rounded to the nearest whole number
                            repTotal = int(round(float(parts[2])+0.001))
                        except (ValueError, IndexError):
                            repTotal = 0
                        try:
                            total = int(parts[1])
                        except (ValueError, IndexError):
                            total = 0
                        subtract = total - repTotal
                        buffer.append(f"{locus}\t{repTotal}\t{subtract}\t{total}\n")

                        if len(buffer) >= 1000:
                            fos.writelines(buffer)
                            buffer.clear()
                if buffer :
                    fos.writelines(buffer)
            output_files.append(str(output_file))
        except Exception as e:
            print(f"[Error] ReadCount process failed: {input_file} -> {e}")
            continue

    return output_files

# Chis Test
def chis_test(observed):
    """
        chis_test:
            Chi square test under large sample size. We constructed a
        2x2 contingency table for testing based on the Rep total and
        Total Mapped Reads indicators of wild-type and mutant types.
        Parameters:
            observed:
                A 2x2 contingency table.
                [   Reptotal_wild, Reptotal_mut
                    TotalMappedReads_wild, TotalMappedReads_mut]
        Returns:
            P-value , Ratio(Fold Change), avg_Mean
    """
    # Extract four observations from the contingency table
    a, b = float(observed[0][0]), float(observed[0][1])
    c, d = float(observed[1][0]), float(observed[1][1])
    # Perform chi square test
    chi2, p_value, dof, expected = st.chi2_contingency([[a,b],
                                                        [c,d]])
    # Calculate fold change(FC) and average mean
    ratio = ((b / (b + d)
              if (b + d) > 0 else 0)
             /
             (a / (a + c)
              if (a + c) > 0 else 1))
    mean = ((b / (b + d)
             if (b + d) > 0 else 0)
            +
            (a / (a + c)
             if (a + c) > 0 else 0)) * 5e5
    return {"p_value": p_value,
            "ratio": ratio,
            "mean": mean}

# Fisher test
def fisher_test(observed,
                havezero):
    """
        fisher_test:
            In cases where one or more of the read counts was five or less,
        Fisher's Exact Test was used instead.If the havezero is True, a value
        of 1 was substituted to avoid undefined ratios.
        Parameters:
            observed:
                A 2x2 contingency table.
                [   Reptotal_wild, Reptotal_mut
                    TotalMappedReads_wild, TotalMappedReads_mut]
        Returns:
            P-value , Ratio(Fold Change), avg_Mean
    """
    # Extract four observations from the contingency table
    a, b = float(observed[0][0]), float(observed[0][1])
    c, d = float(observed[1][0]), float(observed[1][1])
    # If havezero is True, add 1 to each value to avoid undefined ratios.
    _, p_value = st.fisher_exact([[a, b],
                                  [c, d]])
    if havezero:
        a, b, c, d = a+1, b+1, c+1, d+1
    # Perform Fisher's Exact Test
    # Calculate fold change and average mean
    ratio = ((b / (b + d)
              if (b + d) > 0 else 0  )
             /
             (a / (a + c)
              if (a + c) > 0 else 1))
    mean = ((b / (b + d)
             if (b + d) > 0 else 0)
            +
            (a / (a + c)
             if (a + c) > 0 else 0)) * 5e5
    return {"p_value": p_value,
            "ratio": ratio,
            "mean": mean}

def read_data_file(filename: Path) -> Dict[str, Tuple[int, int]]:
    """
        read_data_file:
            Read the data file and parse the site information
        and its corresponding two integer values in each line.
        Parameters:
            filename: Path
                Path to the data file.
        Returns:
            Dict[str, Tuple[int, int]]:
                Dictionary where keys are site names and values
            are tuples (rep_total, subtract).
    """
    data: Dict[str, Tuple[int, int]] = {}
    with filename.open("r") as f:
        next(f)
        for line in f:
            parts = line.strip().split()
            # Skip lines with insufficient number of columns
            if len(parts) < 3:
                continue
            locus = parts[0]
            try:
                rep_total = int(parts[1])
                subtract = int(parts[2])
                data[locus] = (rep_total, subtract)
            except ValueError:
                # Skip lines with invalid values
                continue
    return data

def process_pair(i: int,
                 j: int,
                 mut_type: str,
                 data1: Dict[str, Tuple[int, int]],
                 mapped1: int,
                 data2: Dict[str, Tuple[int, int]],
                 mapped2: int,
                 base_lines: List[List[str]],
                 output_dir: Path,
                 adjpvalue_threshold: float = 0.05):
    """
        process_pair:
            Process data from a pair of wild-type and mutant-type
        samples and perform statistical tests.
        Parameters:
            i : int
                Index of the currently processed wild-type sample.
            j : int
                Index of the currently processed mutant-type sample.
            data1 : Dict[str, Tuple[int, int]]
                Data dictionary for the first sample.
            mapped1 : int
                Total number of mappings for the first sample.
            data2 : Dict[str, Tuple[int, int]]
                Data dictionary for the second sample.
            mapped2 : int
                Total number of mappings for the second sample.
            base_lines : List[List[str]]
                List of base positions.
            output_dir : Path
                Path to the output directory.
            adjpvalue_threshold: float
                Significance level(Default: 0.05)
        Returns:
            str: Path string of the output result file.
    """
    output_file = output_dir / f"StaticsResult_wt_{i}_{mut_type}_{j}.txt"
    p_filtered_file = output_dir / f"StaticsResult_wt_{i}_{mut_type}_{j}_p_{adjpvalue_threshold}.txt"
    error_file = output_dir / f"check_wt{i}_{mut_type}_{j}.txt"

    with output_file.open("w") as out_f, \
         p_filtered_file.open("w") as p_f, \
         error_file.open("w") as err_f:

        out_f.write("#Locus\tP-value\tFC\tMean\n")
        p_f.write("#Locus\tP-value\tFC\tMean\n")

        for parts in base_lines:
            locus = parts[0]
            try:
                rep_total1, _ = data1.get(locus, (0, 0))
                rep_total2, _ = data2.get(locus, (0, 0))

                a, b = rep_total1, rep_total2
                c, d = mapped1 - a, mapped2 - b
                # Construct a statistical test matrix
                observed = np.array([[a, b], [c, d]], dtype=float)
                # Experience condition: Use Fisher for extremely small grids or containing 0
                havezero = (a == 0 or b == 0 or c == 0 or d == 0)
                # For small samples, use Fisher's exact test.
                if min(a, b, c, d) <= 5 or havezero:
                    result = fisher_test(observed, havezero)
                else:
                    result = chis_test(observed)

                out_line = f"{locus}\t{result['p_value']:.7g}\t{result['ratio']:.7g}\t{result['mean']:.7g}\n"
                out_f.write(out_line)

                if result['p_value'] < adjpvalue_threshold:
                    p_f.write(out_line)

            except Exception as e:
                err_f.write(f"Error processing locus {locus}: {str(e)}\n")

    return i, j, str(output_file), str(p_filtered_file)

def generate_statics_test(wt_files: List[str],
                          mut_files: List[str],
                          mut_type:str,
                          wt_mapped: List[int],
                          mut_mapped: List[int],
                          base_file: str,
                          outdir: Path,
                          adjpvalue_threshold: float = 0.05,
                          ):
    """
        generate_statics_test:
            Read the preprocessed file and construct a 2x2 contingency
        table for each small RNA cluster. Use chi square or Fisher's
        method for statistical testing.
        Parameters:
            wt_files: List[str],
                List of the wild type file
            mut_files: List[str],
                List of the mutant type file
            mut_type:str,
                Mutant type
            wt_mapped: List[int],
                List of the wild type mapped reads
            mut_mapped: List[int],
                List of the mutant type mapped reads
            base_file: str,
                Pass the Intersection txt
            outdir: Path,
            adjpvalue_threshold: float,
                Significance level(Default: 0.05)
        Returns:
            Tuple[List[str], List[str]]:
                List of raw statistical files List of files filtered at p< adjpvalue_threshold and unfiltered
    """
    n_wt = len(wt_files)
    n_mut=len(mut_files)
    # Store the unfiltered files
    results_matrix = [[None for _ in range(n_mut)] for __ in range(n_wt)]
    # Store the filtered files with p < 0.05
    p_filtered_matrix = [[None for _ in range(n_mut)] for __ in range(n_wt)]

    if len(wt_files) != len(wt_mapped):
        raise ValueError("Wild type error: The number of files does not match the number of mappings.")
    if len(mut_files) != len(mut_mapped):
        raise ValueError("Mutant type error: The number of files does not match the number of mappings.")

    output_dir = outdir / "StaticsRst"
    output_dir.mkdir(parents=True, exist_ok=True)
    # Read all wild-type and mutant data files in parallel.
    wt_data_list = [read_data_file(Path(f)) for f in wt_files]
    mut_data_list = [read_data_file(Path(f)) for f in mut_files]
    # Read basic site information
    with open(base_file, "r") as pf:
        next(pf)
        base_lines = [line.strip().split() for line in pf]

    tasks = []
    outputs: List[str] = []
    p_filtered_outputs: List[str] = []

    # Use a process pool to execute tasks concurrently
    with ProcessPoolExecutor(max_workers=8) as executor:
        # Parallel statistical testing between analysis files in wild-type and mutant strains
        for i, (data1, mapped1) in enumerate(zip(wt_data_list, wt_mapped),
                                             start=1):
            for j, (data2, mapped2) in enumerate(zip(mut_data_list, mut_mapped),
                                                 start=1):
                tasks.append(
                    # Submit statistical testing tasks to the executor
                    executor.submit(process_pair,
                                    i,
                                    j,
                                    mut_type,
                                    data1,
                                    mapped1,
                                    data2,
                                    mapped2,
                                    base_lines,
                                    output_dir,
                                    adjpvalue_threshold)
                )
        # Collect the results of completed tasks
        for fut in tqdm(as_completed(tasks),
                        total=len(tasks),
                        desc="Statistical Test",
                        position=0):

            i_res, j_res, outpath, p_filtered_path = fut.result()
            results_matrix[i_res - 1][j_res - 1] = outpath
            p_filtered_matrix[i_res - 1][j_res - 1] = p_filtered_path

    # Organized Results List
    for i in range(n_wt):
        for j in range(n_mut):
            outputs.append(results_matrix[i][j])
            p_filtered_outputs.append(p_filtered_matrix[i][j])

    return outputs, p_filtered_outputs

def estimate_adjp_values(PV, m=None, pi=1.0):
    """
        estimate_adjp_values:
            estimate q vlaues from a list of Pvalues
            this algorihm is taken from Storey, significance testing for genomic ...
        Parameters:
            m: number of tests, (if not len(PV)),
            pi: fraction of expected true null (1.0 is a conservative estimate)
    """
    if m is None:
        m = len(PV) * 1.0
    else:
        m *= 1.0
    lPV = len(PV)

    # 1. sort pvalues
    PV = PV.squeeze()
    IPV = PV.argsort()
    PV = PV[IPV]

    # 2. estimate lambda
    if pi is None:
        lrange = np.linspace(0.05, 0.95, max(lPV / 100.0, 10))
        pil = np.double((PV[:, np.newaxis] > lrange).sum(axis=0)) / lPV
        pilr = pil / (1.0 - lrange)
        # ok, I think for SNPs this is pretty useless, pi is close to 1!
        pi = 1.0
        # if there is something useful in there use the something close to 1
        if pilr[-1] < 1.0:
            pi = pilr[-1]

    # 3. initialise q values
    QV_ = pi * m / lPV * PV
    QV_[-1] = min(QV_[-1], 1.0)
    # 4. update estimate
    for i in range(lPV - 2, -1, -1):
        QV_[i] = min(pi * m * PV[i] / (i + 1.0), QV_[i + 1])
    # 5. invert sorting
    QV = np.zeros_like(PV)
    QV[IPV] = QV_
    return QV

def compute_adjpvalues_from_statics(statics_files: List[str],
                                    outdir: Path):
    """
        compute_adjpvalues_from_statics:
            Perform multiple hypothesis testing correction on
        statistical test result files, calculating the q-value
        (False Discovery Rate) corresponding to each p-value to
        reduce the number of false positives in a large number
        of concurrent statistical tests.
        Parameters:
            statics_files : list of str
                List of file paths for statistical test results or filtered statistical test results at p < 0.05.
            Each file should contain at least four columns(#Locus	P-value	FC	Mean).
        Returns:
            List of file paths for adjust P-value results.
    """

    out_dir = outdir / "adjPvalue"
    out_dir.mkdir(parents=True, exist_ok=True)
    outputFileList = []

    for sFile in statics_files:
        # Skip empty filter files (if p<0.05 yields no results)
        if Path(sFile).stat().st_size == 0:
            print(f"[Warning] Empty filtered file: {sFile}, skip adjp calculation")
            continue

        base_name = Path(sFile).stem
        qFile_all = out_dir / f"adjPvalue_{base_name[14:]}.txt"
        pFile = out_dir / f"oneColumn_parsingResult_{base_name[14:]}.txt"

        # Store four metrics
        loci, pvals, fcs, means = [], [], [], []
        '''
            Parse the input file to extract locus, p-value, 
        ratio, and mean information
        '''
        with open(sFile, "r") as sf:
            next(sf)
            for line in sf:
                parts = line.strip().split()
                if len(parts) < 4:
                    continue
                loci.append(parts[0])
                try:
                    pvals.append(float(parts[1]))
                except ValueError:
                    pvals.append(1.0)
                fcs.append(parts[2])
                means.append(parts[3])

        if not pvals:
            print(f"[Warning] No valid data in {sFile}, skip adjp calculation")
            continue

        pvals = np.array(pvals, dtype=float)
        adjPvalues = estimate_adjp_values(pvals)

        with open(pFile, "w") as pf:
            for p in pvals:
                pf.write(f"{p:.12g}\n")

        with open(qFile_all, "w") as of:
            of.write("#Locus\tadjPvalue\tFoldChange\tMean\n")
            for locus, qv, fc, mean_val in zip(loci, adjPvalues, fcs, means):
                line = f"{locus}\t{qv:.6g}\t{fc}\t{mean_val}\n"
                of.write(line)

        outputFileList.append(str(qFile_all))

    return outputFileList

def process_file(i: int,
                 filepath: str):
    """
        process_file:
                Read and parse a single qvalue file to extract locus-related
            attribute values.
        Parameters:
            i : int
                File index identifying the current input file being processed.
            filepath : str
                Path to the file to be processed.
        Returns:
            tuple: Contains three elements:
                - i:
                    File index.
                - mean_maps:
                    A list of three dictionaries storing the qvalue,
                    ratio, and mean values for each locus.
                - common_map:
                    A dictionary recording all loci that have appeared.
    """
    mean_maps = [dict(), dict(), dict()]
    common_map = dict()

    with open(filepath, "r") as f:
        next(f)
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 4:
                continue
            locus = parts[0]
            adjpvalue = float(parts[1])
            fc = float(parts[2])
            mean_val = float(parts[3])

            mean_maps[0][locus] = adjpvalue
            mean_maps[1][locus] = fc
            mean_maps[2][locus] = mean_val
            common_map[locus] = 1

    return i, mean_maps, common_map

def separate_groups(inputs: List[str],
                   mut_type: str,
                   output_base: Path,
                   adjpvalue_threshold) :
    """
        separate_groups:
                Perform aggregate calculations on all threshold files
            for a given mutation type, generating average values and
            classification result files.
        Parameters:
            inputs : List[str]
                List of input file paths.
            mut_type : str
                Name of the mutation type being processed.
            output_base : Path
                Base directory path for output files.
            adjpvalue_threshold :
                Significance level threshold (Default : 0.05)
        Returns:
            tuple: Contains two elements:
                - str: Path to the different type file lists.
                - str: Name of the mutation type.
            Output files:
                - Separate/mean/meanR_wt_vs_<mut_type>.txt:
                    Average statistics for all loci.
                - Separate/up/up_meanR_wt_vs_<mut_type>.txt:
                    Loci meeting the upregulation threshold (1-fold).
                - Separate/up/up2fold_meanR_wt_vs_<mut_type>.txt:
                    Loci meeting the significant upregulation threshold (2-fold).
                - Separate/down/down_meanR_wt_vs_<mut_type>.txt:
                    Loci meeting the downregulation threshold (1-fold).
                - Separate/down/down2fold_meanR_wt_vs_<mut_type>.txt:
                    Loci meeting the significant downregulation threshold (0.5-fold).
                - Separate/sig/nosig_meanR_wt_vs_<mut_type>.txt:
                    If the adjusted p-value exceeds the adjpvalue_threshold,
                    the small RNA cluster is classified as non-significant.
                - Separate/sig/sig_meanR_wt_vs_<mut_type>.txt:
                    If the adjusted p-value is less than the adjpvalue_threshold,
                    the small RNA cluster is classified as significant.
    """
    # Increase the threshold for extremely small quantities to avoid taking the logarithm of extremely small values.
    epsilon = 1e-10
    # Set the number of parallel separate_groups
    inner_jobs = min(len(inputs), 4)

    out_dir = output_base / "Separate"
    # 1x and 2x fold Up
    out_dir.mkdir(parents=True,
                  exist_ok=True)
    out_dir_up = out_dir / "up"
    out_dir_up.mkdir(parents=True,
                     exist_ok=True)
    # Calculation results for adjP, fold change (FC), and average abundance of small RNA clusters across all samples.
    out_dir_mean = out_dir / "mean"
    out_dir_mean.mkdir(parents=True,
                       exist_ok=True)
    # 1x and 2x fold Down
    out_dir_down = out_dir / "down"
    out_dir_down.mkdir(parents=True,
                       exist_ok=True)
    # significant or not significant
    out_dir_sig = out_dir / "sig"
    out_dir_sig.mkdir(parents=True,
                      exist_ok=True)

    output_files = {
        'mean': out_dir_mean / f"meanR_wt_vs_{mut_type}.txt",
        'up': out_dir_up / f"up_meanR_wt_vs_{mut_type}.txt",
        'up2fold': out_dir_up / f"up2fold_meanR_wt_vs_{mut_type}.txt",
        'down': out_dir_down / f"down_meanR_wt_vs_{mut_type}.txt",
        'down2fold': out_dir_down / f"down2fold_meanR_wt_vs_{mut_type}.txt",
        'nosig': out_dir_sig / f"nosig_meanR_wt_vs_{mut_type}.txt",
        'sig': out_dir_sig / f"sig_meanR_wt_vs_{mut_type}.txt"
    }

    lens = len(inputs)
    if lens == 0:
        return {'up':[str(output_files['up'])],
                'up2fold':[str(output_files['up2fold'])],
                'down':[str(output_files['down'])],
                'down2fold':[str(output_files['down2fold'])]
                }

    mean_maps = [defaultdict(dict) for _ in range(3)]
    common_maps = [defaultdict(int) for _ in range(lens)]
    global_common = defaultdict(int)
    # Use multiple processes to process all input files in parallel.
    with ProcessPoolExecutor(max_workers=inner_jobs) as executor:
        futures = [executor.submit(process_file, i, f) for i, f in enumerate(inputs)]
        for fut in tqdm(as_completed(futures), total=len(futures),
                        desc=f"Read {mut_type} file", position=1, leave=False):
            i, tmp_maps, tmp_common = fut.result()
            for k in range(3):
                mean_maps[k][i] = tmp_maps[k]
            common_maps[i] = tmp_common
            for locus in tmp_common:
                global_common[locus] = 1

    # Extract the locus present in all samples
    common_loci = [
        locus for locus in global_common
        if all(common_maps[i].get(locus, 0) for i in range(lens))
    ]

    with open(output_files['mean'], 'w') as f_mean, \
         open(output_files['up'], 'w') as f_up, \
         open(output_files['up2fold'], 'w') as f_up2, \
         open(output_files['down'], 'w') as f_down, \
         open(output_files['down2fold'], 'w') as f_down2, \
         open(output_files['nosig'], 'w') as f_nosig ,\
         open(output_files['sig'], 'w') as f_sig :

        f_mean.write(f"#Locus\tAvg_AdjP\tlog2(Avg_FC)\tAvg_Mean\n")
        f_up.write(f"#Locus\tAvg_AdjP\tlog2(Avg_FC)\tAvg_Mean\n")
        f_down.write(f"#Locus\tAvg_AdjP\tlog2(Avg_FC)\tAvg_Mean\n")
        f_up2.write(f"#Locus\tAvg_AdjP\tlog2(Avg_FC)\tAvg_Mean\n")
        f_down2.write(f"#Locus\tAvg_AdjP\tlog2(Avg_FC)\tAvg_Mean\n")
        f_sig.write(f"#Locus\tAvg_AdjP\tlog2(Avg_FC)\tAvg_Mean\n")
        f_nosig.write(f"#Locus\tAvg_AdjP\tlog2(Avg_FC)\tAvg_Mean\n")

        for locus in common_loci :
            adjPvalues = []
            fcs = []
            log2fcs = []
            abundances = []

            for i in range(lens):
                adjp_val = mean_maps[0][i].get(locus, 0)
                fc_val = mean_maps[1][i].get(locus, 1.0)
                log2fc_val = math.log2(float(fc_val) + epsilon)
                abundance_val = mean_maps[2][i].get(locus, 0)

                adjPvalues.append(float(adjp_val))
                fcs.append(float(fc_val))
                log2fcs.append(log2fc_val)
                abundances.append(float(abundance_val))
            # Screen samples with adjP < 0.05
            sig = all(q < adjpvalue_threshold for q in adjPvalues)
            # All genes with log2FC >= 1 within the same comparison group are labeled as the 2fold upregulated cluster.
            up2 = sig and all(log2fc >= 1 for log2fc in log2fcs)
            # All genes within the same comparison group with 1 > log2FC ≥ 0 are labeled as the 1fold upregulation cluster.
            up1 = sig and all(log2fc >= 0 for log2fc in log2fcs) and (not up2)
            # All genes with log2FC ≤ -1 within the same comparison group are labeled as the 2fold downregulated cluster.
            down2 = sig and all(log2fc <= -1 for log2fc in log2fcs)
            # All samples within the same comparison group where -1 < log2FC < 0 are labeled as the 1fold downregulated cluster.
            down1 = sig and all(log2fc < 0 for log2fc in log2fcs) and (not down2)

            # Calculate the average
            avg_adjP = sum(q for q in adjPvalues) / lens
            avg_fc = sum(fc for fc in fcs) / lens
            avg_mean = sum(mean for mean in abundances) / lens

            details = "\t".join(f"{adjPvalues[i]:.12g}:{log2fcs[i]:.7g}" for i in range(lens))
            f_mean.write(f"{locus}\t{avg_adjP:.6g}\t{math.log2(avg_fc+epsilon):.6g}\t{avg_mean:.6g}\t|{details}\n")

            if sig:
                f_sig.write(f"{locus}\t{avg_adjP:.6g}\t{math.log2(avg_fc+epsilon):.6g}\t{avg_mean:.6g}\t|{details}\n")
                if up1:
                    f_up.write(f"{locus}\t{avg_adjP:.6g}\t{math.log2(avg_fc+epsilon):.6g}\t{avg_mean:.6g}\t|{details}\n")
                if up2:
                    f_up2.write(f"{locus}\t{avg_adjP:.6g}\t{math.log2(avg_fc+epsilon):.6g}\t{avg_mean:.6g}\t|{details}\n")
                if down1:
                    f_down.write(f"{locus}\t{avg_adjP:.6g}\t{math.log2(avg_fc+epsilon):.6g}\t{avg_mean:.6g}\t|{details}\n")
                if down2:
                    f_down2.write(f"{locus}\t{avg_adjP:.6g}\t{math.log2(avg_fc+epsilon):.6g}\t{avg_mean:.6g}\t|{details}\n")
            else:
                f_nosig.write(f"{locus}\t{avg_adjP:.6g}\t{math.log2(avg_fc+epsilon):.6g}\t{avg_mean:.6g}\t|{details}\n")

    # List of files to be returned for subsequent visualization purposes
    result_files = dict()
    result_files['up'] = [output_files['up']]
    result_files['up2fold'] = [output_files['up2fold']]
    result_files['down'] = [output_files['down']]
    result_files['down2fold'] = [output_files['down2fold']]

    return result_files

def process_mut_type(mut_type,
                     mut_mapped,
                     pickUp_files_wt,
                     filtered_files,
                     iFile,
                     BASE_OUTPUT,
                     thresold,
                     args,
                     wt_mapped):
    """
        process_mut_type:
                Encapsulate the processing logic for individual mutation types as
            task functions for parallel execution.
        Parameters:
            mut_type:str
                Types of mutations to be processed
            mut_mapped:List[int]
                List of Total Mapped Reads for Mutant Samples
            pickUp_files_wt:List[Path]
                List of pickUp files
            filtered_files:List[Path]
                List of filtered files
            iFile:Path
                Intersection result file
            BASE_OUTPUT:
                Global working dictionary
            thresold:float
                Significant level (Default : 0.05)
            args,
            wt_mapped:
                List of Total Mapped Reads for Wild Samples
    """
    try:
        # pickUp_readCount Phase
        pickUp = pickUp_readCount(filtered_files[mut_type],
                                  iFile,
                                  BASE_OUTPUT)
        # Statistical Testing Phase
        statics, statics_p = generate_statics_test(
            pickUp_files_wt, pickUp, mut_type, wt_mapped, mut_mapped,
            iFile, BASE_OUTPUT, thresold
        )
        # Adjust P values Phase
        if args.strict:
            qFile = compute_adjpvalues_from_statics(statics, BASE_OUTPUT)
        else:
            qFile = compute_adjpvalues_from_statics(statics_p, BASE_OUTPUT)
        # Separate Groups Phase
        separate = separate_groups(qFile, mut_type, BASE_OUTPUT, thresold)
        # Visualization Phase
        visualizationDir = BASE_OUTPUT / "Visualization"
        visualizationDir.mkdir(parents=True, exist_ok=True)
        vFile = Visualization(separate, outdir=visualizationDir, mut_type=mut_type)

        # Return all results of a single mut_type
        return {
            "mut_type": mut_type,
            "pickUp": pickUp,
            "statics": statics,
            "statics_p": statics_p,
            "qFile": qFile,
            "separate": separate,
            "vFile": vFile,
            "success": True,  # Task marked as successful
            "error": None
        }
    except Exception as e:
        # Catch exceptions for individual tasks and return error messages
        return {
            "mut_type": mut_type,
            "pickUp": None,
            "statics": None,
            "statics_p": None,
            "qFile": None,
            "separate": None,
            "vFile": None,
            "success": False,
            "error": str(e)
        }

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='FastsRNAdiff',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage='FastsRNAdiff [-h] --wt-files WT_FILES [WT_FILES ...] --wt-mapped WT_MAPPED [WT_MAPPED ...]  '
              '--mut-type MUT_TYPE  --mut-files MUT_FILES [MUT_FILES ...] --mut-mapped MUT_MAPPED [MUT_MAPPED ...] '
              '[--mut-type MUT_TYPE  --mut-files MUT_FILES [MUT_FILES ...] --mut-mapped MUT_MAPPED [MUT_MAPPED ...] ...]'
              '[--output-dir OUTPUT_DIR] '
              '[--significance-level SIGNIFICANCE_LEVEL]',
        epilog='''
    Examples:
      FastsRNAdiff --wt-files wt1.txt wt2.txt --wt-mapped 1000000 1500000 --mut-type mutA --mut-files mut1.txt mut2.txt --mut-mapped 1200000 1300000
      FastsRNAdiff --wt-files control*.txt --wt-mapped 1000000 1500000 2000000 --mut-type mut1 --mut-files mut1_*.txt --mut-mapped 900000 1100000 --mut-type mut2 --mut-files mut2_*.txt --mut-mapped 800000 950000
            '''
    )

    # Mandatory: WT
    parser.add_argument('--wt-files',
                        nargs='+',
                        required=True,
                        help='List of wild-type file paths , use spaces to separate files.')
    parser.add_argument('--wt-mapped',
                        nargs='+',
                        type=int,
                        required=True,
                        help='The list of total mapping counts corresponding to wild-type files is in accordance with the file order.')

    # Mandatory: Mut
    parser.add_argument('--mut-type',
                        action='append',
                        help='Mutant type, can be used multiple times to specify multiple types.')
    parser.add_argument('--mut-files',
                        action='append',
                        nargs='+',
                        help='List of mutated file paths, with each -- mut type corresponding to a file list.')
    parser.add_argument('--mut-mapped',
                        action='append',
                        nargs='+',
                        type=int,
                        help='List of mutated mapping numbers, with each -- mut type corresponding to a mapping number list.')
    parser.add_argument('--output-dir',
                        default='OutputRst',
                        help='Specify the output path')
    parser.add_argument('--significance-level',
                        type=float,
                        default=0.05,
                        help='Specify the significance level')
    parser.add_argument('--strict',
                        action='store_true',
                        help='Use strict filtering (all samples must meet thresholds). Default: False (majority rule)')

    return parser.parse_args()

def validate_arguments(args):
    # The number of files in the same mutation sample is equal to the total number of mappings in the corresponding mutation sample.
    if len(args.wt_files) != len(args.wt_mapped):
        raise ValueError(f"The number of wild-type files does not match the number of mappings !")
    # The mutation type, mutation sample file set, and corresponding total mutation type mapping set must be provided.
    if not args.mut_type or not args.mut_files or not args.mut_mapped:
        raise ValueError("Mutational parameters must be provided: --mut-type / --mut-files / --mut-mapped")
    # The number of mutation types must match both the number of mutation sample files and the total mapped reads of mutation sample mappings.
    if len(args.mut_type) != len(args.mut_files) or len(args.mut_type) != len(args.mut_mapped):
        raise ValueError("The number of mutation types, files, and mappings does not match !")

    for i, (mut_files, mut_mapped) in enumerate(zip(args.mut_files, args.mut_mapped)):
        if len(mut_files) != len(mut_mapped):
            raise ValueError(f"The number of mutation type files does not match the number of mappings !")

def main():
    args = parse_arguments()
    validate_arguments(args)
    # Set base output dictionary
    set_base_output(args.output_dir)

    print("Parameter validation passed, analysis begins ...")
    start = time.time()

    thresold = args.significance_level
    # Preprocessing work
    filtered_files = dict()
    filtered_files["wt"] = filterDicer(args.wt_files, "wt", BASE_OUTPUT)
    if not filtered_files:
        raise ValueError("No valid wild-type files after filtering")

    for mut_idx, (mut_type, mut_files) in enumerate(zip(
        args.mut_type, args.mut_files
    ), 1):
        filtered_files[mut_type] = filterDicer(mut_files,mut_type,BASE_OUTPUT)

    iFile = intersection(filtered_files, BASE_OUTPUT)
    pickUp_files, statics_files, statics_p_files, qFiles, separateFiles, vFiles = dict(), dict(), dict(), dict(), dict(), dict()
    pickUp_files["wt"] = pickUp_readCount(filtered_files["wt"], iFile, BASE_OUTPUT)
    wt_mapped = args.wt_mapped

    task_params = [
        (mut_type, mut_mapped, pickUp_files["wt"], filtered_files, iFile, BASE_OUTPUT, thresold, args, wt_mapped)
        for mut_type, mut_mapped in zip(args.mut_type, args.mut_mapped)
    ]
    total_tasks = len(task_params)

    max_workers = os.cpu_count() or 4
    pbar = tqdm(total=total_tasks, desc="Handling mutation types ", unit="task", ncols=100)
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_param = {executor.submit(process_mut_type, *params): params for params in task_params}
        # Traverse completed tasks and collect results
        for future in as_completed(future_to_param):
            try:
                # Retrieve task results
                result = future.result()
                mut_type = result["mut_type"]
                # Populate the original dictionary with parallel results
                if result["success"]:
                    pickUp_files[mut_type] = result["pickUp"]
                    statics_files[mut_type] = result["statics"]
                    statics_p_files[mut_type] = result["statics_p"]
                    qFiles[mut_type] = result["qFile"]
                    separateFiles[mut_type] = result["separate"]
                    vFiles[mut_type] = result["vFile"]
                    tqdm.write(f"[Success] Mutant type {mut_type} has been completed!")
                else:
                    tqdm.write(f"[Error] Mutant type {mut_type} has been failed : {result['error']}")

                # Update progress bar (Complete one task, progress +1)
                pbar.update(1)
            except Exception as e:
                # Capture exceptions for individual tasks to prevent the entire system from crashing.
                print(f"Error occurred while processing mutation type: {e}")
    
    # Close the progress bar and output the final statistics.
    pbar.close()
    end = time.time()
    interval = end - start
    print(f"Analysis completed! Total number of tasks: {total_tasks}, Success Count: {len(pickUp_files)-1}, Number of failures: {total_tasks - (len(pickUp_files)-1)}. The results have been saved to {BASE_OUTPUT}. Total time: {interval:.3f}s.")

if __name__ == "__main__":
    main()