#! /gpfs/user/menghaowei/anaconda2/bin/python
# _*_ coding: UTF-8 _*_

import argparse
import sys
import os
import random
import string
import multiprocessing
import pysam
import logging

# Version information START --------------------------------------------------
VERSION_INFO = \
    """
    Author: MENG Howard

    Version-01:
        2021-11-05
            count .bed like format table from BAM file, and run normalization

    E-Mail: meng_howard@126.com
    """
# Version information END ----------------------------------------------------


# Learning Part START --------------------------------------------------------
LEARNING_PART = \
    """
    Design pipeline:

    Input:
        1. bed like format table

    Output:
        1. bed table
        2. additional count info
    """


# Learning Part END-----------------------------------------------------------


# ---------------------------------------------------------------------------->>>>>>>>>>
# function part
# ---------------------------------------------------------------------------->>>>>>>>>>
# FUN1
def check_input_bam_file(input_bam_filename):
    """
    HELP

        1. check exist

        2. check sort state

        3. check .bai index file

    RETURN
        0 alright

        Not 0:
            1 not exist
            2 not sorted
            3 index not exist
    """
    # check exist
    if not os.path.exists(input_bam_filename):
        return 1

    # check sort
    bam_file = pysam.AlignmentFile(input_bam_filename, "rb")
    if bam_file.header.get("HD").get("SO") != "coordinate":
        return 2

    # check bai index
    bai_filename_1 = input_bam_filename + ".bai"
    bai_filename_2 = os.path.splitext(input_bam_filename)[0] + ".bai"
    if (not os.path.exists(bai_filename_1)) and (not os.path.exists(bai_filename_2)):
        return 3

    return 0


# FUN2 get bam total align reads
def get_BAM_total_align_count(bam_filename):
    """
    INPUT:
        <bam_filename>
            BAM file path with .bai index file in a same dir

    RETURN
        <total_align_count>
            int, BAM file total align count

    """
    bam_file = pysam.AlignmentFile(bam_filename, "rb")

    total_align_count = 0

    for bam_index_info in bam_file.get_index_statistics():
        ref_name = bam_index_info[0]
        map_count = bam_index_info[1]

        total_align_count += map_count

    bam_file.close()

    return total_align_count


# FUN3
def query_region_align_count(
        bam_obj,
        region_chr,
        region_start,
        region_end,
        strand_info=".",
        MAPQ_cutoff=20,
        strand_select_method="b",
        extend_length=1):
    """
    INPUT
        <bam_obj>
            obj, pysam.AlignmentFile sorted by coordinate

        <region_chr>
            str, like chr1 chr2 ...

        <region_start> <region_end>
            int

        <strand_info>
            str, + - or .

        <MAPQ_cutoff>
            int, align with MAPQ lower than this will not be counted

        <strand_select_method>
            str,
                b: both strand
                f: foward strand
                r: reverse strand
                s: specific strand same to strand info

        <extend_length>
            int, region extend length

    RETURN:
        <query_align_count>
            int

    """

    query_align_count = 0

    for align in bam_obj.fetch(contig=region_chr,
                               start=region_start - extend_length,
                               end=region_end + extend_length):

        if align.mapq < MAPQ_cutoff:
            continue

        if strand_select_method == "s":
            if align.is_reverse and strand_info == "+":
                continue

            if not align.is_reverse and strand_info == "-":
                continue

        elif strand_select_method == "f":
            if align.is_reverse:
                continue

        elif strand_select_method == "r":
            if not align.is_reverse:
                continue

        elif strand_select_method == "b":
            pass

        else:
            raise ValueError("<strand_select_method> should be s, f, r or b!")

        query_align_count += 1

    return query_align_count


# FUN4 query with region file
def region_file_BAM_count(region_filename,
                          bam_filename,
                          MAPQ_cutoff=20,
                          strand_select_method="b",
                          extend_length=1,
                          norm_method="CPM",
                          log_verbose=3,
                          **args):
    """
    INPUT:
        <region_filename>
            str, bed-like format file

        <bam_filename>
            str, BAM file, indexed and sorted

        <norm_method>
            str, can be Raw, RPKM, CPM

        <**args>
           other settings can be generated from FUN <query_region_align_count>

    RETURN:
        <region_count_list>
            list, include align raw count
    """

    # --------------------------------------------------->>>>>>
    # log setting
    # --------------------------------------------------->>>>>>
    logging.basicConfig(level=(4 - log_verbose) * 10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    region_base_filename = os.path.basename(region_filename)
    bam_base_filename = os.path.basename(bam_filename)

    # --------------------------------------------------->>>>>>
    # calculate scale factor
    # --------------------------------------------------->>>>>>
    scale_factor = 1.0

    if norm_method == "CPM" or norm_method == "RPKM":
        bam_total_align_count = get_BAM_total_align_count(bam_filename)
        scale_factor = bam_total_align_count / 1.0 / 1e6

    # --------------------------------------------------->>>>>>
    # open file
    # --------------------------------------------------->>>>>>
    region_file = open(region_filename, "r")
    bam_file = pysam.AlignmentFile(bam_filename, "rb")

    # count
    region_count_list = []
    for index, line in enumerate(region_file):
        if index % 1000 == 0:
            logging.info("%s %s counting:%s" % (region_base_filename, bam_base_filename, index + 1))

        line_list = line.strip().split("\t")
        query_chr_name = line_list[0]
        query_start = int(line_list[1])
        query_end = int(line_list[2])

        if line_list[5] in ["+", "-", "."]:
            query_strand_info = line_list[5]
        else:
            query_strand_info = "."

        query_count = query_region_align_count(
            bam_obj=bam_file,
            region_chr=query_chr_name,
            region_start=query_start,
            region_end=query_end,
            strand_info=query_strand_info,
            MAPQ_cutoff=MAPQ_cutoff,
            strand_select_method=strand_select_method,
            extend_length=extend_length)

        if norm_method == "Raw":
            region_count_list.append(query_count)

        elif norm_method == "CPM":
            region_count_list.append(query_count / 1.0 / scale_factor)

        elif norm_method == "RPKM":
            region_len = query_end - query_start + 1 + 2 * extend_length
            region_count_list.append(query_count / 1.0 / scale_factor / region_len * 1e3)

    # close file
    region_file.close()
    bam_file.close()

    return region_count_list


# FUN5 query with region file
def multi_region_file_BAM_count(region_filename_list,
                                bam_filename,
                                thread_num=1,
                                MAPQ_cutoff=20,
                                strand_select_method="b",
                                extend_length=1,
                                norm_method="CPM",
                                log_verbose=3,
                                **args):
    """
    INPUT:
        <region_filename_list>
            list, each item is a filename

        <bam_filename>
            str, BAM file, indexed and sorted

        <threads_num>
            int, CPU threads used in this count part

        <**args>
           other settings can be generated from FUN <query_region_align_count>

    RETURN:
        <region_count_list>
            list, include align raw count
    """

    # --------------------------------------------------->>>>>>
    # log setting
    # --------------------------------------------------->>>>>>
    logging.basicConfig(level=(4 - log_verbose) * 10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    # --------------------------------------------------->>>>>>
    # run part
    # --------------------------------------------------->>>>>>
    logging.info("Starting to count: %s" % bam_filename)

    pool = multiprocessing.Pool(processes=thread_num)

    # record info
    input_BAM_count_result = []

    for query_region_filename in region_filename_list:
        input_BAM_count_result.append(
            pool.apply_async(
                func=region_file_BAM_count,
                args=(
                    query_region_filename,
                    bam_filename,
                    MAPQ_cutoff,
                    strand_select_method,
                    extend_length,
                    norm_method,
                    log_verbose,
                )
            )
        )

    pool.close()
    pool.join()

    # out pool result
    merge_total_count_list = []
    for temp_index, res in enumerate(input_BAM_count_result):
        run_res = res.get()
        merge_total_count_list += run_res

    logging.info("%s Calculation done!" % bam_filename)

    return merge_total_count_list


# FUN6 split file
def split_file_and_make_temp(input_filename, n_part=1, temp_dir=None, force_temp_dir=True, log_verbose=3):
    """
    """
    # --------------------------------------------------->>>>>>
    # log setting
    # --------------------------------------------------->>>>>>
    logging.basicConfig(level=(4 - log_verbose) * 10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    # --------------------------------------------------->>>>>>
    # set temp dir
    # --------------------------------------------------->>>>>>
    if temp_dir is None:
        temp_dir = os.path.abspath(os.path.dirname(input_filename))
    else:
        temp_dir = os.path.abspath(temp_dir)

    # add back slash
    if temp_dir[-1] != "/":
        temp_dir += "/"

    # temp dir check and create
    if not os.path.exists(temp_dir):
        if force_temp_dir:
            logging.warning("<temp_dir> setting is not exist \t %s " % temp_dir)
            logging.warning("<force_temp_dir> set as True, try to create temp dir \t %s" % temp_dir)

            try:
                os.makedirs(os.path.abspath(temp_dir))
            except:
                logging.warning("Temp dir creating error: \t %s" % temp_dir)
                logging.warning("set <temp_dir> as the same dir with <input_mpmat_filename>")
                temp_dir = os.path.abspath(os.path.dirname(input_filename))

        else:
            temp_dir = os.path.abspath(os.path.dirname(input_filename))
            logging.warning("<temp_dir> setting is not exist, set <temp_dir> as %s" % temp_dir)
    else:
        temp_dir = os.path.abspath(temp_dir)

    # --------------------------------------------------->>>>>>
    # get input basename
    # --------------------------------------------------->>>>>>
    input_file_basename = os.path.basename(input_filename)

    # --------------------------------------------------->>>>>>
    # create temp file name and open file
    # --------------------------------------------------->>>>>>
    logging.info("Making temp files...")
    temp_file_list = []
    temp_filename_list = []

    for index in range(n_part):
        # make filename
        temp_file_basename = "temp_" + \
                             input_file_basename + \
                             "." + \
                             str(index) + \
                             "." + \
                             "".join(random.sample(string.ascii_letters + string.digits, 16))
        temp_file_name = os.path.join(temp_dir, temp_file_basename)
        temp_filename_list.append(temp_file_name)

        # open file
        temp_file_list.append(open(temp_file_name, "w"))

    # --------------------------------------------------->>>>>>
    # counting input file line number
    # --------------------------------------------------->>>>>>
    logging.info("Counting input file...")

    total_input_line_num = 0
    with open(input_filename, "r") as input_file:
        for line in input_file:
            total_input_line_num += 1

        if total_input_line_num % n_part == 0:
            each_file_line_num = (total_input_line_num // n_part)
        else:
            each_file_line_num = (total_input_line_num // n_part) + 1

    logging.info("Done!")
    logging.info("input file total line count: %s" % total_input_line_num)

    # --------------------------------------------------->>>>>>
    # split temp files
    # --------------------------------------------------->>>>>>
    # write into output filename
    with open(input_filename, "r") as input_file:
        for index, line in enumerate(input_file):
            file_index = index // each_file_line_num
            temp_file_list[file_index].write(line)

    # close output filename
    input_file.close()
    [temp_file.close() for temp_file in temp_file_list]
    logging.info("Make temp files done!")

    return temp_filename_list


# FUN7 check file
def _check_file_exist(args):
    """
    INPUT:
        <args>
            obj, argparse obj

    RETURN:
        <check_all_state>
            bool, True means all files are exist, False means at least one of files can't pass file check step.

        <check_exist_dict>
            dict, each item contain 3 elements:
                1.input filename
                2.check state
                3. check reason
    """
    # make log
    logging.basicConfig(level=10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    # init list
    check_all_state = True
    input_filename_list = args.bam_file_list.split(",")

    for input_filename in input_filename_list:
        check_res = check_input_bam_file(os.path.abspath(input_filename))

        if check_res == 1:
            logging.error("BAM not exist! %s" % input_filename)
            check_all_state = False

        elif check_res == 2:
            logging.error("BAM not sorted by coordinate! %s" % input_filename)
            check_all_state = False

        elif check_res == 3:
            logging.error("BAM doesn't contain index file! %s" % input_filename)
            check_all_state = False

        else:
            logging.error("Checked OK! %s" % input_filename)

        sys.stderr.write("-" * 80 + "\n")

    # return part
    return check_all_state


# ---------------------------------------------------------------------------->>>>>>>>>>
#  main part
# ---------------------------------------------------------------------------->>>>>>>>>>
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Query BAM count info")

    # ==========================================================================================>>>>>
    # Input and output params
    # ==========================================================================================>>>>>
    parser.add_argument("-i", "--input",
                        help="Input bed", required=True)

    parser.add_argument("-o", "--output",
                        help="Output default=stdout", default="stdout", type=str)

    parser.add_argument("-b", "--bam_file_list",
                        help="BAM file list, like 1.bam,2.bam,3.bam...", required=True)

    parser.add_argument("-t", "--bam_file_tag",
                        help="BAM tag, like H3K4me3,H4K36me2,H3K9me3...", required=True)

    parser.add_argument("-p", "--thread_num",
                        help="thread num, default=1", default=1, type=int)

    parser.add_argument("-s", "--strand_selection",
                        help="strand selection method, can be f, r, b, or s, default=b", default="b", type=str)

    parser.add_argument("--MAPQ_cutoff",
                        help="MAPQ cutoff, default=20", default=20, type=int)

    parser.add_argument("--region_extend_len",
                        help="region extend length, default=1", default=1, type=int)

    parser.add_argument("--count_norm_method",
                        help="Can be RPKM, Raw or CPM, default=CPM", default="CPM", type=str)

    parser.add_argument("--verbose",
                        help="Larger number means out more log info, can be 0,1,2,3 default=3",
                        default=3, type=int)

    # ---------------------------------------------------------------------------->>>>>>>>>>
    # * load the parameters
    # ---------------------------------------------------------------------------->>>>>>>>>>
    ARGS = parser.parse_args()

    bam_filename_list = ARGS.bam_file_list.split(",")
    bam_tag_list = ARGS.bam_file_tag.split(",")

    if len(bam_filename_list) != len(bam_tag_list):
        raise ValueError("<bam_file_list> doesn't match <bam_file_tag>!")

    # ---------------------------------------------------------------------------->>>>>>>>>>
    # check file exist
    # ---------------------------------------------------------------------------->>>>>>>>>>
    sys.stderr.write("\n" + "-" * 80 + "\n")
    check_res = _check_file_exist(args=ARGS)
    if not check_res:
        raise IOError("Please make sure each BAM file is exist!")

    # ---------------------------------------------------------------------------->>>>>>>>>>
    # split region file
    # ---------------------------------------------------------------------------->>>>>>>>>>
    region_split_filename_list = split_file_and_make_temp(input_filename=ARGS.input,
                                                          n_part=ARGS.thread_num,
                                                          temp_dir=None,
                                                          force_temp_dir=True,
                                                          log_verbose=ARGS.verbose)
    sys.stderr.write("\n" + "-" * 80 + "\n")

    # ---------------------------------------------------------------------------->>>>>>>>>>
    # count BAM by order
    # ---------------------------------------------------------------------------->>>>>>>>>>
    final_count_dict = {}
    dict_key_list = []

    for bam_index, run_bam_filename in enumerate(bam_filename_list):
        bam_tag = bam_tag_list[bam_index]
        count_list = multi_region_file_BAM_count(region_filename_list=region_split_filename_list,
                                                 bam_filename=run_bam_filename,
                                                 thread_num=ARGS.thread_num,
                                                 MAPQ_cutoff=ARGS.MAPQ_cutoff,
                                                 strand_select_method=ARGS.strand_selection,
                                                 extend_length=ARGS.region_extend_len,
                                                 norm_method=ARGS.count_norm_method,
                                                 log_verbose=ARGS.verbose)

        dict_key = "%s_%s" % (bam_tag, bam_index)
        final_count_dict[dict_key] = count_list
        dict_key_list.append(dict_key)

    # ---------------------------------------------------------------------------->>>>>>>>>>
    # merge file
    # ---------------------------------------------------------------------------->>>>>>>>>>
    # check output file dir
    if ARGS.output == "stdout":
        output = sys.stdout
    else:
        output = open(ARGS.output, "w")

    # make header
    with open(ARGS.input, "r") as input_file:
        for line in input_file:
            line_list = line.strip().split("\t")
            header_list = ["col_%s" % (i+1) for i in range(len(line_list))]
            header_list += bam_tag_list
            output.write("\t".join(header_list) + "\n")
            break

    # output
    with open(ARGS.input, "r") as input_file:
        for index, line in enumerate(input_file):
            line_list = line.strip().split("\t")
            for dict_key in dict_key_list:
                line_list.append(final_count_dict[dict_key][index])
            output.write("\t".join(map(str, line_list)) + "\n")

    # ---------------------------------------------------------------------------->>>>>>>>>>
    # close and clean temp file
    # ---------------------------------------------------------------------------->>>>>>>>>>
    sys.stderr.write("-" * 80 + "\n")
    logging.info("Removing temp files...")

    for temp_filename in region_split_filename_list:
        if os.path.exists(temp_filename):
            if os.path.isfile(temp_filename):
                try:
                    os.remove(temp_filename)
                except:
                    logging.error("Removing error on \n\t%s" % temp_filename)

    if ARGS.output != "stdout":
        output.close()

    logging.info("Everything done!")
