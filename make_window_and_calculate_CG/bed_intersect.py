#!/home/caogaoxiang/anaconda3/bin/python
# _*_ coding: UTF-8 _*_

import argparse
import os
import logging
import sys

# Version information START --------------------------------------------------
VERSION_INFO = """
Author: Gaoxiang Cao

Version-01:
    2023-03-21
        Find bed file intersect regions, support multi-region overlap

E-Mail: gaoxiang_cao@163.com
"""
# Version information END ----------------------------------------------------


# Learning Part START --------------------------------------------------------
LEARNING_PART = """

"""

# Learning Part END-----------------------------------------------------------


############################################################################################
# FUN Part
############################################################################################
def cmp_region(in_line_a, in_line_b, ref_order_dict):
    """
    INPUT:
    <in_line_a> <in_line_b>
     str, BED like format, tab seperated

     <ref_order_dic>
      dict, key is chr name, value is order index

    RETURN:
    int,
        -1, in_line_a is upstream of in_line_b
        0, in_line_a has overlap with in_line_b
        1, in_line_a is downstream of in_line_b
    """

    line_a_list = in_line_a.strip().split('\t')
    line_b_list = in_line_b.strip().split('\t')

    order_a = ref_order_dict.get(line_a_list[0])
    order_b = ref_order_dict.get(line_b_list[0])

    if order_a < order_b:
        return -1

    if order_a > order_b:
        return 1

    start_a = int(line_a_list[1])
    end_a = int(line_a_list[2])

    start_b = int(line_b_list[1])
    end_b = int(line_b_list[2])

    if end_a < start_b:
        return -1
    elif end_b < start_a:
        return 1
    else:
        return 0

############################################################################################
# main part
############################################################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Find bed file intersect regions, support multi-region overlap")

    parser.add_argument("-a", "--input_a",
                        help="BED file a", required=True)

    parser.add_argument("-b", "--input_b",
                        help="BED file a", default="out.pmat")

    parser.add_argument("-r", "--reference_index",
                        help="reference index file, can use samtools faidx to generate ", default="out.pmat")

    parser.add_argument("-o", "--output",
                        help="output file, default=stdout", default="stdout")

    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    # load args
    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    ARGS = parser.parse_args()

    # open file
    bed_a = open(ARGS.input_a, "r")
    bed_b = open(ARGS.input_b, "r")

    # output setting
    if ARGS.output == "stdout":
        out_file = sys.stdout
    else:
        out_file = open(ARGS.output, "w")

    # make ref ord dict
    ref_order_dict = {}

    ref_fai = open(ARGS.reference_index, "r")

    index = 0
    for line in ref_fai:
        line_list = line.strip().split("\t")
        ref_order_dict[line_list[0]] = index
        index += 1

    ref_fai.close()

    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    # run part
    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>

    line_a = bed_a.readline()
    line_b = bed_b.readline()

    temp_list = []
    new_line_a_state = True

    while line_a and line_b:

        if len(temp_list) > 0 and new_line_a_state:
            for temp_line_b in temp_list:
                temp_cmp_res = cmp_region(line_a, temp_line_b, ref_order_dict)

                if temp_cmp_res == 0:
                    out_str = line_a.strip() + "\t" + temp_line_b.strip()
                    out_file.write(out_str + "\n")

        cmp_res = cmp_region(line_a, line_b, ref_order_dict)

        if cmp_res == -1:
            line_a = bed_a.readline()
            new_line_a_state = True

        elif cmp_res == 1:
            new_line_a_state = False
            line_b = bed_b.readline()

        else:
            new_line_a_state = False
            out_str = line_a.strip() + "\t" + line_b.strip()
            out_file.write(out_str + "\n")

            temp_list.append(line_b)
            line_b = bed_b.readline()

    # test temp_list if empty
    file_a_end_state = False

    if line_a == '':
        file_a_end_state = True

    while not file_a_end_state:
        line_a = bed_a.readline()
        new_line_a_state = True

        if line_a == '':
            break
        if len(temp_list) > 0 and new_line_a_state:
            for temp_line_b in temp_list:
                temp_cmp_res = cmp_region(line_a, temp_line_b, ref_order_dict)

                if temp_cmp_res == 0:
                    out_str = line_a.strip() + "\t" + temp_line_b.strip()
                    out_file.write(out_str + "\n")
                    # Final edit date: 2023-03-21


