#! /gpfs/user/menghaowei/anaconda2/bin/python
# _*_ coding: UTF-8 _*_

import argparse
import sys
import pysam

# Version information START --------------------------------------------------
VERSION_INFO = \
    """
    Author: MENG Howard

    Version-01:
        2022-11-13
            count .bed like format table from BAM file, and run normalization related to the reference point

    E-Mail: meng_howard@126.com
    """
# Version information END ----------------------------------------------------


# ---------------------------------------------------------------------------->>>>>>>>>>
# function part
# ---------------------------------------------------------------------------->>>>>>>>>>
# FUN1
def query_multi_region_signal(
        signal_bam_obj,
        query_region_list,
        scale_method="CPM"):
    """
    INPUT:
        <signal_bam_obj>
            obj, pysam AlignmentFile obj

        <query_region_list>
            list, each item contains three info, chrome, start, end

        <scale_method>
            str, cpm, rpkm, raw

    RETURN:
        <signal_value_list>
            list, same length as input <query_region_list>

    """

    # count total
    total_count = 0
    for info in signal_bam_obj.get_index_statistics():
        total_count += info.mapped

    # query value
    signal_value_list = []

    for region_info in query_region_list:
        region_chr_name = region_info[0]
        region_start = int(region_info[1])
        region_end = int(region_info[2])

        region_count = 0
        for align in signal_bam_obj.fetch(contig=region_chr_name, start=region_start, end=region_end):
            region_count += 1

        scale_value = 0
        if scale_method == "CPM":
            scale_value = region_count / 1.0 / (total_count / 1e6)

        elif scale_method == "RPKM":
            region_length = region_end - region_start
            scale_value = region_count / 1.0 / (total_count / 1e6) / (region_length / 1e3)

        elif scale_method == "RAW":
            scale_value = region_count

        signal_value_list.append(round(scale_value, 6))

    return signal_value_list


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

    parser.add_argument("-b", "--input_bam",
                        help="Input bam file", required=True)

    parser.add_argument("-o", "--output",
                        help="Output default=stdout", default="stdout", type=str)

    parser.add_argument("--up_extend_length",
                        help="upstream region extend length, default=3000", default=3000, type=int)

    parser.add_argument("--down_extend_length",
                        help="downstream region extend length, default=3000", default=3000, type=int)

    parser.add_argument("--extend_binsize",
                        help="region extend binsize, default=100", default=100, type=int)

    parser.add_argument("--count_norm_method",
                        help="Can be RPKM, Raw or CPM, default=RPKM", default="RPKM", type=str)

    parser.add_argument("--verbose",
                        help="Larger number means out more log info, can be 0,1,2,3 default=3",
                        default=3, type=int)

    # ---------------------------------------------------------------------------->>>>>>>>>>
    # load the parameters
    # ---------------------------------------------------------------------------->>>>>>>>>>
    ARGS = parser.parse_args()

    # load bam
    bam_file_obj = pysam.AlignmentFile(ARGS.input_bam, "r")

    # load bed
    in_bed_file = open(ARGS.input, "r")

    # set output
    if ARGS.output == "stdout":
        output_file = sys.stdout
    else:
        output_file = open(ARGS.output, "w")

    # ---------------------------------------------------------------------------->>>>>>>>>>
    # count
    # ---------------------------------------------------------------------------->>>>>>>>>>
    all_info_list = []

    for region in in_bed_file:
        region_list = region.strip().split("\t")

        if region_list[5] == "+":
            ref_site = int(region_list[1])
            upstream_extend_len = ARGS.up_extend_length
            downstream_extend_len = ARGS.down_extend_length
        else:
            ref_site = int(region_list[2])
            upstream_extend_len = ARGS.down_extend_length
            downstream_extend_len = ARGS.up_extend_length

        chr_name = region_list[0]

        # make sub region
        upstream_region_list = []
        query_region_list = []
        downstream_region_list = []

        # ref site region
        ref_site_start = ref_site - int(ARGS.extend_binsize / 2)
        ref_site_end = ref_site + int(ARGS.extend_binsize / 2)

        # upstream
        for up_region_start in range(ref_site_start - upstream_extend_len, ref_site_start, ARGS.extend_binsize):
            up_region_end = min(up_region_start + ARGS.extend_binsize, ref_site_start)

            upstream_region_list.append(
                [chr_name, up_region_start, up_region_end]
            )

        # downstream
        for down_region_start in range(ref_site_end, ref_site_end + downstream_extend_len, ARGS.extend_binsize):
            down_region_end = min(down_region_start + ARGS.extend_binsize, ref_site_end + downstream_extend_len)

            downstream_region_list.append(
                [chr_name, down_region_start, down_region_end]
            )

        # query region
        query_region_list.append([chr_name, ref_site_start, ref_site_end])

        # get signal value
        up_signal_val_list = query_multi_region_signal(
            signal_bam_obj=bam_file_obj,
            query_region_list=upstream_region_list,
            scale_method=ARGS.count_norm_method
        )

        query_signal_val_list = query_multi_region_signal(
            signal_bam_obj=bam_file_obj,
            query_region_list=query_region_list,
            scale_method=ARGS.count_norm_method
        )

        down_signal_val_list = query_multi_region_signal(
            signal_bam_obj=bam_file_obj,
            query_region_list=downstream_region_list,
            scale_method=ARGS.count_norm_method
        )

        # merge signal
        if region_list[5] == "+":
            merge_val_list = up_signal_val_list + query_signal_val_list + down_signal_val_list

        elif region_list[5] == "-":
            merge_val_list = down_signal_val_list[::-1] + query_signal_val_list + up_signal_val_list[::-1]

        all_info_list.append(
            region_list + merge_val_list
        )

    # sort signal
    all_info_list_sort = sorted(all_info_list, key=lambda x: sum(x[6:]) / (len(x) - 6), reverse=True)

    # output header
    header_list = ["chrom", "tss", "tes", "gene_id", "exon_num", "strand"]
    value_header_list = [
        "col_%s" % i for i in range(1, len(all_info_list_sort[0]) - 6 + 1)
    ]

    out_header_list = header_list + value_header_list
    out_header_str = ",".join(out_header_list)
    output_file.write(out_header_str + "\n")

    # output value
    for signal_info in all_info_list_sort:
        out_str = ",".join(map(str, signal_info))
        output_file.write(out_str + "\n")

    # close file
    if ARGS.output != "stdout":
        output_file.close()
