import logging
import os
from multiprocessing import Pool
from pysam import AlignmentFile
from bamtools._error import ParsingError


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
                               end=region_end + extend_length
                               ):
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
            raise ParsingError(
                "<strand_select_method> should be s, f, r or b!")

        query_align_count += 1

    return query_align_count


def get_BAM_total_align_count(bam_filename):
    """
    INPUT:
        <bam_filename>
            BAM file path with .bai index file in a same dir

    RETURN
        <total_align_count>
            int, BAM file total align count

    """
    bam_file = AlignmentFile(bam_filename, "rb")

    total_align_count = 0

    for bam_index_info in bam_file.get_index_statistics():
        ref_name = bam_index_info[0]
        map_count = bam_index_info[1]

        total_align_count += map_count

    bam_file.close()

    return total_align_count


def region_file_BAM_count(region_filename,
                          bam_filename,
                          MAPQ_cutoff=20,
                          strand_select_method="b",
                          extend_length=1,
                          norm_method="CPM"):
    """
     INPUT:
         <region_filename>
             str, bed-like format file

         <bam_filename>
             str, BAM file, indexed and sorted

         <norm_method>
             str, can be Raw, RPKM, CPM

     RETURN:
         <region_count_list>
             list, include align raw count
     """
    region_base_filename = os.path.basename(region_filename)
    bam_base_filename = os.path.basename(bam_filename)

    # -------------------------------------------------------->>>>>>
    # calculate scale factor
    # -------------------------------------------------------->>>>>>
    scale_factor = 1.0

    if norm_method == "CPM" or norm_method == "RPKM":
        bam_total_align_count = get_BAM_total_align_count(bam_filename)
        scale_factor = bam_total_align_count / 1.0 / 1e6

    # -------------------------------------------------------->>>>>>
    # open file
    # -------------------------------------------------------->>>>>>
    region_file = open(region_filename, "rt")
    bam_file = AlignmentFile(bam_filename, "rb")

    # count
    region_count_list = []

    for index, line in enumerate(region_file):
        if index % 1000 == 0:
            logging.info(
                f"{region_base_filename} {bam_base_filename} counting: {index+1}")

        line_list = line.rstrip().split("\t")
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
            region_count_list.append(
                query_count / 1.0 / scale_factor / region_len / 1e3)

    # close file
    region_file.close()
    bam_file.close()

    return region_count_list


def multi_region_file_BAM_count(
    region_filename_list,
    bam_filename,
    process_num,
    MAPQ_cutoff=20,
    strand_select_method="b",
    extend_length=1,
    norm_method="CPM",
):
    """
    INPUT:
        <region filename list>
            list, each file is a list

        <BAM filename>
            str, BAM file, indexed and sorted

        <thread num>
            int, CPU thread in this count part

    RETURN:
        <region count list>
            list, include align raw count
    """
    # -------------------------------------------------------->>>>>>
    # run part
    # -------------------------------------------------------->>>>>>
    logging.info("Starting to count: %s" % bam_filename)

    pool = Pool(processes=process_num)

    # record infor
    input_BAM_count_result = []

    for query_region_filename in region_filename_list:
        input_BAM_count_result.append(
            pool.apply_async(
                func=region_file_BAM_count,
                args=(query_region_filename,
                      bam_filename,
                      MAPQ_cutoff,
                      strand_select_method,
                      extend_length,
                      norm_method
                      )
            )
        )

    pool.close()
    pool.join()

    # output pool result
    merge_total_count_list = []

    for temp_index, res in enumerate(input_BAM_count_result):
        run_res = res.get()
        merge_total_count_list += run_res

    logging.info(f"{bam_filename} Calculation done!")

    return merge_total_count_list
