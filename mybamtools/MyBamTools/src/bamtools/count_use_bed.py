import gzip
import logging
import sys
import os
# import time
from bamtools import RemoveTempError
from bamtools import set_logging_level
from bamtools._spliter import split_file_and_make_temp
from bamtools._counter import multi_region_file_BAM_count


def main(
        input,
        output,
        bams,
        tags,
        method,
        process,
        mapq,
        strand,
        extend_length,
        temp_dir,
        verbose,
):
    # set log level
    set_logging_level(level=verbose)

    # parse params
    bam_filename_list = bams.split(",") if not isinstance(
        bams, tuple) else list(bams)
    bam_tag_list = tags.split(",") if not isinstance(
        tags, tuple) else list(tags)

    # print(bam_filename_list)
    # print(bam_tag_list)
    # print(type(bam_filename_list))
    # print(type(bam_tag_list))

    # split region file (bed)
    region_split_filename_list = split_file_and_make_temp(input_filename=input,
                                                          n_part=process,
                                                          temp_dir=temp_dir
                                                          )
    # print(f"region_split_filename_list:{region_split_filename_list}")

    # -------------------------------------------------------->>>>>>
    # Count BAM by order !
    # -------------------------------------------------------->>>>>>
    final_count_dict = {}
    tag_res_list = []

    for bam_index, run_bam_filename in enumerate(bam_filename_list):
        bam_tag = bam_filename_list[bam_index]
        count_list = multi_region_file_BAM_count(
            region_filename_list=region_split_filename_list,
            bam_filename=run_bam_filename,
            process_num=process,
            MAPQ_cutoff=mapq,
            strand_select_method=strand,
            extend_length=extend_length,
            norm_method=method)

        dict_key = f"{bam_tag}_{bam_index}"
        # print(dict_key)
        final_count_dict[dict_key] = count_list
        # print(f"fina_count_dict {final_count_dict}")
        tag_res_list.append(dict_key)
        # print(f"tag_res_list {tag_res_list}")
    # print(f"fina_count_dict {final_count_dict}")

    # -------------------------------------------------------->>>>>>
    # merge file
    # -------------------------------------------------------->>>>>>
    # check output file dir
    if not output:
        output_file = sys.stdout

    else:
        output_file = open(
            output,
            "wt") if ".gz" not in output else gzip.open(
            output,
            "wt")

    input_file = open(
        input,
        "rt") if not".gz" in input else gzip.open(
        input,
        "rt")
    line = input_file.readline().strip()
    input_file.close()

    line_list = line.strip().split("\t")
    header_list = [f"col_{i+1}" for i in range(len(line_list))]
    header_list += [f"{tag}_{method}" for tag in bam_tag_list]
    # print(header_list)

    output_file.write("\t".join(header_list) + "\n")

    # output
    input_file = open(
        input,
        "rt") if ".gz" not in input else gzip.open(
        input,
        "rt")
    for index, line in enumerate(input_file):
        line_list = line.strip().split("\t")
        # print(line_list)

        for dict_key in tag_res_list:
            # print(dict_key)
            # print(tag_res_list)
            line_list.append(final_count_dict[dict_key][index])
            # print(final_count_dict[dict_key][index])
        output_file.write("\t".join(map(str, line_list)) + "\n")

    input_file.close()

    # -------------------------------------------------------->>>>>>
    # close and clear temp file
    # -------------------------------------------------------->>>>>>
    logging.info("Removing temp files ...")

    # time.sleep(20)

    for temp_filename in region_split_filename_list:
        try:
            os.remove(temp_filename)
        except Exception:
            logging.error(f"Removing error on \n\t{temp_filename}")
            raise RemoveTempError

    logging.info("Everything Done!")



if __name__ == '__main__':
    # main(1, 2)
    pass
