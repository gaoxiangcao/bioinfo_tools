import os
import logging
import gzip
import string
import random


def split_file_and_make_temp(input_filename,
                             n_part=1,
                             temp_dir=None
                             ):
    """Doc.
    """
    # -------------------------------------------------------->>>>>>
    # set temp_dir
    # -------------------------------------------------------->>>>>>

    if temp_dir:
        # 已经指定
        temp_dir = os.path.abspath(temp_dir)
    else:
        # 未指定， 设置为 input path
        temp_dir = os.path.abspath(os.path.dirname(input_filename))

    if not os.path.exists(temp_dir):
        # 路径不存在
        os.makedirs(os.path.abspath(temp_dir))
        logging.info(f"temp_dir: {temp_dir} not exists, please create ! ")

    input_file_basename = os.path.basename(input_filename)


    logging.info("Make temp files...")
    temp_file_list = []
    temp_filename_list = []

    for index in range(n_part):
        # make filename
        rnd_str = "".join(
            random.sample(
                (string.ascii_letters + string.digits),
                16))
        temp_file_basename = f"temp_{input_file_basename}.{index}.{rnd_str}"
        # print(temp_file_basename)
        temp_file_name = os.path.join(temp_dir, temp_file_basename)
        # print(temp_file_name)
        temp_filename_list.append(temp_file_name)

        # open file
        temp_file_list.append(open(temp_file_name, "wt"))

    # -------------------------------------------------------->>>>>>
    # counting input file line number
    # -------------------------------------------------------->>>>>>
    logging.info("Counting input file...")

    total_input_line_num = 0

    # input_file = open(
    #     input_filename,
    #     "rt") if not input_filename.endwith(".gz") else gzip.open(
    #     input_filename,
    #     "rt")
    input_file = open(input_filename, 'rt') \
        if not input_filename.endswith('.gz') else gzip.open(input_filename, 'rt')

    for line in input_file:
        total_input_line_num += 1

    input_file.close()

    # each_file_line_num = (total_input_line_num // n_part)
    if total_input_line_num % n_part == 0:
        each_file_line_num = (total_input_line_num // n_part)
    else:
        each_file_line_num = (total_input_line_num // n_part) + 1

    logging.info("Done!")
    logging.info("input file total line count: %s" % total_input_line_num)

    # -------------------------------------------------------->>>>>>
    # split temp file
    # -------------------------------------------------------->>>>>>
    # write into output filename
    input_file = open(
        input_filename,
        "rt") if not input_filename.endswith(".gz") else gzip.open(
        input_filename,
        "rt")

    for index, line in enumerate(input_file):
        # print(line.rstrip())
        # print(line)
        file_index = index // each_file_line_num
        temp_file_list[file_index].write(line)

    input_file.close()

    [temp_file.close() for temp_file in temp_file_list]

    logging.info("Make temp files done!")

    return temp_filename_list

