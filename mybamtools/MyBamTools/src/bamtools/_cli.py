from bamtools.count_use_bed import main as run_count_use_bed
from fire import Fire
# import sys
# import os
# sys.path.append("/home/caogaoxiang/python/project_5/MyBamTools/src/")
# sys.path.append("/home/caogaoxiang/python/project_5/MyBamTools/src/bamtools")
# sys.path.append("/home/caogaoxiang/python/project_5/MyBamTools/src/bamtools")
# sys.path.append("/home/caogaoxiang/python/project_5/MyBamTools")
# sys.path.append("/home/caogaoxiang/python/project_5/MyBamTools/src")
# sys.path.append(os.path.dirname(sys.path[0]))
# print(sys.path)
# __dir__ = os.path.dirname(os.path.abspath("/home/caogaoxiang/python/project_5/MyBamTools/src/bamtools/_cli"))
# sys.path.append(__dir__)
# sys.path.append(os.path.abspath(os.path.join(__dir__, '../..')))


class Cli(object):
    """Cli interface of python package <bamtools>.

    - bamtools is a commandline tool and a python package.
    """

    def count_use_bed(
            self,
            input: str,
            bams: str,
            tags: str,
            output: str = None,
            method: str = "CPM",
            process: int = 1,
            mapq: int = 20,
            strand: str = "b",
            extend_length: int = 1,
            temp_dir: str = None,
            verbose: str = "ERROR",
    ):
        """Doc of count_use_bed

        :param str input: path of input <bed|bed.gz>
        :param str output: path of output tsv table
        :param list bams: BAM files. bam paths joint by comma, e.g. <1.bam,2.bam,3.bam>
        :param list tags: BAM tags. tags joint by comma, e.g. <H3K4me3,H4K36me2,H3K9me3>
        :param str method: RPKM, Raw or CPM
        :param int process: Process to use, must ≤ number of CPUs
        :param int mapq: MAPQ cutoff
        :param str strand: strand selection method, can be f, r, b, or s
        :param int extend_length: region extend length
        :param str temp_dir: folder to put temp files, use input dir if not defined
        :param str verbose: 'CRITICAL', 'FATAL', 'ERROR', 'WARNING', 'WARN', 'INFO', 'DEBUG', 'NOTSET'
        """
        return run_count_use_bed(
            input=input,
            output=output,
            bams=bams,
            tags=tags,
            method=method,
            process=process,
            mapq=mapq,
            strand=strand,
            extend_length=extend_length,
            temp_dir=temp_dir,
            verbose=verbose)

    def sub_cmd2(self):
        """Description of sub_cmd2
        """
        return "python"


def main():
    cli = Cli()
    Fire(cli, name="bamtools")


if __name__ == "__main__":
    main()
