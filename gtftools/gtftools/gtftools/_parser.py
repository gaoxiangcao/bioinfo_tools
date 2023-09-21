"""Doc.
GFF3, GTF
col 9s
col 9
col 1-8

- class Father # 解析文件
- class Son  # 解析文件细节
"""
import gzip
import gtftools


class Annotation:
    """Father class of GFF3 or GTF.
    Loader of files.
    """
    REQUIRED_COLUMNS = [
        "seq_name", "source", "feature",
        "start", "end", "score",
        "strand", "frame", "attribute"
    ]

    def __init__(self, file):
        # 装载文件
        self.__file = file
        self._lines = None
        self.__load_file()

    def __load_file(self):
        self._lines = open(
            self.__file,
            "rt") if ".gz" not in self.__file else gzip.open(
            self.__file,
            "rt")


class Gff3(Annotation):
    """Gff3 class.
    Parse Gff3 files.
    """
    DROP_ATTRIBUTE = (
        # chromosome
        'chromosome', 'scaffold',
        # gene
        'gene', 'pseudogene', 'ncRNA_gene', 'V_gene_segment',
        'D_gene_segment', 'J_gene_segment', 'C_gene_segment',
        #
        'exon',
        #
        'biological_region', 'pseudogenic_transcript', 'unconfirmed_transcript',
    )

    def to_gtf(self, file):
        print("Convert Gff3 obj to Gtf obj...")
        f = open(file, "wt") if ".gz" not in file else gzip.open(file, "wt")

        counter_line = 0

        for line in self._lines:
            line = line.strip()  # 1 based line
            counter_line += 1

            if line:
                if line.startswith("#"):
                    print(f"Skip # annotation{counter_line}")  # logging
                    continue
                else:
                    try:
                        seq_name, source, feature, start, end, score, strand, frame, attribute = line.split(
                            "\t")
                    except Exception as e:
                        raise gtftools.ParsingError(str(e))

                        attribute2 = {}

                        for i in attribute.split(";"):
                            _ls = i.split("=")
                            if len(_ls) != 2:
                                continue
                            elif len(_ls) == 2:
                                k, v = _ls
                                attribute2[k] = v

                        attribute = attribute2

                        _dt = {
                            k: v for k,
                            v in zip(
                                self.REQUIRED_COLUMNS,
                                (seq_name,
                                 source,
                                 feature,
                                 start,
                                 end,
                                 score,
                                 strand,
                                 frame,
                                 attribute))}
                        # print(_dict)
                        # 逻辑判断，写个方法，处理_dict, 保持代码整洁
                        final_dict = self.__to_gtf(dt=_dt)
                        # 处理成功
                        if final_dict:
                            # print(final_dict)
                            line = ""

                            for k, v in final_dict.items():
                                if k != self.REQUIRED_COLUMNS[8]:
                                    line += f"{v}\t"
                                else:
                                    attri = ""

                                    for k2, v2 in v.items():
                                        attri += f'{k2} "{v2}"; '
                                    attri = attri.rstrip()
                                    line += attri + "\n"

                            f.write(line)

                        # final_dict  = None
                        else:
                            continue

            else:
                continue

        f.close()

    def __to_gtf(self, dt: dict) -> None:

        feature = dt[self.REQUIRED_COLUMNS[2]]

        if feature in self.DROP_ATTRIBUTE:
            return
        else:
            attributes = dt[self.REQUIRED_COLUMNS[8]]

            # five_prime_UTR Parent=transcript:ENST00000673477
            # three_prime_UTR Parent=transcript:ENST00000673477
            # CDS
            # ID=CDS:ENSP00000500094;Parent=transcript:ENST00000673477;protein_id=ENSP00000500094

            # mRNA
            # ID=transcript:ENST00000673477;Parent=gene:ENSG00000160072;Name=ATAD3B-206;biotype=protein_coding;

            # UTR ID没有， transcript_id 在 parent里，gene_id也没有
            # CDS ID有， transcript_id 在 parent里，gene_id也没有
            # RNA snRNA, miRNA... ID有， transcript_id 在 parent里，gene_id在parent里
            # print("-" * 10)
            # print("DEBUG", feature, attributes)

            _id = attributes["ID"] if "ID" in attributes else None
            parent = attributes["Parent"] if "Parent" in attributes else None

            gene_id = ""
            transcript_id = ""

            if feature in ["five_prime_UTR", "three_prime_UTR", "CDS"]:
                _dt = {
                    "five_prime_UTR": "five_prime_utr",
                    "three_prime_UTR": "three_prime_utr",
                    "CDS": "CDS"}
                feature = _dt[feature]
                transcript_id = parent.split(
                    ":")[1] if "transcript" in parent else ""
                gene_id = ""

            elif "RNA" in feature:
                transcript_id = _id.split(
                    ":")[1] if "transcript" in _id else ""
                gene_id = parent.split(
                    ":")[1] if parent.startswith("gene") else ""

            # print(f"{feature}, {transcript_id}, {gene_id}")
            dt_info = {
                "gene_id": gene_id,
                "transcript_id": transcript_id
            }
            for k, v in attributes.items():
                if k[0].islower():
                    dt_info[k] = v

            dt[self.REQUIRED_COLUMNS[2]] = feature
            dt[self.REQUIRED_COLUMNS[8]] = dt_info
            return dt


class Gtf(Annotation):
    """Gtf class.
    Parse Gtf files.
    """
    pass
