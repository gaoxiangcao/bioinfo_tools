import gzip


class GeneRecord(object):
    def __init__(self, gene_record: dict):
        self.transcript_id = gene_record['name']
        self.chromsome = gene_record['chrom']
        self.strand = gene_record['strand']

        self.transcript_start = int(gene_record['txStart']) + 1
        self.transcript_end = int(gene_record['txEnd'])

        self.cds_start = int(gene_record['cdsStart']) + 1
        self.cds_end = int(gene_record['cdsEnd'])

        self.exon_count = int(gene_record['exonCount'])
        self.symbol = gene_record['name2']

        self.exon_starts = self._fix_str_coord_list(gene_record['exonStarts'])
        # print(self.exon_starts)
        # 41847188, 41848870, 41918412, 42035806,
        self.exon_starts = [i + 1 for i in self.exon_starts]
        self.exon_ends = self._fix_str_coord_list(gene_record['exonEnds'])
        self.exon_frames = self._fix_str_coord_list(gene_record['exonFrames'])

    @staticmethod
    def _fix_str_coord_list(record: str):
        if record:
            record = record.strip()
            if record.endswith(','):
                return list(map(int, record[:-1].split(",")))
            else:
                return list(map(int, record.split(",")))
        else:
            return None

    def get_exon_lengths_of_one_transcript(self):
        transcript_lengths = []
        exon_5 = 0
        exon_3 = 0
        exon_lengths = []

        if self.exon_count:
            if self.transcript_start < self.exon_starts[0]:
                exon_5 = self.exon_starts[0] - self.transcript_start

            if self.transcript_end > self.exon_ends[-1]:
                exon_3 = self.transcript_end - self.exon_ends[-1]

            for idx in range(self.exon_count):
                exon_length = self.exon_ends[idx] - self.exon_starts[idx] + 1
                # 2 2 1
                exon_lengths.append(exon_length)
            # length of: exon1 exon2 exon3 ...
            # print(f'exon_lengths = {exon_lengths}')

            transcript_lengths.append(exon_5)
            transcript_lengths = transcript_lengths + exon_lengths
            transcript_lengths.append(exon_3)
            # length of: exon5(contains 5utr) + exon1 exon2 ... exon3(contains 3utr)
            # print(f'transcript_lengths = {transcript_lengths}')

            if self.strand == "-":
                # print('!')
                transcript_lengths = transcript_lengths[::-1]
                # print(f'rev transcript_lengths = {transcript_lengths}')
            # print(transcript_lengths)
            return transcript_lengths
        else:
            return None


def parse_database(database: str) -> dict:
    """Parse NCBI gene table.
    Parse tsv table like
    ucsc_hg38_genes-and-gene-predictions_NCBI-ref_seq_refGene_knownGene.tsv

    params
    :param str database: path of database file <tsv|tsv.gz>
    :return: database information
    :rtype: dict
    """
    f = open(database, "rt") if not database.endswith(
        '.gz') else gzip.open(database, 'rt')

    headers = f.readline().rstrip().replace('#', '').split('\t')
    # print(headers)
    # get name_idx
    name_idx = None

    for idx, name in enumerate(headers):
        if name == 'name':
            name_idx = idx

    assert name_idx is not None

    dt_annotation = {}
    # print(dt_annotation)

    for line in f:
        line = line.rstrip().split("\t")
        # print(line)
        name = line[name_idx]
        dt_annotation[name] = {k: v for k, v in zip(headers, line)}
        # print(dt_annotation)
        # break
    f.close()
    return dt_annotation


def parse_aim_transcript(filepath):
    """Parse aim transcript information.
    Parse tsv table like
    NM_1234 123
    NM_1235 13

    params
    :param str filepath: path of aim list <tsv|tsv.gz>
    :return: aim transcript information
    :rtype: dict
    """
    dt_aim_list = {}
    f = open(filepath, "rt") if not filepath.endswith(
        '.gz') else gzip.open(filepath, 'rt')

    for line in f:
        if line.startswith('#'):
            continue
        else:
            line = line.rstrip().split("\t")
            nm_id = line[0]
            rel_coord = line[1]
            key = f'{nm_id}_{rel_coord}'
            dt_aim_list[key] = int(rel_coord)
    f.close()
    # print(dt_aim_list)
    return dt_aim_list


def query_transcripts(
        transcript_id: str,
        transcript_rel_coord: int,
        database: dict):
    """Get exact positions of given transcripts

    params
    :param str transcript_id: like 'NM_123'
    :param int transcript_rel_coord: relative coordinate, must be 1-based
    :param dict database: database dict from function [parse_database]
    :return:
    """
    # fake from transcript
    # transcript_id = 'NM_000027'
    # transcript_rel_coord = 951

    if transcript_id in database:
        pass
    else:
        raise KeyError(f'{transcript_id} not in this database!')

    gene_annotation = database[transcript_id]
    # print(f'gene_annotation = {gene_annotation}')
    gr = GeneRecord(gene_annotation)
    # print(transcript_id)
    gene_relative_cord = gr.get_exon_lengths_of_one_transcript()

    # print(f'gene_relative_cord = {gene_relative_cord}')
    # accumulate_sum
    gene_relative_cord_accumulate = accumulate_sum(gene_relative_cord)
    # print(f'gene_relative_cord_accumulate = {gene_relative_cord_accumulate}')

    # short name
    rel = transcript_rel_coord
    accu_rels = gene_relative_cord_accumulate
    # print(f'accu_rels = {accu_rels}')

    for idx, accumulate_value in enumerate(accu_rels):
        chrom = gr.chromsome
        exact_pos = None

        if accumulate_value < rel:
            # 要查询的相对位置比当前累加长度要长,去下一个循环
            continue
        else:
            """
            要查询的相对位置比当前累加长度要短,开始计算
                exact_pos:  ↓22222        length=34            ↓22256
                ref:        AAAAAAAAGGGGGGGGGCCCCC[C]CCTTTTTTTTT
                0           123456789              ↑22244 (aim ~!)
                1                    0123456789
                2                              012 3 456789
                3                                          01234
                rel:  23                           ↑
            """
            if gr.strand == "+":
                # accu_rels idx 0    1      2      3     4           n-1    n      n + 1
                # length of: exon5 exon1 + exon2 exon3 exon4 ...exon(n-1) exon n  exon3
                # 用上一个 exon 的起始位点开始算 gr.exon_starts[idx - 1]
                # 22222 + 23 - 1 = 22244
                """
                NM_001293562	chr1	+
                transcript start: 33546713  # 	33546710  rel 4 rel2
                transcript end:   33586132
                cds start:        33547850
                cds end:          33585783
                exon start:       33546713,33546991,33547778, ... 33563667,33583502,33585644,
                exon end:         33546895,33547109,33547955, ... 33563780,33583717,33586132,
                """
                if idx == 0:
                    # [3, 185, ...]
                    # rel 2 idx=0 accumulate_value=3 exact_pos=33546711
                    exact_pos = gr.transcript_start + rel - 1
                else:
                    # [3, 185, ...]
                    # rel 4 idx=0 accumulate_value=185 exact_pos=33546714
                    interval = rel - accu_rels[idx - 1]
                    exact_pos = gr.exon_starts[idx - 1] + interval - 1
            elif gr.strand == "-":
                """
                要查询的相对位置比当前累加长度要短,开始计算
                    exact_pos:  ↓22222        length=34            ↓22256
                    ref:        AAAAAAAAGGGGGGGGGCCCCC[C]CCTTTTTTTTT
                    0           123456789              ↑22244 (aim ~!)
                    1                    0123456789
                    2                              012 3 456789
                    3                                          01234
                    rel:  23                           ↑
                    rel_rev: 12
                """
                # accu_rels idx 0    1      2        3           n-1    n       n + 1
                # length of: exon3 exon(n),exon(n-1),exon(n-2)...exon2  exon 1  exon5
                # exon        [exon]                          []         []
                """
                NM_001323573 chr1 -
                transcript start: 48998526
                transcript end:   50489626  # 50489629
                cds start:        48999844
                cds end:          50489468
                exon start:       48998526, 49005313, 49052675, ... 50162948, 50317067, 50489434,
                exon end:         48999965, 49005410, 49052838, ... 50163109, 50317190, 50489626,
                """
                if idx == 0:
                    # [3, 442, ...]
                    # rel 2 idx=0 accumulate_value=3 exact_pos=50489628
                    exact_pos = gr.transcript_end - rel + 1
                else:
                    # [3, 442, ...]
                    # rel 4 idx=0 accumulate_value=185 exact_pos=50489626
                    interval = rel - accu_rels[idx - 1] + 1
                    exact_pos = gr.exon_ends[-idx] - interval
            break
    return chrom, exact_pos, gr.strand, idx


def accumulate_sum(numbers):
    accumulates = [0] * len(numbers)
    for idx, value in enumerate(numbers):
        if idx == 0:
            accumulates[idx] = value
        else:
            accumulates[idx] = accumulates[idx - 1] + value
    return accumulates


def main(aim_transcripts, database, out_path='successful_query_table.csv'):
    f = open(out_path, 'wt')
    f.write('id,transcript_position,chromsome,exact_position\n')
    # parse database
    database = parse_database(database)
    # print(str(database)[:1000])
    # parse aim list
    transcripts = parse_aim_transcript(aim_transcripts)
    # print(transcripts)

    transcripts_failed = []

    # near_seq = 10
    # index = 1
    # 实现核心功能的函数
    for nm_id, relative_cord in transcripts.items():
        nm_id = '_'.join(nm_id.split('_')[:-1])
        # print(f'nm_id = {nm_id}')
        # print(f'relative_cord = {relative_cord}')
        try:
            query_result = query_transcripts(
                transcript_id=nm_id,
                transcript_rel_coord=relative_cord,
                database=database)
            print(f'{nm_id}, {query_result}')
            print(f'for_IGV_check: {query_result[0]}\t{query_result[1] - 20}\t{query_result[1] + 20}')
            # f.write(
            #     f'{nm_id},{relative_cord},{query_result[0]},{query_result[1]}\n'
            # )
            # break

        except KeyError:
            transcripts_failed.append(nm_id)
            continue
    print()
    print('Failed ids:')
    for i in transcripts_failed:
        print(f'\t{i}')

    f.close()


if __name__ == '__main__':
    # 给相对的位置
    # 返回绝对坐标
    # 201283508
    # 201283531
    # 8817626
    main(
        aim_transcripts='pus7_dependent_pseudo_u.table',
        database='ucsc/ucsc_hg38_genes-and-gene-predictions_NCBI-refseq_refGene_knownGene.tsv.gz'
    )

