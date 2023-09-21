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
        print(self.exon_starts)
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
            # print(f'for_IGV_check: {query_result[0]}\t{query_result[1] - 20}\t{query_result[1] + 20}')
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

