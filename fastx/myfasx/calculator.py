"""
1.读取FASTQ文件并且以FASTA格式输出
2.解析FASTQ的质量值，并计算Q30的比例
"""
import gzip
import random
from coden import CODEN
from tqdm import tqdm
from loader import load_fastx,load_fastx_generator

# 面向对象 ！ ~ fasta Genome Transcriptome
class Genome(object):
    """Genome class
    """
    PAIR = {k:v for k, v in zip('ATGCUN','TACGAN')}


    def __init__(self, file):
        self.__FILE_PATH = file
        # load genome
        self.GENOME = {}
        self.__parse_genome()
        self.__genome_length = {}
        self.__gc_ratio = {}
        self.__effective_length = {}
        # print(self.GENOME.keys())
        # print(self.GENOME)


    def __parse_genome(self):
        print('parse genome...')

        for chrom, seq in tqdm(load_fastx(file=self.__FILE_PATH)):
            chrom = chrom[1:] # >chr
            self.GENOME[chrom] = seq.upper()
                # .replace('U', 'T')


    def replace_base(self, *args, convert:dict):
        # {'A':'T', 'U':'T'}
        # U to T
        # C to T
        print('replace base...')
        if not args:
            args=self.GENOME.keys()
            print(f'chromosomes are not specified, use all: {list(args)}')
        else:
            print(f'chromosomes are  specified, use : {sorted(list(args))}')

        for chrom in args:
            for k,v in convert.items():
                self.GENOME[chrom] = self.GENOME[chrom].replace(k,v)


    def reverse_complement(self):
        print('Reverse complement')
        # reverse
        for chrom in tqdm(self.GENOME.keys()):
            self.GENOME[chrom] = self.GENOME[chrom][::-1]
        # complemented
            self.GENOME[chrom] = ''.join([self.PAIR[base] for base in self.GENOME[chrom]])
            print(self.GENOME[chrom])


    def calculate_length(self):
        print('calculate_length......')
        for chrom in self.GENOME.keys():
            self.__genome_length[chrom] = len(self.GENOME[chrom])


    def get_length(self):
        self.calculate_length()
        return self.__genome_length


    def __calculate_gc_ratio(self):
        print('calculate_gc_ratio......')
        for chrom in self.GENOME.keys():
            seq = self.GENOME[chrom]
            g = seq.count('G')
            c = seq.count('C')
            total = len(seq)
            self.__gc_ratio[chrom] = (g + c) / total


    def get_gc_ratio(self):
        self.__calculate_gc_ratio()
        return self.__gc_ratio


    def __calculate_effective_length(self):
        print('calculate_effective_length......')
        for chrom in self.GENOME.keys():
            seq = self.GENOME[chrom]
            n = seq.count('N')
            total = len(seq)
            self.__effective_length[chrom] =  total - n


    def get_effective_length(self):
        self.__calculate_effective_length()
        return self.__effective_length

# 1.读取FASTQ文件并且以FASTA格式输出
def fastq_to_fasta(file, out_name, use_iter=True):
    """
    # Convert FASTQ to FASTA
    :param str file: path of a <fastq | fastq.gz>
    :param str out_name: path of a <fasta | fasta.gz>
    :param bool use_iter: True -> iterator, False -> list to loop
    :return: None
    """
    if use_iter:
        #  use iterator
        reads = load_fastx_generator(file=file)
    else:
        #  use list
        reads = load_fastx(file=file)
    # print(reads)
    # print(type(reads))

    f = open(out_name,'wt') if not '.gz' in out_name else gzip.open(out_name,'wt')
    # f = open(out_name, 'wt') if '.gz' not in out_name else gzip.open(out_name, 'wt')
    # for read in reads:
    #     print(read)
    for header, seq, _, _, in reads:
        # print(header, seq)
        header = '>' + header[1:]
        # print(header,seq)
        f.write(
            f'{header}\n{seq}\n'
        )

    f.close()

    print("Covert done!")


# 2.解析FASTQ的质量值，并计算Q30的比例, Qwhat?, >Q?, <Q?
def get_aim_quality_ratio(file, quality=30, method='>Q', use_iter=True):
    """
    # Calculator aim quality ratio > quality or < quality
    :param str file: path of a <fastq | fastq.gz>
    :param int quality: the quality to compare with (should > 0)
    :param str method: how to compare
    :param bool use_iter: True -> iterator, False -> list to loop
    :return: aim quality ratio
    :rtype: float
    """
    if use_iter:
        #  use iterator
        reads = load_fastx_generator(file=file)
    else:
        #  use list.
        reads = load_fastx(file=file)

    # Phred (quality) = -10 * log10(errorP)  # 以illumina为例，0.001， -3 * -10 = 30
    # 30 + 33 -> ASCII
    # print(next(reads))
    # Phred(quality) = -10 * log10(errorP)  # errorP，以 illumina 为例, 0.001, -3 * - 10 = 30
    # 30 + 33 -> ASCII  # 0~ 127 = 128
    # Q/Phred + 33 -> ASCII
    # Q/Phred ASCII->value - 33
    total_base = 0
    total_base_aim = 0

    for _, _, _, r_quality in reads:
        q = [ord(base) - 33 for base in r_quality]

        if method == '>Q':
            q_aim = [i for i in q if i > quality]
        elif method == '<Q':
            q_aim = [i for i in q if i < quality]
        else:
            raise ValueError('Param method is wrong!')

        # print(q)
        # print(q_aim)
        total_base += len(q)
        total_base_aim += len(q_aim)

    return total_base_aim / total_base


# 3. trim fastq
def trim_fastq(file, out_name, trim_start:int=0, trim_end:int=None, use_iter:bool=True):
    """
    # trim_fastq
    :param file:
    :param out_name:
    :param trim_start:
    :param trim_end:
    :param use_iter:
    :return:
    """
    if use_iter:
        #  use iterator
        reads = load_fastx_generator(file=file)
    else:
        #  use list
        reads = load_fastx(file=file)

    try:
        assert 0 <= trim_start < trim_end
    except AssertionError:
        print('Must follow this: 0 <= trim_start < trim_end')


    f = open(out_name, 'wt') if not '.gz' in out_name else gzip.open(out_name, 'wt')
    # f = open(out_name, 'wt') if '.gz' not in out_name else gzip.open(out_name, 'wt')
    # for read in reads:
    #     print(read)
    for header, seq, info, r_quality, in reads:
        # print(header, seq, info, r_quality,)
        # break
        # AGCTACTAAACCCCC
        # 012345678910
        seq = seq[trim_start:trim_end] # [) step1
        r_quality = r_quality[trim_start:trim_end] # [) step1

        f.write(
            f'{header}\n{seq}\n{info}\n{r_quality}\n'
        )

    f.close()
    print("Trim done!")


# 4. filter fastq
def filter_fastq(file, out_name, tiles_to_drop:list, use_iter:bool=True):
    """

    :param file:
    :param out_name:
    :param tiles_to_drop:
    :param use_iter:
    :return:
    """


    if use_iter:
        #  use iterator
        reads = load_fastx_generator(file=file)
    else:
        #  use list
        reads = load_fastx(file=file)


    f = open(out_name, 'wt') if not '.gz' in out_name else gzip.open(out_name, 'wt')

    tile = None
    counter_dropped_reads = 0
    counter_all_reads = 0

    # header: illumina style! mgi no! fake no!
    for header, seq, info, r_quality, in reads:
        # print(header)
        counter_all_reads += 1

        try:
            tile = int(header.split('\t')[0].split(':')[-3])
        except IndexError:
            print('Parse <header> failed\n'
                  'Please make sure this is an illumina NGS file!\n'
                  f'\t<header>:{header}\n')
        # print(tile)
        if tile not in tiles_to_drop:
            # 写入
            f.write(
                f'{header}\n{seq}\n{info}\n{r_quality}\n'
            )
        else:
            # 写出去
            counter_dropped_reads += 1

    f.close()

    print('filter done')
    print(f'{counter_dropped_reads}/{counter_all_reads} reads\n')
    print(f'{counter_dropped_reads/counter_all_reads:.3%} reads were dropped!')


# 5. translate
# fasta -> AA
def translate(file, out_name, use_iter:bool=True):
    """

    :param file:
    :param out_name:
    :param use_iter:
    :return:
    """


    if use_iter:
        #  use iterator
        reads = load_fastx_generator(file=file)
    else:
        #  use list
        reads = load_fastx(file=file)

    dt_start_coden = {k:v for k,v in CODEN.items() if '#' in v}
    dt_stop_coden = {k:v for k, v in CODEN.items() if '$' in v}
    print(dt_start_coden)
    print(dt_stop_coden)
    print(CODEN)
    # seq.upper()
    # T -> U
    # start coden
    # stop coden

    f = open(out_name, 'wt') if not '.gz' in out_name else gzip.open(out_name, 'wt')

    for header, seq in reads:
        # 只考虑AUG
        # fix AUG
        seq = seq.upper().replace('T','U')
        start = seq.find('AUG')
        seq = seq[start:]
        seq_aa = ''
        while True:
            if len(seq) >= 3:
                aa = CODEN[seq[:3]].replace('#', '').replace('$', '')
                seq = seq[3:]

                if aa:
                    # 正常密码子
                    seq_aa += aa
                else:
                    # stop coden!
                    print('Stop!')
                    seq_aa += '\n'
                    break
            else:
                seq_aa += '\n'
                break
        f.write(f'{header} translate to AA\n{seq_aa}\n')
    f.close()


def down_sampling(file, out_name, ratio: float = None, number: int = None):
    """

    :param file:
    :param out_name:
    :param ratio:
    :param number:
    :return:
    """

    # reads = load_fastx(file=file)
    # count_all_reads = len(reads)
    reads = load_fastx_generator(file=file)
    count_all_reads = 0

    for read in reads:
        count_all_reads += 1

    reads = load_fastx_generator(file=file)


    if ratio and number:
        raise ValueError('Only one of ratio and number can be defined')
    elif ratio:
        assert 0 <= ratio <=1
        number = int(count_all_reads * ratio)
    elif number:
        if number > count_all_reads:
            raise ValueError(f'number must <= total: {count_all_reads}')
    else:
        raise ValueError('Only/Must one of ratio and number can be defined')

    # number
    # random select
    to_be_select = random.sample(range(count_all_reads), number)
    # print(to_be_select)
    to_be_select.sort(reverse=True)
    # print(to_be_select)

    f = open(out_name, 'wt') if not '.gz' in out_name else gzip.open(out_name, 'wt')

    idx = 0
    select = to_be_select.pop()

    try:
        for header, seq, info, r_quality, in reads:
            # print(idx, select,to_be_select)
            if idx == select:
                # write
                f.write(f'{header}\n{seq}\n{info}\n{r_quality}\n')
                idx += 1
                select = to_be_select.pop()
            else:
                idx += 1
                continue

    except IndexError:
        print('Down sampling done!')

    f.close()
    print(f'\t{number} / {count_all_reads} reads {number/count_all_reads: .3%} were selected')


if __name__ == '__main__':
    # ------------------------------------------------------------------->>>>>>>>>>
    # files
    # ------------------------------------------------------------------->>>>>>>>>>
    FQ_TEST = '/home/caogaoxiang/python/Prepared_data/FASTQ/fake_fq.fastq'
    FA_TEST = '/home/caogaoxiang/python/Prepared_data/FASTA/fake_fa.fasta'
    FQ_TEST2 = '/home/caogaoxiang/python/Prepared_data/FASTQ/from_illumina/from_illumina_R1.fastq.gz'
    FQ_TEST3 = '/home/caogaoxiang/python/Prepared_data/FASTQ/from_mgi/from_mgi_R1.fastq.gz'
    FA_TEST2 = '/home/caogaoxiang/python/Prepared_data/FASTA/mRNA_CTCF.fasta'
    GENOME_TEST = '/home/caogaoxiang/python/Prepared_data/FASTA/genome_XY_for_test.fa.gz'
    GENOME = '/home/caogaoxiang/python/Prepared_data/FASTA//genome_ucsc_mm39.fa.gz'
    GENOME_FAKE = '/home/caogaoxiang/python/Prepared_data/FASTA/genome_fake.fa'
    # ------------------------------------------------------------------->>>>>>>>>>
    # fastq_to_fasta
    # ------------------------------------------------------------------->>>>>>>>>>
    # fastq_to_fasta(file=FQ_TEST, out_name='test_fastq_to_fasta.fasta.gz')
    # fastq_to_fasta(file=FQ_TEST, out_name='test_fastq_to_fasta.fasta', use_iter=False)
    # fastq_to_fasta(file=FQ_TEST, out_name='test_fastq_to_fasta.fasta', use_iter=True)
    # TODO Fix MGI line4 @@
    # fastq_to_fasta(file=FQ_TEST3, out_name='test_fastq_to_fasta3.fasta', use_iter=True)
    # fastq_to_fasta(file=FQ_TEST3, out_name='test_fastq_to_fasta3.fasta', use_iter=False)
    # ------------------------------------------------------------------->>>>>>>>>>
    # get_aim_quality_ratio
    # ------------------------------------------------------------------->>>>>>>>>>
    # print(get_aim_quality_ratio(file=FQ_TEST, use_iter=False))
    # print(get_aim_quality_ratio(file=FQ_TEST2, quality=30, method='>Q', use_iter=True))
    # print(get_aim_quality_ratio(file=FQ_TEST2, quality=30, method='>Q', use_iter=False))
    # print(get_aim_quality_ratio(file=FQ_TEST2, quality=20, use_iter=True))
    # print(get_aim_quality_ratio(file=FQ_TEST2, quality=20, method='<Q', use_iter=True))
    # print(get_aim_quality_ratio(file=FQ_TEST3, quality=20, method='<Q', use_iter=True))
    # TODO Fix MGI line4 @@
    # print(get_aim_quality_ratio(file=FQ_TEST3, quality=20, method='<Q', use_iter=True))
    # print(get_aim_quality_ratio(file=FQ_TEST3, quality=20, method='<Q', use_iter=True))
    # print(get_aim_quality_ratio(file=FQ_TEST3, quality=20, method='<Q', use_iter=False))
    # ------------------------------------------------------------------->>>>>>>>>>
    # trim_fastq
    # ------------------------------------------------------------------->>>>>>>>>>
    # trim_fastq(file=FQ_TEST2, out_name='test_trim_fastq2.fastq', trim_start=0, trim_end=109)
    # fake
    # trim_fastq(file=FQ_TEST, out_name='test_trim_fastq.fastq', trim_start=1, trim_end=10)
    # mgi
    # trim_fastq(file=FQ_TEST3, out_name='test_trim_fastq3.fastq', trim_start=0, trim_end=109)
    # ------------------------------------------------------------------->>>>>>>>>>
    # filter_fastq
    # ------------------------------------------------------------------->>>>>>>>>>
    # illumina
    # filter_fastq(file=FQ_TEST2, out_name='test_filter_fastq2.fastq.gz', tiles_to_drop=[2201,1116])
    # fix bug (header has no tile info!)
    # filter_fastq(file=FQ_TEST, out_name='test_filter_fastq.fastq.gz', tiles_to_drop=[2201,1116])
    # filter_fastq(file=FQ_TEST3, out_name='test_filter_fastq3.fastq.gz', tiles_to_drop=[2201,1116])
    # ------------------------------------------------------------------->>>>>>>>>>
    # translate
    # ------------------------------------------------------------------->>>>>>>>>>
    # translate(file=FA_TEST2, out_name='test_translate.fa')
    # ------------------------------------------------------------------->>>>>>>>>>
    # Genome Ojb!
    # ------------------------------------------------------------------->>>>>>>>>>
    # genome = Genome(file=GENOME)
    # genome = Genome(file=GENOME_TEST)
    # genome = Genome(file=GENOME)
    # genome = Genome(file=GENOME_FAKE)
    # print(genome.GENOME)
    # genome.replace_base(convert={'C':'T','U':'T'})
    # print(genome.GENOME)

    # genome = Genome(file=GENOME_FAKE)
    # print(genome.GENOME)
    # genome.replace_base('chrM', 'chr2', convert={'C':'T','U':'T'})
    # print(genome.GENOME)
    # 访问隐藏的变量
    # print(genome._Genome__FILE_PATH)
    # genome.reverse_complement()
    # print(genome.get_length())
    # print(genome.get_gc_ratio())
    # print(genome.get_effective_length())
    # down_sampling
    # down_sampling(file=FQ_TEST3, out_name='test_down_sampling.fastq.gz', ratio=0.1, number=20)
    down_sampling(file=FQ_TEST, out_name='test_down_sampling.fastq', ratio=0.6)

