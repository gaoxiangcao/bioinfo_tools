"""
FASTQ 读取
FASTA 读取
...
"""
import gzip


def load_fastx(file: str) -> list:
    """

    :param str file: path of input <fastq/fastaq.gz/fasta/fasta.gz>
    :return:a list for all reads
        i.e.[['header','seq','info','quality'],[...]]
    :rtype: list
    """

    f = open(file, 'rt') if '.gz' not in file else gzip.open(file, 'rt')
    # fastq @ fasta >
    symbol = f.read(1)
    # print(symbol)
    # print(f.read())
    f.close()
    # 当文件很小的时候，直接载入内存就好
    f = open(file, 'rt') if '.gz' not in file else gzip.open(file, 'rt')
    if symbol == '@':
        # FASTQ
        # print('is fastq!')
        # print(f.readlines()) # .readlines方法只能使用一次
        raw_info = [i.rstrip() for i in f.readlines()]
        ls = []
        read = []
        line_fix = None

        for line in raw_info:
            # header line!
            if line.startswith('@'):

                n = len(read)
                # line == []
                if n == 0:
                    read.append(line)
                # line == ['header','seq','info','quality']
                elif n == 4:
                    ls.append(read)
                    read = []
                    read.append(line)
                elif n == 3:
                    # TODO
                    """
                    @@IIEIBCE>IC<IBIIIIEAIEIEB<IDECCD6ICBCED< # 期望是header，但它是quality
                    [
                    '@Beta12AdemL1C001R00100001768/1',
                    'ATCCCCGTATCTTCACCCCACCACAAACTA',
                    '+']
                    @@IIEIBCE>IC<IBIIIIEAIEIEB<IDECCD6ICBCED<
                    """
                    # read.append(line)
                    if line_fix:
                        read.append(line_fix)
                        ls.append(read)
                        read = []
                        read.append(line)
                        line_fix = None
                    else:
                        line_fix = line
                else:
                    f.close()
                    raise ValueError('The file may be not incomplete')
            # not header line!
            else:
                read.append(line)
        f.close()
        # read.append(line)
        ls.append(read)
        return ls
    elif symbol == '>':
        # print('is fasta!')
        ls = []
        read = []
        seq = ''
        raw_info = [i.rstrip() for i in f.readlines()]
        # print(raw_info)
        for line in raw_info:
            if line.startswith('>'):
                # 读取header
                n = len(read)
                if n == 0:
                    read.append(line)
                elif n == 1:
                    # 已经有一个header了，需要添加seq。
                    read.append(seq)  # add seq line
                    ls.append(read)
                    read = []  # 重置 read 这个 list
                    read.append(line)
                    seq = ''  # 重置 seq 这个 str
                else:
                    f.close()
                    raise ValueError('The file may be not incomplete')
            else:
                # 读取并添加seq
                seq += line
        read.append(seq)
        ls.append(read)
        return ls
        # FASTA
    else:
        raise ValueError(
            'Input line one must start with "@" for FQ or ">" for FA')
    f.close()


# 当文件很大的时候，使用生成器函数，产生一个迭代器来进行文件读取
def load_fastx_generator(file):
    """

    :param str file: path of input <fastq/fastaq.gz/fasta/fasta.gz>
    :return:a generator for all reads
        i.e.print(next(obj)) -> ['header','seq','info','quality'],[...]
    :rtype: generator
    """

    f = open(file, 'rt') if not'.gz' in file else gzip.open(file, 'rt')
    # fastq @ fasta >
    symbol = f.read(1)
    # print(symbol)
    # print(f.read())
    f.close()
    # 当文件很小的时候，直接载入内存就好
    f = open(file, 'rt') if '.gz' not in file else gzip.open(file, 'rt')
    if symbol == '@':
        # FASTQ
        # print('is fastq!')
        read = []
        line = f.readline().rstrip()
        while True:
            if not line:
                break
            else:
                # header line!
                if line.startswith('@'):

                    n = len(read)
                    # line == []
                    if n == 0:
                        read.append(line)
                        line = f.readline().rstrip()
                    # line == ['header','seq','info','quality']
                    elif n == 3:
                        # TODO
                        """
                        @@IIEIBCE>IC<IBIIIIEAIEIEB<IDECCD6ICBCED< # 期望是header，但它是quality
                        [
                        '@Beta12AdemL1C001R00100001768/1',
                        'ATCCCCGTATCTTCACCCCACCACAAACTA',
                        '+']
                        @@IIEIBCE>IC<IBIIIIEAIEIEB<IDECCD6ICBCED<
                        """
                        read.append(line)
                        line = f.readline().rstrip()

                        if not line.startswith('@'):
                            raise ValueError('The file may be not incomplete')

                    elif n == 4:
                        yield read
                        read = []
                        read.append(line)
                        line = f.readline().rstrip()
                    else:
                        f.close()
                        raise ValueError('The file may be not incomplete')
                # not header line!
                else:
                    read.append(line)
                    line = f.readline().rstrip()

        yield read
        f.close()
    elif symbol == '>':
        # FASTA
        # print('is fasta!')
        ls = []
        read = []
        seq = ''
        line = f.readline().rstrip()
        while True:
            if not line:
                break
            else:
                if line.startswith('>'):
                    # 读取header
                    n = len(read)
                    if n == 0:
                        read.append(line)
                        line = f.readline().rstrip()
                    elif n == 1:
                        # 已经有一个header了，需要添加seq。
                        read.append(seq)  # add seq line
                        yield read
                        read = []  # 重置 read 这个 list
                        read.append(line)
                        line = f.readline().rstrip()
                        seq = ''  # 重置 seq 这个 str
                    else:
                        f.close()
                        raise ValueError('The file may be not incomplete')
                else:
                    # 读取并添加seq
                    seq += line
                    line = f.readline().rstrip()
        read.append(seq)
        yield read
    else:
        f.close()
        raise ValueError(
            'Input line one must start with "@" for FQ or ">" for FA')
    f.close()


if __name__ == '__main__':
    FQ_TEST = '/home/caogaoxiang/python/Prepared_data/FASTQ/fake_fq.fastq.gz'
    FA_TEST = '/home/caogaoxiang/python/Prepared_data/FASTA/fake_fa.fasta'
    # ls = load_fastx(file=FQ_TEST)
    # print(ls)
    # print(ls[0])
    # print(ls[-1][0])  # 访问header
    # ls = load_fastx(file=FA_TEST)
    # print(ls)
    # it = load_fastx_generator(file=FQ_TEST)
    # print(it)
    # print(next(it))
    # print(next(it))
    # print(next(it))
    # print(next(it))
    # print(next(it))

    # while True:
    #     try:
    #        print(next(it))
    #     except StopIteration:
    #         print('迭代完了')
    #         break

    # it = load_fastx_generator(file=FA_TEST)
    # print(it)
    # print(next(it))
    # print(next(it))
    # print(next(it))
    # print(next(it))
    # print(next(it))
