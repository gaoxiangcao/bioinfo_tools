import gtftools  # sys.path
from gtftools import Gtf, Gff3


def test_info():
    # 测试基本信息
    print(gtftools.__version__)
    print(gtftools.__author__)
    print(gtftools.__email__)


def test_gff3():
    gff3 = Gff3(file = "../data/test.gff3.gz")
    gff3.to_gtf(file = "./test_to_gtf.gtf")


def test_gtf():
    pass


def main():
    # test_info()
    # test_to_gft()
    pass


if __name__ == "__main__":
    main()
