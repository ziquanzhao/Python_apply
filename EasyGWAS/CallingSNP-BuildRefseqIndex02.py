import argparse
import os


def main(ReferenceGenome, ExecutingScript):

    current_directory = os.getcwd()
    if not os.path.exists(f'{current_directory}/BuildRefseqIndex02'):
        os.mkdir(f'{current_directory}/BuildRefseqIndex02')

    os.system(f"cp {ReferenceGenome} {current_directory}/BuildRefseqIndex02/YourReferenceGenome.fa")

    with open(f'{current_directory}/BuildRefseqIndex02/BuildRefseqIndex02.sh', 'w') as refseq_index:
        refseq_index.write('#!/bin/bash\n')
        refseq_index.write(f'bwa index -a bwtsw {current_directory}/BuildRefseqIndex02/YourReferenceGenome.fa -p {current_directory}/BuildRefseqIndex02/YourReferenceGenome.fa\n')
        refseq_index.write(f'samtools faidx {current_directory}/BuildRefseqIndex02/YourReferenceGenome.fa --fai-idx YourReferenceGenome.fa.fai\n')
        refseq_index.write(f'mv {current_directory}/YourReferenceGenome.fa.fai {current_directory}/BuildRefseqIndex02/\n')
        refseq_index.write(f'gatk CreateSequenceDictionary -REFERENCE {current_directory}/BuildRefseqIndex02/YourReferenceGenome.fa -OUTPUT {current_directory}/BuildRefseqIndex02/YourReferenceGenome.dict')

    if ExecutingScript == "no" or ExecutingScript == "NO" or ExecutingScript == "No" or ExecutingScript == "N" or ExecutingScript == "n":
        print("")
        print(f"\033[1;34m在{current_directory}/BuildRefseqIndex02/目录下有个BuildRefseqIndex02.sh脚本，你可以先检查一下其内容是否正确，然后手动执行它! 也可以在原有命令上添加'-es y'参数直接执行。\033[0m")
        print("")
        pass
    elif ExecutingScript == "yes" or ExecutingScript == "YES" or ExecutingScript == "Yes" or ExecutingScript == "Y" or ExecutingScript == "y":
        os.system(f"bash {current_directory}/BuildRefseqIndex02/BuildRefseqIndex02.sh")
        print("")
        print("\033[1;34mThis mission has been accomplished!\033[0m")
        print("")
    else:
        print("")
        print("\033[1;34m-es参数是用于判断你是否直接执行脚本的，你只能输入no/NO/No/N/yes/YES/Yes/Y,不能输入其他的内容，请重新设置参数!\033[0m")
        print("")


if __name__ == "__main__":
    # 创建ArgumentParser对象
    description = ("\033[1;31m"
        "1.此步是建立calling SNP过程中所需要的参考基因组的各种索引。\n"
        "2.需要使用-rg参数输入参考基因组文件的绝对路径，该脚本会自动生成calling SNP过程中所需要的基因组索引文件。\n"
        "3.该脚本是调用bwa、samtools、gatk软件执行命令的，你需要安装好这些软件，并设置好环境变量。"
        "\033[0m")

    parser = argparse.ArgumentParser(description=description)

    # 接受的输入参数
    parser.add_argument("-rg", metavar="ReferenceGenome", type=str, required=True,
                        help="输入参考基因组文件的绝对路径 (例如: /mnt/d/cleandata/xso_genome.fa)")
    parser.add_argument("-es", metavar="ExecutingScript", type=str, default="no",
                        help="是否直接执行脚本,选择yes或者no，默认是no (默认: %(default)s)")

    # 解析命令行参数
    args = parser.parse_args()

    # 调用主函数进行处理
    main(args.rg, args.es)
