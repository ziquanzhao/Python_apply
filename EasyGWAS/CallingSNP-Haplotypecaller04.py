import argparse
import os
import re


def main(ChromosomeLabel, ChromosomeNumber, ExecutingScript):

    current_directory = os.getcwd()
    if not os.path.exists(f'{current_directory}/Haplotypecaller04'):
        os.mkdir(f'{current_directory}/Haplotypecaller04')

    with open(f'{current_directory}/Haplotypecaller04/Haplotypecaller04.sh', 'w') as hap:
        hap.write('#!/bin/bash\n')
        for filename in os.listdir(f'{current_directory}/SamToBamSortDuplicateIndex03/'):
            if '.bai' not in filename and '.metrics' not in filename and '_SortDup.bam' in filename:
                filename = re.sub('_SortDup.bam', '', filename)
                for chrnum in range(1, ChromosomeNumber + 1, 1):
                    hap.write(
                        f'gatk HaplotypeCaller -R {current_directory}/BuildRefseqIndex02/YourReferenceGenome.fa -I {current_directory}/SamToBamSortDuplicateIndex03/{filename}_SortDup.bam -ERC GVCF -L {ChromosomeLabel}{chrnum} -O {current_directory}/Haplotypecaller04/{filename}_chr{chrnum}.g.vcf\n')

    if ExecutingScript == "no" or ExecutingScript == "NO" or ExecutingScript == "No" or ExecutingScript == "N" or ExecutingScript == "n":
        print("")
        print(f"\033[1;34m在{current_directory}/Haplotypecaller04/目录下有个Haplotypecaller04.sh脚本，你可以先检查一下其内容是否正确，然后手动执行它! 也可以在原有命令上添加'-es y'参数直接执行。\033[0m")
        print("")
        pass
    elif ExecutingScript == "yes" or ExecutingScript == "YES" or ExecutingScript == "Yes" or ExecutingScript == "Y" or ExecutingScript == "y":
        os.system(f"bash {current_directory}/Haplotypecaller04/Haplotypecaller04.sh")
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
        "1.此步骤是calling SNP。\n"
        "2.需要使用-cl参数输入参考基因组中染色体命名标签。例如:基因组中染色体以Chr1、Chr2、Chr3命名，则输入：Chr\n"
        "3.需要使用-cn参数输入参考基因组中染色体条数。\n"
        "4. -ERC GVCF(默认)是指GVCF文件记录非突变位点的时候，以块的形式来记录； -ERC BP_RESOLUTION 是指GVCF文件对非突变和突变位点一视同仁,都以点的形式记录"
        "\033[0m")

    parser = argparse.ArgumentParser(description=description)

    # 接受的输入参数
    parser.add_argument("-cl", metavar="ChromosomeLabel", type=str, required=True,
                        help="输入你的参考基因组中的染色体标签 (例如基因组中染色体以Chr1、Chr2、Chr3命名，则输入：Chr)")
    parser.add_argument("-cn", metavar="ChromosomeNumber", type=int, required=True,
                        help="输入你的参考基因组中的染色体条数 (例如：15)")
    parser.add_argument("-es", metavar="ExecutingScript", type=str, default="no",
                        help="是否直接执行脚本,选择yes或者no，默认是no (默认: %(default)s)")

    # 解析命令行参数
    args = parser.parse_args()

    # 调用主函数进行处理
    main(args.cl, args.cn, args.es)
