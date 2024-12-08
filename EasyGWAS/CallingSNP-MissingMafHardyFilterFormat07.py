import argparse
import os
import re


def main(ResultFileName, geno, maf, mind, hwe, min_alleles, max_alleles, thread, format, ExecutingScript):

    current_directory = os.getcwd()
    if not os.path.exists(f'{current_directory}/MissingMafHardyFilterFormat07'):
        os.mkdir(f'{current_directory}/MissingMafHardyFilterFormat07')

    with open(f'{current_directory}/MissingMafHardyFilterFormat07/MissingMafHardyFilterFormat07.sh', 'w') as missMafHardy:
        missMafHardy.write('#!/bin/bash\n')
        missMafHardy.write(
            f'plink1.9 --vcf {current_directory}/MergeHardFilter06/AllSample_HardFilter.vcf.gz --export vcf --allow-extra-chr '
            f'--geno {geno} --mind {mind} --maf {maf} --hwe {hwe} '
            f'--min-ac {min_alleles} --max-ac {max_alleles} --out {current_directory}/MissingMafHardyFilterFormat07/{ResultFileName} --threads {thread}\n')
        if format == "vcf":
            missMafHardy.write(f'bgzip {current_directory}/MissingMafHardyFilterFormat07/{ResultFileName}.vcf\n')
        elif format == "bed":
            missMafHardy.write(f'plink1.9 --vcf {current_directory}/MissingMafHardyFilterFormat07/{ResultFileName}.vcf --vcf-half-call missing  --allow-extra-chr --make-bed --out {current_directory}/MissingMafHardyFilterFormat07/{ResultFileName}\n')
        elif format == "hapmap":
            missMafHardy.write(f'run_pipeline.pl -Xms20g -Xmx100g -vcf {current_directory}/MissingMafHardyFilterFormat07/{ResultFileName}.vcf -sortPositions -export {current_directory}/MissingMafHardyFilterFormat07/{ResultFileName} -exportType HapmapDiploid\n')

    if ExecutingScript == "no" or ExecutingScript == "NO" or ExecutingScript == "No" or ExecutingScript == "N" or ExecutingScript == "n":
        print("")
        print(f"\033[1;34m在{current_directory}/MissingMafHardyFilterFormat07/MissingMafHardyFilterFormat07.sh脚本，你可以先检查一下其内容是否正确，然后手动执行它! 也可以在原有命令上添加'-es y'参数直接执行。\033[0m")
        print("")
        pass
    elif ExecutingScript == "yes" or ExecutingScript == "YES" or ExecutingScript == "Yes" or ExecutingScript == "Y" or ExecutingScript == "y":
        os.system(f"bash {current_directory}/MissingMafHardyFilterFormat07/MissingMafHardyFilterFormat07.sh")
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
        "1.此步骤是vcf进行missing、maf、hwe、等位基因个数过滤。\n"

        "\033[0m")

    parser = argparse.ArgumentParser(description=description)

    # 接受的输入参数
    parser.add_argument("-rfn", metavar="ResultFileName", type=str, required=True,
                        help="最终结果文件的前缀名 (例如: AllSample_FinalResult)")
    parser.add_argument("-geno", metavar="geno", type=float, required=True,
                        help="缺失率超过多少的SNP将会被过滤掉 (例如: 0.1~0.2, 数值越小过滤掉的越多)")
    parser.add_argument("-maf", metavar="maf", type=float, required=True,
                        help="如果一个位点的最小等位基因频率小于多少将被过滤掉 (例如: 0.05~0.1, 一般0.05, 数值越大过滤掉的越多)")
    parser.add_argument("-mind", metavar="mind", type=float, default=1,
                        help="缺失率超过多少的样本将会被过滤掉 (默认: %(default)s, 意思是不过滤样本, 数值越小过滤掉的越多)")
    parser.add_argument("-hwe", metavar="hwe", type=float, default=1e-04,
                        help="哈温过滤,根据检验P值过滤，P值越小越不符合哈温平衡，要过滤掉。注意数值越大过滤掉的越多 (默认: %(default)s)")
    parser.add_argument("-min_ac", metavar="min_alleles", type=int, default=2,
                        help="最小等位基因个数,最小最大等位基因都设置为2意味着只要双等位基因的位点 (默认: %(default)s)")
    parser.add_argument("-max_ac", metavar="max_alleles", type=int, default=2,
                        help="最大等位基因个数,最小最大等位基因都设置为2意味着只要双等位基因的位点 (默认: %(default)s)")
    parser.add_argument("-t", metavar="thread", type=int, default=4,
                        help="使用多少个线程计算 (默认: %(default)s)")
    parser.add_argument("-format", metavar="format", type=str, default="vcf",
                        help="想要输出什么格式的结果, 可以选择vcf, bed, hapmap (默认: %(default)s)")
    parser.add_argument("-es", metavar="ExecutingScript", type=str, default="no",
                        help="是否直接执行脚本,选择yes或者no，默认是no (默认: %(default)s)")

    # 解析命令行参数
    args = parser.parse_args()

    # 调用主函数进行处理
    main(args.rfn, args.geno, args.maf, args.mind, args.hwe, args.min_ac, args.max_ac, args.t, args.format, args.es)
