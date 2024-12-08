import argparse
import os
import re


def main(ChromosomeNumber, ChromosomeLabel, ExecutingScript):

    current_directory = os.getcwd()
    if not os.path.exists(f'{current_directory}/GvcfToVcfSplitSnpIndel05'):
        os.mkdir(f'{current_directory}/GvcfToVcfSplitSnpIndel05')

    with open(f'{current_directory}/GvcfToVcfSplitSnpIndel05/GvcfToVcfSplitSnpIndel05.sh', 'w') as gvcftovcf:
        gvcftovcf.write('#!/bin/bash\n')
        for chr_num in range(1, ChromosomeNumber + 1, 1):
            lachesis_group = []
            for filename in os.listdir(f'{current_directory}/Haplotypecaller04/'):
                if '.idx' not in filename and ('_chr1.g.vcf') in filename:
                    filename = re.sub('_chr1.g.vcf', '', filename)
                    lachesis_group.append(f'-V {current_directory}/Haplotypecaller04/{filename}_chr{chr_num}.g.vcf')
            lachesis_group = ' '.join(lachesis_group)
            gvcftovcf.write(
                f'gatk CombineGVCFs -R {current_directory}/BuildRefseqIndex02/YourReferenceGenome.fa {lachesis_group} -O {current_directory}/GvcfToVcfSplitSnpIndel05/Chr{chr_num}.g.vcf\n')
            gvcftovcf.write(
                f'gatk GenotypeGVCFs -R {current_directory}/BuildRefseqIndex02/YourReferenceGenome.fa -V {current_directory}/GvcfToVcfSplitSnpIndel05/Chr{chr_num}.g.vcf -O {current_directory}/GvcfToVcfSplitSnpIndel05/Chr{chr_num}.vcf\n')
            gvcftovcf.write(f'gatk SelectVariants -R {current_directory}/BuildRefseqIndex02/YourReferenceGenome.fa -V {current_directory}/GvcfToVcfSplitSnpIndel05/Chr{chr_num}.vcf -O {current_directory}/GvcfToVcfSplitSnpIndel05/Chr{chr_num}_snp.vcf -L {ChromosomeLabel}{chr_num} --select-type-to-include SNP\n')
            gvcftovcf.write(f'gatk SelectVariants -R {current_directory}/BuildRefseqIndex02/YourReferenceGenome.fa -V {current_directory}/GvcfToVcfSplitSnpIndel05/Chr{chr_num}.vcf -O {current_directory}/GvcfToVcfSplitSnpIndel05/Chr{chr_num}_indel.vcf -L {ChromosomeLabel}{chr_num} --select-type-to-include INDEL\n')
            gvcftovcf.write(f'rm {current_directory}/GvcfToVcfSplitSnpIndel05/Chr{chr_num}.g.vcf {current_directory}/GvcfToVcfSplitSnpIndel05/Chr{chr_num}.g.vcf.idx {current_directory}/GvcfToVcfSplitSnpIndel05/Chr{chr_num}.vcf {current_directory}/GvcfToVcfSplitSnpIndel05/Chr{chr_num}.vcf.idx\n')
            gvcftovcf.write(f'gatk VariantsToTable -V {current_directory}/GvcfToVcfSplitSnpIndel05/Chr{chr_num}_snp.vcf -F QD -F QUAL -F SOR -F FS -F MQ -O {current_directory}/GvcfToVcfSplitSnpIndel05/Chr{chr_num}_snp.recode.table\n')
            gvcftovcf.write(f'gatk VariantsToTable -V {current_directory}/GvcfToVcfSplitSnpIndel05/Chr{chr_num}_indel.vcf -F QD -F SOR -F FS -O {current_directory}/GvcfToVcfSplitSnpIndel05/Chr{chr_num}_indel.recode.table\n')



    if ExecutingScript == "no" or ExecutingScript == "NO" or ExecutingScript == "No" or ExecutingScript == "N" or ExecutingScript == "n":
        print("")
        print(f"\033[1;34m在{current_directory}/OriginalQualityControl01/目录下有个OriginalQualityControl01.sh脚本，你可以先检查一下其内容是否正确，然后手动执行它! 也可以在原有命令上添加'-es y'参数直接执行。\033[0m")
        print("")
        pass
    elif ExecutingScript == "yes" or ExecutingScript == "YES" or ExecutingScript == "Yes" or ExecutingScript == "Y" or ExecutingScript == "y":
        os.system(f"bash {current_directory}/GvcfToVcfSplitSnpIndel05/GvcfToVcfSplitSnpIndel05.sh")
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
        "1.此步骤是将每个样本的gvcf文件转换为vcf文件，并按照染色体合并，最后将SNP和Indel分开。\n"
        "2.需要使用-cn参数输入参考基因组中染色体条数。\n"
        "3.snp.recode.table和indel.recode.table文件是snp和indel硬过滤参数统计文件，可以先用R绘图，看看统计结果，然后再进行下一步硬过滤"
        "\033[0m")

    parser = argparse.ArgumentParser(description=description)

    # 接受的输入参数
    parser.add_argument("-cn", metavar="ChromosomeNumber", type=int, required=True,
                        help="输入你的参考基因组中的染色体条数 (例如：15)")
    parser.add_argument("-cl", metavar="ChromosomeLabel", type=str, required=True,
                        help="输入你的参考基因组中的染色体标签 (例如基因组中染色体以Chr1、Chr2、Chr3命名，则输入：Chr)")
    parser.add_argument("-es", metavar="ExecutingScript", type=str, default="no",
                        help="是否直接执行脚本,选择yes或者no，默认是no (默认: %(default)s)")

    # 解析命令行参数
    args = parser.parse_args()

    # 调用主函数进行处理
    main(args.cn, args.cl, args.es)
