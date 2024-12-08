import argparse
import os
import re


def main(snp_QD_filter, snp_QUAL_filter, snp_SOR_filter, snp_FS_filter, snp_MQ_filter, indel_QD_filter, indel_SOR_filter, indel_FS_filter, ExecutingScript):

    current_directory = os.getcwd()
    if not os.path.exists(f'{current_directory}/MergeHardFilter06'):
        os.mkdir(f'{current_directory}/MergeHardFilter06')

    with open(f'{current_directory}/MergeHardFilter06/MergeHardFilter06.sh', 'w') as mergehf:
        mergehf.write('#!/bin/bash\n')
        vcfHFlist = []
        for filename in os.listdir(f'{current_directory}/GvcfToVcfSplitSnpIndel05/'):
            if '_snp.vcf' in filename and '.idx' not in filename:
                filename = re.sub('_snp.vcf', '', filename)
                mergehf.write(f'gatk VariantFiltration -V {current_directory}/GvcfToVcfSplitSnpIndel05/{filename}_snp.vcf -R {current_directory}/BuildRefseqIndex02/YourReferenceGenome.fa '
                                 f'-filter "QD < {snp_QD_filter}" --filter-name "QD{snp_QD_filter}" '
                                 f'-filter "QUAL < {snp_QUAL_filter}" --filter-name "QUAL{snp_QUAL_filter}" '
                                 f'-filter "SOR > {snp_SOR_filter}" --filter-name "SOR{snp_SOR_filter}" '
                                 f'-filter "FS > {snp_FS_filter}" --filter-name "FS{snp_FS_filter}" '
                                 f'-filter "MQ < {snp_MQ_filter}" --filter-name "MQ{snp_MQ_filter}" '
                                 f'--cluster-window-size 10  --cluster-size 3 --missing-values-evaluate-as-failing -O {current_directory}/MergeHardFilter06/{filename}_snp_filter.vcf\n')
                mergehf.write(f'grep -E \'#|PASS\' {current_directory}/MergeHardFilter06/{filename}_snp_filter.vcf > {current_directory}/MergeHardFilter06/{filename}_snp_pass.vcf\n')
                mergehf.write(f'gatk VariantFiltration -V {current_directory}/GvcfToVcfSplitSnpIndel05/{filename}_indel.vcf -R {current_directory}/BuildRefseqIndex02/YourReferenceGenome.fa '
                                 f'-filter "QD < {indel_QD_filter}" --filter-name "QD{indel_QD_filter}" '
                                 f'-filter "SOR > {indel_SOR_filter}" --filter-name "SOR{indel_SOR_filter}" '
                                 f'-filter "FS > {indel_FS_filter}" --filter-name "FS{indel_FS_filter}" '
                                 f'--missing-values-evaluate-as-failing -O {current_directory}/MergeHardFilter06/{filename}_indel_filter.vcf\n')
                mergehf.write(f'grep -E \'#|PASS\' {current_directory}/MergeHardFilter06/{filename}_indel_filter.vcf > {current_directory}/MergeHardFilter06/{filename}_indel_pass.vcf\n')
                mergehf.write(f'gatk MergeVcfs -I {current_directory}/MergeHardFilter06/{filename}_snp_pass.vcf -I {current_directory}/MergeHardFilter06/{filename}_indel_pass.vcf -O {current_directory}/MergeHardFilter06/{filename}_HardFilter.vcf\n')
                vcfHFlist.append(f'{current_directory}/MergeHardFilter06/{filename}_HardFilter.vcf')
        vcfHFlist = ' '.join(vcfHFlist)
        mergehf.write(f'bcftools concat {vcfHFlist} -o {current_directory}/MergeHardFilter06/AllSample_HardFilter.vcf.gz\n')
        mergehf.write(f'rm {current_directory}/MergeHardFilter06/*.vcf.idx {current_directory}/MergeHardFilter06/*_HardFilter.vcf {current_directory}/MergeHardFilter06/*_snp_filter.vcf {current_directory}/MergeHardFilter06/*_indel_filter.vcf {current_directory}/MergeHardFilter06/*_snp_pass.vcf {current_directory}/MergeHardFilter06/*_indel_pass.vcf\n')
        mergehf.write(f'plink1.9 --vcf {current_directory}/MergeHardFilter06/AllSample_HardFilter.vcf.gz --allow-extra-chr --make-bed --max-ac 2 --min-ac 2 --out {current_directory}/MergeHardFilter06/AllSample_HF_statistics --threads 4\n')
        mergehf.write(f'plink1.9 -bfile {current_directory}/MergeHardFilter06/AllSample_HF_statistics --missing --out {current_directory}/MergeHardFilter06/AllSample_HF_statistics_miss --allow-extra-chr --threads 4\n')
        mergehf.write(f'plink1.9 -bfile {current_directory}/MergeHardFilter06/AllSample_HF_statistics --freq --out {current_directory}/MergeHardFilter06/AllSample_HF_statistics_freq --allow-extra-chr --threads 4\n')
        mergehf.write(f'plink1.9 -bfile {current_directory}/MergeHardFilter06/AllSample_HF_statistics --hardy --out {current_directory}/MergeHardFilter06/AllSample_HF_statistics_hardy --allow-extra-chr --threads 4\n')
        mergehf.write(f'rm {current_directory}/MergeHardFilter06/AllSample_HF_statistics.* {current_directory}/MergeHardFilter06/*.log {current_directory}/MergeHardFilter06/*.nosex')




    if ExecutingScript == "no" or ExecutingScript == "NO" or ExecutingScript == "No" or ExecutingScript == "N" or ExecutingScript == "n":
        print("")
        print(f"\033[1;34m在{current_directory}/OriginalQualityControl01/目录下有个OriginalQualityControl01.sh脚本，你可以先检查一下其内容是否正确，然后手动执行它! 也可以在原有命令上添加'-es y'参数直接执行。\033[0m")
        print("")
        pass
    elif ExecutingScript == "yes" or ExecutingScript == "YES" or ExecutingScript == "Yes" or ExecutingScript == "Y" or ExecutingScript == "y":
        os.system(f"bash {current_directory}/MergeHardFilter06/MergeHardFilter06.sh")
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
        "1.此步骤是对snp和indel进行硬过滤,然后将所有vcf文件合并起来，最后统计一下Maf、Missing、hardy,为下一步缺失过滤做好准备。\n"
        "2..lmiss是基因型缺失统计(第5列),.imiss是样本缺失统计(第5列),.frq是等位基因频率统计(第5列),.hwe是哈温平衡统计(第9列)。\n"
        "\033[0m")

    parser = argparse.ArgumentParser(description=description)

    # 接受的输入参数
    parser.add_argument("-qd_snp", metavar="snp_QD_filter", type=int, default=2.0)
    parser.add_argument("-qual_snp", metavar="snp_QUAL_filter", type=int, default=30.0)
    parser.add_argument("-sor_snp", metavar="snp_SOR_filter", type=int, default=3.0)
    parser.add_argument("-fs_snp", metavar="snp_FS_filter", type=int, default=60.0)
    parser.add_argument("-mq_snp", metavar="snp_MQ_filter", type=int, default=40.0)
    parser.add_argument("-qd_indel", metavar="indel_QD_filter", type=int, default=2.0)
    parser.add_argument("-sor_indel", metavar="indel_SOR_filter", type=int, default=10.0)
    parser.add_argument("-fs_indel", metavar="indel_FS_filter", type=int, default=200.0)
    parser.add_argument("-es", metavar="ExecutingScript", type=str, default="no",
                        help="是否直接执行脚本,选择yes或者no，默认是no (默认: %(default)s)")

    # 解析命令行参数
    args = parser.parse_args()

    # 调用主函数进行处理
    main(args.qd_snp, args.qual_snp, args.sor_snp, args.fs_snp, args.mq_snp, args.qd_indel, args.sor_indel, args.fs_indel, args.es)
