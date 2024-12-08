import argparse
import os
import re


def main(ResequencingFileDir, FilenameCharacteristics1, FilenameCharacteristics2, RefGenome, thread, TmpDir, ValidationStringency1, SortOrder, MaxFileHandles, RemoveDuplicates, ValidationStringency2, ExecutingScript):

    current_directory = os.getcwd()
    if not os.path.exists(f'{current_directory}/SamToBamSortDuplicateIndex03'):
        os.mkdir(f'{current_directory}/SamToBamSortDuplicateIndex03')

    with open(f'{current_directory}/SamToBamSortDuplicateIndex03/SamToBamSortDuplicateIndex03.sh', 'w') as trans_sort:
        trans_sort.write('#!/bin/bash\n')
        for filename in os.listdir(f'{ResequencingFileDir}'):
            if FilenameCharacteristics1 in filename:
                filename = re.sub(f'{FilenameCharacteristics1}', '', filename)
                trans_sort.write(
                    f'bwa mem -t {thread} -M -R "@RG\tID:{filename}\tLB:{filename}\tPL:illumina\tPU:{filename}\tSM:{filename}\" {RefGenome} {ResequencingFileDir}{filename}{FilenameCharacteristics1} {ResequencingFileDir}{filename}{FilenameCharacteristics2} > {current_directory}/SamToBamSortDuplicateIndex03/{filename}.sam' + "\n")
                trans_sort.write(
                    f"PicardCommandLine SortSam -TMP_DIR {TmpDir} -VALIDATION_STRINGENCY {ValidationStringency1} -INPUT {current_directory}/SamToBamSortDuplicateIndex03/{filename}.sam -OUTPUT {current_directory}/SamToBamSortDuplicateIndex03/{filename}_sort.bam -SORT_ORDER {SortOrder}\n")
                trans_sort.write(
                    f'PicardCommandLine MarkDuplicates -MAX_FILE_HANDLES_FOR_READ_ENDS_MAP {MaxFileHandles} -REMOVE_DUPLICATES {RemoveDuplicates} -INPUT {current_directory}/SamToBamSortDuplicateIndex03/{filename}_sort.bam -OUTPUT {current_directory}/SamToBamSortDuplicateIndex03/{filename}_SortDup.bam -METRICS_FILE {current_directory}/SamToBamSortDuplicateIndex03/{filename}_SortDup.metrics -VALIDATION_STRINGENCY {ValidationStringency2}\n')
                trans_sort.write(f'samtools index -b {current_directory}/SamToBamSortDuplicateIndex03/{filename}_SortDup.bam > {current_directory}/SamToBamSortDuplicateIndex03/{filename}_SortDup.bam.bai\n')
                trans_sort.write(f'rm {current_directory}/SamToBamSortDuplicateIndex03/{filename}.sam {current_directory}/SamToBamSortDuplicateIndex03/{filename}_sort.bam\n')
    os.system(f"sed -i 's/\\t/\\\\t/g' {current_directory}/SamToBamSortDuplicateIndex03/SamToBamSortDuplicateIndex03.sh")

    if ExecutingScript == "no" or ExecutingScript == "NO" or ExecutingScript == "No" or ExecutingScript == "N" or ExecutingScript == "n":
        print("")
        print(f"\033[1;34m在{current_directory}/SamToBamSortDuplicateIndex03/目录下有个SamToBamSortDuplicateIndex03.sh脚本，你可以先检查一下其内容是否正确，然后手动执行它! 也可以在原有命令上添加'-es y'参数直接执行。\033[0m")
        print("")
        pass
    elif ExecutingScript == "yes" or ExecutingScript == "YES" or ExecutingScript == "Yes" or ExecutingScript == "Y" or ExecutingScript == "y":
        os.system(f"bash {current_directory}/SamToBamSortDuplicateIndex03/SamToBamSortDuplicateIndex03.sh")
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
        "1.此步是对质控后的重测序数据进行比对、排序、去重和建立索引文件。\n"
        "2.需要使用-rfd参数输入经过原始质控后的重测序数据文件所在文件夹的绝对路径，该脚本会根据-fc1/fc2参数输入的文件后缀自动检测该目录下重测序文件。\n"
        "3.需要使用-fc1/-fc2输入经过原始质控后的重测序数据文件的后缀\n"
        "4.需要使用-rg参数输入参考基因组及其索引所在绝对路径，索引包含.dict/.fai/.amb/.ann/.bwt/.pac/.sa\n"
        "5.该脚本是调用bwa、picard、samtools软件执行命令的，你需要安装好这些软件，并设置好环境变量。"
        "\033[0m")

    parser = argparse.ArgumentParser(description=description)

    # 接受的输入参数
    parser.add_argument("-rfd", metavar="ResequencingFileDir", type=str, required=True,
                        help="输入经过原始质控后的重测序数据文件所在文件夹的绝对路径 (例如: /mnt/d/GWAS/OriginalQualityControl01/)")
    parser.add_argument("-fc1", metavar="FilenameCharacteristics1", type=str, required=True,
                        help="输入经过原始质控后的重测序数据文件名特征 (例如: _CleanData.R1.fq.gz)")
    parser.add_argument("-fc2", metavar="FilenameCharacteristics2", type=str, required=True,
                        help="输入经过原始质控后的重测序数据文件名特征 (例如: _CleanData.R2.fq.gz)")
    parser.add_argument("-rg", metavar="RefGenome", type=str, required=True,
                        help="输入参考基因组的绝对路径 (例如：/mnt/d/GWAS/BuildRefseqIndex02/YourReferenceGenome.fa)")
    parser.add_argument("-t", metavar="thread", type=int, default=4,
                        help="使用多少个线程计算 (默认: %(default)s)")
    parser.add_argument("-td", metavar="TmpDir", type=str, default="./",
                        help="设置临时文件存放路径 (默认: %(default)s)")
    parser.add_argument("-vs1", metavar="ValidationStringency1", type=str, default="SILENT",
                        help="检验Sam文件的格式正确与否，设置为'SILENT'可提高性能 (默认: %(default)s)")
    parser.add_argument("-so", metavar="SortOrder", type=str, default="coordinate",
                        help="输出文件的排序方式, coordinate的意思是按照染色体上的顺序排序, (默认: %(default)s)")
    parser.add_argument("-mfh", metavar="MaxFileHandles", type=int, default=800,
                        help="将读取结束溢出到磁盘时保持打开状态的最大文件句柄数 (默认: %(default)s)")
    parser.add_argument("-rd", metavar="RemoveDuplicates", type=str, default="true",
                        help="是否丢弃duplicated序列, 对结果应该没有什么影响 (默认: %(default)s)")
    parser.add_argument("-vs2", metavar="ValidationStringency2", type=str, default="LENIENT",
                        help="在BWA比对时,没有map到基因组上的read归到了ref以外的区域，Picard认为这些read是不应该出现的,会报错. 如果想忽略报错的话,就使用这行代码 (默认: %(default)s)")
    parser.add_argument("-es", metavar="ExecutingScript", type=str, default="no",
                        help="是否直接执行脚本,选择yes或者no, 默认是no (默认: %(default)s)")

    # 解析命令行参数
    args = parser.parse_args()

    # 调用主函数进行处理
    main(args.rfd, args.fc1, args.fc2, args.rg, args.t, args.td, args.vs1, args.so, args.mfh, args.rd, args.vs2, args.es)
