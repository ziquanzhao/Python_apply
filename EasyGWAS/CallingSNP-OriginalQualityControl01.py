import argparse
import os
import re


def main(InputFileDir, FilenameCharacteristics1, FilenameCharacteristics2, thread, QualifiedQualityPhred, UnqualifiedPercentLimit, LengthRequired, TrimFront, ExecutingScript):

    current_directory = os.getcwd()
    if not os.path.exists(f'{current_directory}/OriginalQualityControl01'):
        os.mkdir(f'{current_directory}/OriginalQualityControl01')

    with open(f'{current_directory}/OriginalQualityControl01/OriginalQualityControl01.sh', 'w') as trim:
        trim.write('#!/bin/bash\n')
        for filename in os.listdir(InputFileDir):
            if FilenameCharacteristics1 in filename:
                filename = re.sub(FilenameCharacteristics1, '', filename)
                trim.write(f'fastp -i {InputFileDir}{filename}{FilenameCharacteristics1} -o {current_directory}/OriginalQualityControl01/{filename}_CleanData{FilenameCharacteristics1} -I {InputFileDir}{filename}{FilenameCharacteristics2} -O {current_directory}/OriginalQualityControl01/{filename}_CleanData{FilenameCharacteristics2} '
                           f'-h {current_directory}/OriginalQualityControl01/{filename}_report.html -j {current_directory}/OriginalQualityControl01/{filename}_report.json '
                           f'--detect_adapter_for_pe -p -P 20 -q {QualifiedQualityPhred} -u {UnqualifiedPercentLimit} -l {LengthRequired} -f {TrimFront} -F {TrimFront} -c -w {thread}\n')

    if ExecutingScript == "no" or ExecutingScript == "NO" or ExecutingScript == "No" or ExecutingScript == "N" or ExecutingScript == "n":
        print("")
        print(f"\033[1;34m在{current_directory}/OriginalQualityControl01/目录下有个OriginalQualityControl01.sh脚本，你可以先检查一下其内容是否正确，然后手动执行它! 也可以在原有命令上添加'-es y'参数直接执行。\033[0m")
        print("")
        pass
    elif ExecutingScript == "yes" or ExecutingScript == "YES" or ExecutingScript == "Yes" or ExecutingScript == "Y" or ExecutingScript == "y":
        os.system(f"bash {current_directory}/OriginalQualityControl01/OriginalQualityControl01.sh")
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
        "1.此步骤是对原始测序数据进行质控，去掉低质量、短小的read。\n"
        "2.需要使用-i参数输入重测序文件所在路径，该脚本会根据-fc1/fc2参数输入的文件后缀自动检测该目录下重测序文件。\n"
        "3.需要输入重测序数据文件后缀命名特征，例如：_R1.fq.gz;_R2.fq.gz，这是为了正确地检测重测序文件的文件名。\n"
        "4._report.html文件是输出结果报告，可以查看，也可以根据这个报告重新设置参数质控。\n"
        "5.该脚本是调用fastp软件执行命令的，你需要安装好fastp软件，并设置好环境变量。"
        "\033[0m")

    parser = argparse.ArgumentParser(description=description)

    # 接受的输入参数
    parser.add_argument("-i", metavar="InputFileDir", type=str, required=True,
                        help="输入重测序数据文件所在文件夹的绝对路径 (例如: /mnt/d/cleandata/)")
    parser.add_argument("-fc1", metavar="FilenameCharacteristics1", type=str, default="_R1.fq.gz",
                        help="输入重测序数据文件名特征 (默认: %(default)s)")
    parser.add_argument("-fc2", metavar="FilenameCharacteristics2", type=str, default="_R2.fq.gz",
                        help="输入重测序数据文件名特征 (默认: %(default)s)")
    parser.add_argument("-q", metavar="QualifiedQualityPhred", type=int, default=25,
                        help="一个碱基被认为合格的质量值。默认20意味着phred质量>=Q20被认为是合格的 (默认: %(default)s)")
    parser.add_argument("-t", metavar="thread", type=int, default=4,
                        help="使用多少个线程计算 (默认: %(default)s)")
    parser.add_argument("-u", metavar="UnqualifiedPercentLimit", type=int, default=10,
                        help="允许多少百分比的碱基不合格（0~100）。默认10意味着10%% (默认: %(default)s)")
    parser.add_argument("-l", metavar="LengthRequired", type=int, default=50,
                        help="短于length_required的读取将被丢弃，默认为50 (默认: %(default)s)")
    parser.add_argument("-tf", metavar="TrimFront", type=int, default=0,
                        help="修剪read前端的碱基数，默认是0 (默认: %(default)s)")
    parser.add_argument("-es", metavar="ExecutingScript", type=str, default="no",
                        help="是否直接执行脚本,选择yes或者no，默认是no (默认: %(default)s)")

    # 解析命令行参数
    args = parser.parse_args()

    # 调用主函数进行处理
    main(args.i, args.fc1, args.fc2, args.t, args.q, args.u, args.l, args.tf, args.es)

"""
# --detect_adapter_for_pe 用于自动检测成对末端（paired-end, PE）测序数据中的接头序列

--overrepresentation_analysis 
-p 
作用: 该参数启用后，fastp 会在测序数据中检测和报告那些出现频率异常高的序列。这些过度表示的序列可能是由于技术性污染（如接头序列、引物二聚体等）或生物学上的原因（如高拷贝数的基因、重复序列等）所导致的。
分析结果: fastp 会生成一个报告，列出这些频繁出现的序列，并提供这些序列在总 reads 中的百分比。这有助于用户评估数据质量，并决定是否需要进一步处理这些过度表示的序列。

--overrepresentation_sampling 20
-P 20
fastp 中与过度表示序列分析（overrepresentation analysis）相关的一个参数，它用于控制在进行过度表示序列分析时采样的数据比例
作用: 这个参数指定了在进行过度表示序列分析时，fastp 从总的测序数据中进行采样的比例。由于对所有序列进行过度表示分析可能计算量较大，因此通过采样可以加快分析速度，同时仍然能获得可靠的结果。
默认值: 通常设置为 20，表示分析时只采样 20% 的数据

#如果是ChIP-seq测序加上下面列出的一些参数 4 20是一对宽松的参数，5 30是一对较为严格的参数
-5 --cut_front_window_size 4 --cut_front_mean_quality 20
-3 --cut_tail_window_size 4 --cut_tail_mean_quality 20
--cut_right --cut_right_window_size 4 --cut_right_mean_quality 20
"""

