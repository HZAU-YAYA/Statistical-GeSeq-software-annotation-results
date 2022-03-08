# -*- coding: utf-8 -*-
import re
import logging
import argparse
import sys


LOG = logging.getLogger(__name__)
__version__ = "1.0.1"    #设置版本信息
__author__ = ("Boya Xu",)   #输入作者信息
__email__ = "xby@bioyigene.com"
__all__ = []


def add_help_args(parser): #帮助函数
    parser.add_argument('-f',  default="", type=str, help='fasta file')
    parser.add_argument('-a',  default="", type=str, help='Annotation file')
    parser.add_argument('--out', '-o', type=str, default="", help="out put file")
    return parser


def process_file(reader):
    '''Open, read,and print a file'''
    names=[]
    index=0
    dict={}

    for line in reader:
        if line.startswith('>'):
           if index >=1:
               names.append(line)
           index =index+1
           name=line[:-1]
           seq = ''
        else:
           seq +=line[:-1]
           dict[name]=seq
    return dict


def statistic(gff, genomic, outfile):
    #统计Gene CDS rRNA tRNA信息的主函数
    sum_g, num_g, sum_r, sum_c, num_r, sum_t, num_t, num_c = 0, 0, 0, 0, 0, 0, 0, 0
    length_g, length_r, length_t, length_c, length = 0, 0, 0, 0, 0
    m = 1
    input_file = open(genomic, "r")
    reader = input_file.readlines()
    items = process_file(reader)
    for key in items:
        length = int(len(items[key]))
    with open(gff) as file_object:
        for line in file_object:
            if line.startswith('#'):
                if m == 0:
                    num_g -= 1
                    sum_g -= length_g
                else:
                    m = 0
                continue
            else:
                w_1 = line
                p1 = re.compile(r'[a-zA-Z]+')
                w_1_1 = p1.findall(w_1)
                p2 = re.compile(r'\d+')
                w_1_2 = p2.findall(w_1)

                if w_1_1[2] == 'CDS':
                    length_c = int(w_1_2[2]) - int(w_1_2[1]) + 1
                    sum_c += length_c
                    num_c += 1
                    m += 1
                elif w_1_1[2] == 'gene':
                    length_g = int(w_1_2[2]) - int(w_1_2[1]) + 1
                    sum_g += length_g
                    num_g += 1
                elif w_1_1[2] == 'rRNA':
                    length_r = int(w_1_2[2]) - int(w_1_2[1]) + 1
                    sum_r += length_r
                    num_r += 1
                elif w_1_1[2] == 'tRNA':
                    length_t = int(w_1_2[2]) - int(w_1_2[1]) + 1
                    sum_t += length_t
                    num_t += 1
                else:
                    continue
    len_g_p = round((sum_g * 1.00) / num_g, 2)
    len_r_p = round((sum_r * 1.00) / num_r, 2)
    len_t_p = round((sum_t * 1.00) / num_t, 2)
    len_c_p = round((sum_c * 1.00) / num_c, 2)
    percentage_g = round((sum_g * 100.00)/length, 2)
    percentage_t = round((sum_t * 100.00)/length, 2)
    percentage_r = round((sum_r * 100.00)/length, 2)
    percentage_c = round((sum_c * 100.00)/length, 2)
    site_g = {"gene": "gene\t{}\t{}\t{}\t{}%".format(sum_g, num_g, len_g_p, percentage_g)}
    site_c = {"CDS": "CDS\t{}\t{}\t{}\t{}%".format(sum_c, num_c, len_c_p, percentage_c)}
    site_r = {"rRNA": "rRNA\t{}\t{}\t{}\t{}%".format(sum_r, num_r, len_r_p, percentage_r)}
    site_t = {"tRNA": "tRNA\t{}\t{}\t{}\t{}%".format(sum_t, num_t, len_t_p, percentage_t)}
    f = open(outfile, 'w+')
    output = 'class\ttotal\tnumber\taverage\tpercentage\n{}\n{}\n{}\n{}\n'
    output = output.format(site_g['gene'], site_c['CDS'], site_r['rRNA'], site_t['tRNA'])
    f.write(output)
    f.close()
    return
 
    
def main():   #主函数，执行函数
    logging.basicConfig(stream=sys.stderr, level=logging.INFO, format="[%(levelname)s] %(message)s")
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=''' 
name:statistic.py -- 统计分析GeSeq软件注释结果
attention: python statistic.py -f genomic.fasta -a genomic.gff3 -o out
version: %s
contact: %s <%s>\ 
''' % (__version__, ' '.join(__author__), __email__))
    args = add_help_args(parser).parse_args()
    statistic(args.a, args.f, args.out)


if __name__ == "__main__":           #固定格式，使 import 到其他的 python 脚本中被调用（模块重用）执行
    main()
