# Statistical-GeSeq-software-annotation-results
统计分析GeSeq软件注释结果

usage: statistic.py [-h] [-f F] [-a A] [--out OUT]

 
name:statistic.py -- 统计分析GeSeq软件注释结果
attention: python statistic.py -f genomic.fasta -a genomic.gff3 -o out
version: 1.0.0
email:xby@bioyigene.com 

optional arguments:
 
 -h, --help         show this help message and exit
 
 -f F               fasta file
 
 -a A               Annotation file
 
 --out OUT, -o OUT  out put file

结果：
|class|total|number|average|percentage|
|-----|-----|-----|-----|-----|
|gene|92390|84|1099.88|54.36%|
|CDS|77934|100|779.34|45.86%|
|rRNA|9052|8|1131.5|5.33%|
|tRNA|10667|37|288.3|6.28%|
