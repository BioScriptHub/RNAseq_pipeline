# 04_reference

把参考基因组和注释文件放在这里。

不需要改成固定文件名。脚本按后缀自动识别。

```text
基因组：*.fa、*.fasta、*.fna、*.fas
注释：*.gtf、*.gff3、*.gff
```

NCBI 下载的文件通常可以直接放进来，例如：

```text
GCF_000000000.1_genomic.fna
genomic.gff
```

要求：

- 同一目录中只能放一个基因组 FASTA 和一个注释文件
- 基因组和注释来自同一参考版本
- 参考和注释不匹配会导致比对率低或 reads 分配率低
