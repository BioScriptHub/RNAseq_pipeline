# 00_rawdata

把最终用于分析的 FASTQ 文件放在这里。

常见格式：

```text
*_R1.fastq.gz
*_R2.fastq.gz
```

`.fastq.gz` 和 `.fq.gz` 是标准输入文件，不需要解压成 `.fastq`。

如果原始数据来自 `00_archives/` 中的交付压缩包，运行 `bash scripts/00_unpack_rawdata.sh` 后，FASTQ 会整理到本目录。
