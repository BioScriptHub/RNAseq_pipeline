# 09_metadata

这里保存样本分组和差异比较信息。

默认只需要编辑：

```text
contrasts.csv
```

创建方法：

```bash
cp 09_metadata/contrasts.example.csv 09_metadata/contrasts.csv
```

`metadata.tsv` 会从 `00_rawdata/` 中的 FASTQ 文件名自动生成。

推荐 FASTQ 命名：

```text
Control_1_R1.fastq.gz
Control_1_R2.fastq.gz
Treatment_1_R1.fastq.gz
Treatment_1_R2.fastq.gz
```

如果文件名无法自动推断分组，再创建兜底表：

```bash
cp 09_metadata/samples.example.csv 09_metadata/samples.csv
```

手动生成：

```bash
bash scripts/00_make_metadata.sh
```

脚本会自动生成程序使用的 TSV：

```text
metadata.tsv
contrasts.tsv
```

不建议新手手动编辑 TSV，因为 tab 容易和空格混淆。
