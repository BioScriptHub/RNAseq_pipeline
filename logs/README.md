# logs

这里保存每一步运行日志。

流程失败时，先看终端中的 `[ERROR]` 和 `[HINT]`，再打开本目录中对应步骤的日志文件。

常见日志：

```text
00_unpack_rawdata.log
00_make_metadata.log
00_preflight_check.log
00_prepare_reference.log
01_fastp_qc.log
02_hisat2_align.log
03_bam_qc.log
04_featureCounts.log
05_DESeq2.log
```
