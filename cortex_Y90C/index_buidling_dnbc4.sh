#!/bin/bash
# 1. 在宿主机上先创建用于存放华大新索引的隔离文件夹
mkdir -p /mnt/windowsdata/qiuzerui/scannotations/mouse/mgi_mm10_v3_index

# 2. 使用 Apptainer 穿透运行容器内的建索引命令
# ⚠️ 必须包含 --bind /mnt:/mnt，且请确保 dnbc4tools.sif 在你当前的操作目录下（否则需写绝对路径）
apptainer exec --bind /mnt:/mnt /mnt/disk1/qiuzerui/apptainer/dnbc4tools.sif dnbc4tools rna mkref \
    --fasta /mnt/windowsdata/qiuzerui/scannotations/mouse/refdata-gex-mm10-2020-A/fasta/genome.fa \
    --ingtf /mnt/windowsdata/qiuzerui/scannotations/mouse/refdata-gex-mm10-2020-A/genes/genes.gtf \
    --genomeDir /mnt/windowsdata/qiuzerui/scannotations/mouse/mgi_mm10_v3_index \
    --species mm10
    --threads 128