#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@Email       :   h2624366594@163.com
@Author      :   Fuchuan Han 
@Time        :   2023/03/21 16:30:49
"""

from Bio import SeqIO

record = SeqIO.read("/home/hanfc/pub/Linux_001/workspace/mitogenome/Salix_wilsonii/Salix_wilsonii.gb", "genbank")

# 提取蛋白编码基因的信息
CDS_features = [feature for feature in record.features if feature.type == "CDS"]
rRNA_features = [feature for feature in record.features if feature.type == "rRNA"]
tRNA_features = [feature for feature in record.features if feature.type == "tRNA"]
intron_features = [feature for feature in record.features if feature.type == "intron"]

# 计算每个特征的总长度，并将它们相加
def get_total_length(features):
    total_length = 0
    for feature in features:
        length = len(feature.location)
        total_length += length
    return total_length

total_coding_length = get_total_length(CDS_features)
total_rRNA_length = get_total_length(rRNA_features)
total_tRNA_length = get_total_length(tRNA_features)

# 计算基因组总长度
total_genome_length = len(record.seq)

# 计算各个特征所占的百分比
def get_percentage(total_length):
    return total_length / total_genome_length * 100

coding_percentage = get_percentage(total_coding_length)
rRNA_percentage = get_percentage(total_rRNA_length)
tRNA_percentage = get_percentage(total_tRNA_length)

# 输出结果
print("Total length of protein-coding genes:", total_coding_length, f"({coding_percentage:.2f}%)")
print("Total length of rRNA genes:", total_rRNA_length, f"({rRNA_percentage:.2f}%)")
print("Total length of tRNA genes:", total_tRNA_length, f"({tRNA_percentage:.2f}%)")

# 计算内含子总长度
total_intron_length = get_total_length(intron_features)

# 输出结果
print("Total intron length:", total_intron_length)
