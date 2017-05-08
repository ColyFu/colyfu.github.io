---
title: GWAS如何确定候选区间和causal gene
date: 2017-04-09 20:33:37
tags:
- GWAS
- causal gene
- LD
- heat map
categories:
- GWAS
---
随着高通量测序成本的急剧下降，越来越多的GWAS研究通过全基因组重测序获得标记。全基因组重测序的高密度标记，使得通过GWAS研究快速找到causal gene甚至causal mutation成为可能。但是通过GWAS得到一个信号后，怎么确定与这个信号关联的基因呢？  

<!--more-->

基本的思路就是找与这个信号强连锁的区域，一般来说r<sup>2</sup>大于0.6的区域视为强连锁的区域。简单粗糙的做法：我们可以通过群体的全基因组LD-decay分析，找到LD decay到r<sup>2</sup>等于0.6时所对应的距离，将GWAS超过阈值的信号前后各延伸这个距离作为候选区间。  
但是基因组不同区域的连锁程度差异很大，上面一刀切的做法可能使我们漏掉一些基因，或者多调查许多关联性并不强的基因。而且很多作物的连锁性很强，导致候选基因很多，这样就大大增加了工作难度。下面介绍一个非常简单高效的方法（[参考文献](http://www.nature.com/ng/journal/v48/n8/abs/ng.3596.html)）。  

<p align="center">
  <img src="http://onzjn6hm6.bkt.clouddn.com/20170409_heatmap.png" width="300"/>
</p>

1. 找到信号后，向前后延伸一段距离（可以根据全基因组的LD-decay水平大概估计），计算这段区域内所有标记pairwise r<sup>2</sup>，将r<sup>2</sup>大于0.6的block作为候选区间。pairwise r<sup>2</sup>可以用[PLINK](http://zzz.bwh.harvard.edu/plink/)计算：

 ```bash
 plink --noweb --bfile <bfile_prefix> \
       --chr 5 --from-bp 13641890 --to-bp 17641890 \
       --matrix --r2 --out <out_prefix>
 ```

 画图用R,输入文件为plink计算的到的r<sup>2</sup>矩阵，以及标记的位置:

 ```R
 #!/usr/bin/env Rscript
 library(LDheatmap)
 argv <- commandArgs(TRUE)
 ldmatrix <- as.matrix(read.table(argv[1],sep=' '))
 pos <- as.numeric(unlist(read.table(argv[2], head=FALSE)))
 pdf(argv[3])
 rgb.palette <- colorRampPalette(rev(c("yellow", "orange", "red")), space = "rgb")
 LDheatmap(ldmatrix, genetic.distances=pos, color=rgb.palette(100), flip=TRUE)
 dev.off()
```

2. 将这个区域内的标记按照其对基因功能的影响程度分为5类：  
 > 1） 标记与性状显著关联（-log<sub>10</sub>P大于阈值），且该标记影响氨基酸编码，或者位于剪接位点；  
 2） 标记与性状显著关联，且位于基因起始密码子上游2 kb内；  
 3） 标记与性状显著关联，且位于基因内，除开1）和2）之外的标记  
 4） 标记与性状显著关联，位于基因间区  
 5） 标记与性状不显著关联  
 
然后，按照这5类的顺序，依次调查，一般来说，属于1）类的可能性很大，而且基因一般就几个，这样就大大减少了工作难度。
