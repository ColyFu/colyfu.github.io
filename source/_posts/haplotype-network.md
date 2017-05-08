---
title: 单体型网络（Haplotype Network）构建及展示
date: 2017-05-08 00:24:35
tags:
- Haplotype Network
- 单体型网络
- PopART
- 谱系地理
categories:
- 群体遗传及进化
---

单体型网络（Haplotype Network）是谱系地理研究的重要手段。通过单体型网络，我们可以推断群体的起源、扩散历史。那么怎样理解单体型网络以及怎样构建单体型网络？

### 什么是单体型网络？
单体型（haplotype）在单体型网络中是指一段遗传连锁的核酸序列。不同的单体型，通过序列中的变异来区分（常用SNP）。一般用来进行单体型网络构建的单体型有：1）线粒体基因组；2）Y染色体；3）叶绿体基因组；4）基因组上一段强连锁的区段。  

<p align="center">
  <img src="http://onzjn6hm6.bkt.clouddn.com/haplotype_network.png" width="300"/>
</p>

<!--more-->

上图是一个单体型网络的例子，图中一个圆圈表示一个单体型，两个圆圈之间的连线表示这两个单体型相关（一个是由另一个突变而来），连线上面的短竖线表示从一个单体型到与其相连的单体型需要经历的碱基替换数，一个竖线表示一个替换。彩色的圆圈表示我们实际取样到的单体型，圆圈大小表示这种单体型的个数。灰色圆圈表示推断出来可能存在的中间单体型，没有被取样到。一种颜色一般表示一个群体，如按地理划分，品种划分等。图中例子一个单体型只存在于一个群体中，实际情况一个单体型往往在多个群体中出现。此时，一个单体型圆圈中填充多种颜色，以饼图的形式展示。  
从上图我们可以猜测一种可能的群体历史：A群体和C群体都起源于B群体，A群体从B群体分化出来后，经历过急剧的群体扩张，导致A群体的单体型演化出多个亚型。当然，为了避免过度解读，推断群体历史需要多种证据结合起来。
### 构建单体型网络的工具

比较经典的软件就是[Network](http://www.fluxus-engineering.com/sharenet.htm)了，但是该软件只支持Windows系统，用起来也比较复杂。这里推荐[PopART](http://popart.otago.ac.nz/downloads.shtml)，该软件支持Windows，Mac，Linux系统，而且用起来也非常方便，支持多种常用的Network构建方法，关键是该软件支持地图的形式展示单体型分布。下面简单介绍一下该软件的使用方法。
PopART的输入文件格式为NEXUS，一般主要用到两个部分DATA和TRAITS。

```
#NEXUS

begin data;
dimensions ntax=4 nchar=30;
format datatype=dna missing=N gap=-;
matrix
seq1 CCACCGTTGCTAAAAATTCATGACACAAGG
seq2 CCACAGTTTCTAAAAATTCGTGATACAAGG
seq3 CCACAGTTGCTACAAATTCATGATACAAGG
seq4 CCACAGGTGCTAAAAATTCATGAAACAAGG
;
end;

BEGIN TRAITS;
    Dimensions NTRAITS=5;
    Format labels=yes missing=? separator=Comma;
    TraitLatitude 53 43.6811 5.4 -25.61 -0;
    TraitLongitude 16.75 87.3311 26.5 134.355 -76;
    TraitLabels Europe Asia Africa Australia America;
    Matrix
    seq1 10,5,0,6,0
    seq2 0,0,5,0,0
    seq3 4,0,10,0,0
    seq4 0,0,0,4,2
    ;
END;
```

DATA部分主要纪录单体型信息，比较好理解。TRAITS部分主要纪录单体型来源的群体。如上所示，例子中取了来自5个大洲的样本，一共4种单体型，TRAITS纪录了每种单体型在不同大洲取样的个数，如`seq1`在`Europe`有10个，在`Asia`有5个等等。关键字TraitLatitude和TraitLongitude纪录5个群体取样地点的经纬度，该信息在单体型网络构建中可以不用，当需要用地图展示单体型地理分布时，需要填该信息。NEXUS文件生成后，打开PopART，通过`File -> Open`输入NEXUS文件，然后通过菜单栏`Network`选择单体型网络构建算法，如常用的Median Joining Network。选择Median Joining Network后，会提示填写Epsilon参数，该参数用来控制推断中间单体型的细节程度，该值越大，会展示更多推断的中间单体型，一般选择默认的0就好。填好该参数后，点击`OK`，就生成了我们需要的单体型网络。然后通过菜单栏`Edit`下的选项，对图中群体的颜色、字体、图例等进行调整。  
下面介绍一下PopART的单体型地理分布展示。

<p align="center">
  <img src="http://onzjn6hm6.bkt.clouddn.com/network_map.png" width="300"/>
</p>
点击`View -> Switch to map view`就可以得到如上所示的单体型地理分布图。通过单体型的地理分布直观展示，我们就可以对群体的扩散迁徙途径进行推断。

### 单体型序列的获得
像线粒体、叶绿体、Y染色体等这些在遗传过程中不发生重组的序列，我们直接把检测到的变异替换到参考序列中，就可以用于单体型网络构建。但是通常研究的二倍体基因组数据，由于存在重组，因此不能简单的拿一段序列就进行单体型网络构建，这样的单体型网络就失去了它本身的意义。如要用基因组上面的序列，可以通过以下步骤：  
1. 找到基因组上强连锁的区段；
2. 对这段区域中的SNP进行phase；

第一步可以用[PLINK](https://www.cog-genomics.org/plink2)，命令如下：
```shell
plink --noweb --bfile bed_prefix --blocks no-pheno-req
```
该命令运行会生成`.blocks.det`文件，该文件纪录了强连锁区段的起止位置、所包含的SNP。接下来，我们需要对感兴趣的强连锁区域的杂合SNP进行phase，构建单体型，该步可以使用软件BEAGLE。  
第二步，使用[BEAGLE](单体型网络（Haplotype Network）构建及展示)进行phase：
```
java -Xss5m -Xmx4g -jar beagle.jar gt=prefix.vcf out=phased.vcf chrom=[chr]:[start]-[end] 
```