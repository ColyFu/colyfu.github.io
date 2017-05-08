---
title: 判断fastq文件质量编码格式及Phred64转Phred33方法
date: 2017-04-25 21:37:56
tags:
- fastq
- phred33
- phred64
- 格式
categories:
- QC
---

前几天有朋友从网上下载了一批fastq文件，为Phred64格式，分析之前没发现，比对的时候BWA报错了。很多人可能没有接触过老Illumina数据，不知道怎么判断编码格式，也不知道怎么转为Phred33格式，所以简单说明一下怎么判断fastq文件的质量编码方式以及怎么将Phred64编码格式转为Phred33编码格式。

<!--more-->

* ### Phred64/33质量编码格式区别  
fastq文件中，碱基质量是用ASCII字符表示。Phred64编码格式，碱基质量值为字符的十进制ASCII码减去64。同理，Phred33编码格式，碱基质量值为字符的十进制ASCII码减去33。一般碱基的质量值范围为[0, 41]，X ten之后，质量值只有(2,7,11,22,27,32,37,42)几个梯度。因此，GATK的BQSR过程对于X ten的数据可能更加重要。  
Phred质量值Q和出错的概率P的关系为Q = -10*lg(P)，如碱基质量值为30，表示出错的概率为0.001，碱基质量值为20，表示出错的概率为0.01。
目前主流软件如BWA, GATK等都识别的是Phred33质量编码格式，如果为Phred64格式，则可能会报错，即使不报错，后续的分析也会有问题。因此，下载的数据如果不清楚编码格式，需要先判断，如果为Phred64，则需要转为Phred33格式。

* ### 判断fastq文件质量值编码格式  
格式判断既可以通过肉眼快速判断，也可以使用[下面的脚本](http://onzjn6hm6.bkt.clouddn.com/CheckFqQualityCode.pl)进行判断。如果文件少，没有编程基础，可以用快速判断的方法。利用程序判断更加准确、快速、可批量处理，因此有编程基础的人，应该尽量使用程序判断。

> **肉眼快速判断**：质量字符有数字`[0~9]`的为Phred33，有小写字母`[a~z]`的为Phred64。  
> **通过程序**：基本思路为将一定数量的reads质量值字符转为ASCII码，然后判断质量值的范围。     
> [脚本](http://onzjn6hm6.bkt.clouddn.com/CheckFqQualityCode.pl)的使用方法如下，第一个参数为fastq文件，自动判断是否为压缩文件，第二个参数为用于判断的reads数，默认为1000。
  ```bash
  perl CheckFqQualityCode.pl prefix.fq[.gz] [1000]
  ```

* ### Phred64格式转Phred33格式
Phred64格式转Phred33格式的原理很简单，只需在原有ASCII码的基础上减去64再加上33既可。大家可以自己写，也可以使用lh写的[seqtk](https://github.com/lh3/seqtk)工具，使用命令：
```bash
seqtk seq -VQ64 prefix.phred64.fq.gz | gzip > prefix.phred33.fq.gz
```