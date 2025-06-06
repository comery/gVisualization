## **THIS DIRECTORY HAS BEEN MOVED TO => comery/Data-Visulization**  
this copy will be forzen!

![example](auto_depth/example/circos.png)
example: Fruhstorferiola huayinensis mitochondrion, complete genome, download from NCBI. NC_031379

#### Summary
circos is very powerful software to illustrate genetic information, especially for mitochondrion and chloroplast genomes. Here I provide a useful perl script to convert genome annotation with genebank format to all features as circos input, and then draw a circos image.

#### Requirement

Overall, circos must be installed in your computer, though it is a little difficult to install.
and some perl modules are required too.

> Cwd, FindBin, Statistics::Descriptive, Bio::Perl


 Actually, there are two versions, one is “standard”, and the other is “auto_depth”. Literally, auto_depth means automatically calculate depth by mapping, so you need give fastq file and associaged softwares.

- for “standard” mode, you just need a genebank input file, and general configure file, like below.

```perl
# where circos you install
circos_path		=	/usr/local/bin/circos
#color
cds				=   102,194,165
rRNA			=	252,141,98
tRNA			=	123,50,148
D-loop        =   67,162,202
# whether draw GC content circle
gc				=	yes
win				=	50
gc_fill			=	146,197,222
# whether draw depth abundance circle
depth			=	no
depth_fill		= 	5,113,176
depth_file		=	../data/example.depth.txt
# whether draw base around circle
base 			=	no
# locus name's clolor showed on center of circle
locus_color		=	black
# gene name label color
label_color		=	black
# image
outdir			=	./example
png				=	yes
svg				=	yes
# color or file
background		=	white
#background		=	./background.png
```
let's get ball running.
```
perl draw_circos_for_mitogenome_auto_depth.pl -gb ../example.gb -conf mitogenome.auto_depth.conf 
```

- for “auto\_depth” mode, you also need a genebank input file and configure file, but what you should notice carefully is path setted in head of configure file - "circos_path", "bwa", and "samtools", as a result, paired fastq file should give in before file.

```
# where circos you install
circos_path		=	circos
#color
cds				=	102,194,165
rRNA			=	252,141,98
tRNA			=	123,50,148
# whether draw GC content circle
gc				=	yes
win				=	50
gc_fill			=	146,197,222
# whether draw depth abundance circle
bwa             =   bwa
samtools        =   samtools
depth_file      =
run_map			=	yes
depth_fill		= 	5,113,176
fq				=	data/test1.fq.gz,data/test2.fq.gz
# threads for baw
threads         =   2
# optional arguments for samtools
opts_samtools   =  -a -a

# whether draw base around circle
base 			=	no
# locus name's clolor showed on center of circle
locus_color		=	black
# gene name label color
label_color		=	black
# image
outdir			=	./example
png				=	yes
svg				=	yes
# color or file
background		=	white
#background		=	./background.png
```
now, you can run it by:

```shell
perl draw_circos_for_mitogenome_auto_depth.pl -gb data/SRR611174_mitoscaf.fa.gbf -conf mitogenome.auto_depth.conf

```
#### Note

##### The logic to use depth information

If you set "depth_file = xxx" and "run_map = yes" in the same time, it will directly use the existing depth file; otherwise, if you set the "run_map = yes", it will retrieve the fastq file path from "fq = xxx_1.fq, xxx_2.fq", and then run bwa mapping to generate depth file.


#### See more
recommend some good learning resources:

- http://qiubio.com/circos教程
- http://circos.ca/
