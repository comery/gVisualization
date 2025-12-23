perl ../bin/Draw_Dot.pl -GeneStr \
-Qinfo Query.len -Tinfo Target.len \
-Qcol "0,8,9" -Tcol "1,11,12" -Scol 16 \
-w 800 -h 800 -YtR \
-Tgff NPY4R.change.gff \
-GeneW 5.0 \
dotplot.svg \
blastn.tab.filter

