#!/bin/bash
#Input Format
#Family.ID       Child   Mother  Father  Flag    Variant.ID Gene.Name

IN=$1

awk '{print $2"\n"$3"\n"$4}' <(sed '1d' $IN) |uniq > $(dirname $IN)/$(basename ${IN%.*}).list

while read line; do  
    tmp=$(find \
        $(mysql --defaults-group-suffix=seqdb --defaults-file=/nfs/goldstein/software/annodb/.my.cnf -e "select AlignSeqFileLoc from seqdbClone where CHGVID = '$line' ;" |head -n2|tail -n1) \
        -name '*final.bam')

    if [ -z "$tmp" ]
        then >&2 echo "ERROR: Did not find BAM for $line"
    elif [ $(du -sh $tmp|wc -l) -gt 1 ]
        then >&2 echo "ERROR: Found more than one BAM for $line"
             >&2 echo "$tmp"
    else
        echo "$(readlink $tmp)"
    fi

done<$(dirname $IN)/$(basename ${IN%.*}).list > $(dirname $IN)/$(basename ${IN%.*}).bamloc

sed 's/\/nfs/\\\\10.73.52.21/' $(dirname $IN)/$(basename ${IN%.*}).bamloc |sed 's/10.73.52.21\/fastq/igm-avere.igm.cumc.columbia.edu\/fastq/' | sed 's/\//\\/g' > $(dirname $IN)/$(basename ${IN%.*}).bamwinloc

while read i
    do PRO_LOC=$(grep $(echo $i|awk '{print $2}') $(dirname $IN)/$(basename ${IN%.*}).bamwinloc|head -n1)
       DAD_LOC=$(grep $(echo $i|awk '{print $3}') $(dirname $IN)/$(basename ${IN%.*}).bamwinloc|head -n1)
       MOM_LOC=$(grep $(echo $i|awk '{print $4}') $(dirname $IN)/$(basename ${IN%.*}).bamwinloc|head -n1)
       VAR_CHR=$(echo $i|awk '{print $6}')
       VAR_LOC=$(echo $i|awk '{print $7}')
       VAR_GENE=$(echo $i|awk '{print $9}')
       PRO=$(echo $i|awk '{print $2}')
       echo "$PRO_LOC $DAD_LOC $MOM_LOC $VAR_CHR $VAR_LOC $((VAR_LOC-40)) $((VAR_LOC+40)) $PRO $VAR_GENE";done < <(sed 's/ /_/g' $IN|sed 's/-/\t/'|sed 's/-/\t/'|sed '1d'|sort -k6,6n) > $(dirname $IN)/$(basename ${IN%.*}).info.txt

       awk -v dir=$(pwd|sed 's/\/nfs/\\\\10.73.52.21/'|sed 's/\//\\\\/g') '{print "#"$8" "$4":"$5"\n""new\ngenome 1kg_v37\nload " $1"\nload " $2"\nload "$3"\nsnapshotDirectory \\"dir"\ngoto "$4":"$6"-"$7"\nsort position\nsnapshot "$8"."$9"."$4"-"$5".png\ncollapse\nsnapshot "$8"."$9"."$4"-"$5".collapsed.png\n"}' $(dirname $IN)/$(basename ${IN%.*}).info.txt >$(dirname $IN)/$(basename ${IN%.*}).batch
