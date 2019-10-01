#!/bin/bash


#Normalise reads and decontam befor a multiple assembly to lower the memory requirments 

##Logging
exec > >(tee -i logfile.txt)
exec 2>&1

##init global variables
TOPD=$(pwd)

aws s3 sync s3://matt-storey-sydney/sponge_reads . --exclude "*" --include "MH-s*" 


##Sort each PE pair into its own directory
for FILE in *fq.gz; do 
    DIR=$( echo $FILE | rev |cut -c9- | rev );
    if [ ! -d "$DIR" ]; then
        mkdir $DIR 
    fi;
    GZ=( ${DIR}*.fq.gz )
    if [ -f $GZ ]; then
        mv $GZ $DIR
    fi;
done

cd $TOPD

echo 'download the references for mapping'
aws s3 cp s3://matt-storey-sydney/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz
aws s3 cp s3://matt-storey-sydney/fusedERPBBmasked2.fa.gz fusedERPBBmasked2.fa.gz
mkdir human_contam
mv hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz human_contam
cd human_contam && bbmap.sh ref=hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz -Xmx50g
cd ..
mkdir ecoli_contam
mv fusedERPBBmasked2.fa.gz ecoli_contam
cd ecoli_contam && bbmap.sh ref=fusedERPBBmasked2.fa.gz -Xmx23g
cd $TOPD



cd $TOPD
for d in ${TOPD}/*; do
        if [ -d $d ]; then cd $d;
            for GZ in $(ls *fq.gz | rev | cut -c 8- | rev | uniq );
                do


 	         #echo "subsampling for testing"
		# reformat.sh in1=${GZ}1.fq.gz in2=${GZ}2.fq.gz out1=${GZ}sub.1.fq.gz out2=${GZ}sub.2.fq.gz reads=200000

		echo "Trim adapters for ${GZ}";
                bbduk.sh -Xmx300g in1=${GZ}1.fq.gz in2=${GZ}2.fq.gz out=trimmed.fq ktrim=r k=23 mink=11 \
                hdist=1 ref=/home/ubuntu/miniconda3/opt/bbmap-38.22-0/resources/adapters.fa tbo tpe maxns=0 trimq=10 qtrim=r maq=12;
                echo ""
                
                echo "Remove small contaminants from ${GZ}"
                bbduk.sh -Xmx300g in=trimmed.fq out=${GZ}.filtered.fq k=31 \
                ref=/home/ubuntu/miniconda3/opt/bbmap-38.22-0/resources/sequencing_artifacts.fa.gz,/home/ubuntu/miniconda3/opt/bbmap-38.22-0/resources/phix174_ill.ref.fa.gz
                echo 

                echo "Remove masked ecoli from ${GZ}"
                bbmap.sh -Xmx300g minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 \
                path=$TOPD/ecoli_contam/ qtrim=rl trimq=10 untrim in=${GZ}.filtered.fq outu=ecoli_cleaned.fq
                echo ""
                
                echo "Remove masked human from ${GZ}"
                bbmap.sh -Xmx300g minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 \
                path=$TOPD/human_contam/ qtrim=rl trimq=10 untrim in=ecoli_cleaned.fq outu=human_cleaned.fq
                echo ""
                
                echo "Error-correct ${GZ} 1st run - bbmerge.sh"
                bbmerge.sh -Xmx300g in=human_cleaned.fq out=ecco.fq ecco mix vstrict adapters=default
                echo ""
                
                echo "Error-correct ${GZ} 2nd run - clumpify.sh" 
                clumpify.sh -Xmx30g in=ecco.fq out=eccc.fq ecc passes=6
                echo ""
                
                echo "Error-correct ${GZ} 3rd run - tadpole.sh"
                tadpole.sh -Xmx300g in=eccc.fq out=ecct.fq ecc
                                
                echo "Merge"
                bbmerge.sh -Xmx300g in=eccc.fq out=${GZ}clean.merged.fq outu=${GZ}clean.unmerged.fq \
                rem k=62 extend2=50 adapters=default
                
                echo "Cleaning up "    
                rm trimmed.fq ecco.fq eccc.fq ecct.fq human_cleaned.fq ecoli_cleaned.fq 

                echo "Spadesing"
                spades.py -k 21,41,71,101,127 -o ${GZ}clean.merged.spades --12 ${GZ}clean.unmerged.fq --merged ${GZ}clean.merged.fq -t 96 -m 700


 	            done
        fi
done




cd $TOPD
for d in ${TOPD}/*; do
        if [ -d $d ]; then cd $d;
            for GZ in $(ls *fq.gz | rev | cut -c 8- | rev | uniq );
                do
			bbnorm.sh in=${GZ}clean.merged.fq out=${GZ}clean.merged.norm.fq target=175
			bbnorm.sh in=${GZ}clean.unmerged.fq out=${GZ}clean.norm.unmerged.fq target=175
                done
        fi
done

#run co-assembly here


