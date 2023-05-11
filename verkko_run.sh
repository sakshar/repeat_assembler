export DI=/Users/sakshar5068/Desktop/repeat_assembler
export REFS=$DI/data/ref/indel5
export READS=$DI/data/HiFi/in_del5
export VERKKO_DIR=$DI/verkko/in_del5


for repeat_size in 15000 20000 25000;
do
    for copy in 5 10;
    do
        for snp in 250 500 1000;
        do
            export CUR_DIR="$repeat_size"_"$copy"_"$snp"
            mkdir $VERKKO_DIR/$CUR_DIR
            for depth in 30;
            do
                mkdir $VERKKO_DIR/$CUR_DIR/$depth
                verkko -d $VERKKO_DIR/$CUR_DIR/$depth --hifi $READS/$CUR_DIR/$depth.fasta
                export FILE=$VERKKO_DIR/$CUR_DIR/$depth/assembly.fasta
                if [ -f $FILE ]; then
                    cp $FILE $VERKKO_DIR/$CUR_DIR/$depth/verkko.fasta
                else
                    echo ">0" >> $VERKKO_DIR/$CUR_DIR/$depth/verkko.fasta
                fi
                #$DROPBOX upload $HICANU_DIR/$CUR_DIR/$depth/canu.contigs.fasta $UPLOAD/$CUR_DIR/$depth/canu.fasta
                #$NUCMER -prefix=out -l 10 $REFS/$CUR_DIR.fasta $HICANU_DIR/$CUR_DIR/$depth/canu.contigs.fasta
                #$MUMMER --postscript -l -p plot out.delta -R $REFS/$CUR_DIR.fasta -Q $HICANU_DIR/$CUR_DIR/$depth/canu.contigs.fasta
                #mv plot.ps $HICANU_DIR/$CUR_DIR/$depth/ref.vs.canu.ps
                #rm plot.*
                #rm out.delta
                #$DROPBOX upload $HICANU_DIR/$CUR_DIR/$depth/ref.vs.canu.ps $UPLOAD/$CUR_DIR/$depth/ref.vs.canu.ps
            done
        done
    done
done