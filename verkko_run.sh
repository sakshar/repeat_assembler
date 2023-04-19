export DI=/Users/sakshar5068/Desktop/repeat_assembler
export REFS=$DI/data/ref
export READS=$DI/data/HiFi
export VERKKO_DIR=$DI/verkko


for repeat_size in 25000;
do
    for copy in 5 10;
    do
        for snp in 250 500 1000;
        do
            export CUR_DIR="$repeat_size"_"$copy"_"$snp"
            for depth in 30;
            do
                verkko -d $VERKKO_DIR/$CUR_DIR/$depth --hifi $READS/$CUR_DIR/$depth.fasta
                cp $VERKKO_DIR/$CUR_DIR/$depth/assembly.fasta $VERKKO_DIR/$CUR_DIR/$depth/verkko.fasta
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