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
      mkdir $VERKKO_DIR/$CUR_DIR
      for depth in 30;
      do
        mkdir $VERKKO_DIR/$CUR_DIR/$depth
      done
    done
  done
done