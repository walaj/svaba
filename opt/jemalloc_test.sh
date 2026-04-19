REGION=chr12:10,000,000-12,000,000
T=/Users/jeremiahwala/Desktop/svaba_compare/SG.wgs.UCLA.2025.01.tumor_cleaned.recal.bam
N=/Users/jeremiahwala/Desktop/svaba_compare/blood.recal.bam

# prime the page cache so cold-I/O doesn't dominate
REF=${HOME}/ref_genome/Homo_sapiens_assembly38.fasta
for f in $REF ${REF}.bwt ${REF}.sa ${REF}.pac
do cat "$f" > /dev/null; done

# baseline (system malloc) 3 trials
for t in 1 2 3; do
  /usr/bin/time -l ~/git/svaba/build/svaba run \
    -t $T -n $N  -G $REF \
    -k $REGION -p 16 -a base_$t 2>&1 \
  | tail -15 | grep -E "real|user|sys|maximum resident"
done

# jemalloc 3 trials
for t in 1 2 3; do
  /usr/bin/time -l env \
    DYLD_INSERT_LIBRARIES=$(brew --prefix jemalloc)/lib/libjemalloc.dylib \
    DYLD_FORCE_FLAT_NAMESPACE=1 \
    MALLOC_CONF=background_thread:true,narenas:24,dirty_decay_ms:10000 \
    ~/git/svaba/build/svaba run \
    -t $T -n $N -G $REF \
    -k $REGION -p 16 -a jem_$t 2>&1 \
  | tail -15 | grep -E "real|user|sys|maximum resident"
done
