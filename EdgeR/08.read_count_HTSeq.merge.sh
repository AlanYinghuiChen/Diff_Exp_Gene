dataPath='/BioII/lulab_b/chenyinghui/project/exRNA/CRC/cfRNA/08.read_count_HTSeq'
outDir='/BioII/lulab_b/chenyinghui/project/exRNA/CRC/cfRNA/08.read_count_HTSeq.merge'
sample_list='/BioII/lulab_b/chenyinghui/project/exRNA/CRC/cfRNA/all_sample_list.txt'
sep="\\\t"

script="$outDir/script"
if [ ! -d $script ]; then mkdir -p $script ; fi

if [ -f "$script/merge.read_count.HTSeq.R" ]; then rm "$script/merge.read_count.HTSeq.R" ; fi

for n in `head -1 $sample_list`
do
i=`echo "$n" |sed 's/-/_/g'`
echo -e """
$i <- read.table('$dataPath/${n}.read_count.HTSeq.union.txt',header=F, quote='',check.names=F, sep='$sep')
$i <- as.data.frame($i)
geneID <- ${i}[,'V1']
rownames($i) <- geneID
all_counts_data <- ${i}['V1']
colnames(all_counts_data) <- c('geneID')
""" >> "$script/merge.read_count.HTSeq.R"
done

for n in `cat $sample_list`
do
i=`echo "$n" |sed 's/-/_/g'`
echo -e """

$i <- read.table('$dataPath/${n}.read_count.HTSeq.union.txt',header=F, row.names=1,  quote='',check.names=F, sep='$sep')
$i <- as.data.frame($i)
${i}_data <- ${i}['V2']
colnames(${i}_data) <- c('$i')

all_counts_data <- cbind(all_counts_data, ${i}_data)

""" >> "$script/merge.read_count.HTSeq.R"
done

echo -e """

write.table(all_counts_data, file='$outDir/all.read_count.HTSeq.union.txt', quote=F, sep='$sep', row.names=F, col.names=T)
""" >> "$script/merge.read_count.HTSeq.R"
