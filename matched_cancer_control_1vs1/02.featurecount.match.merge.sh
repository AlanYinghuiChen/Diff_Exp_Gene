dataPath='/BioII/lulab_b/chenyinghui/project/exRNA/result_yunfan/CRC/tissue/01.featurecount'
outDir='/BioII/lulab_b/chenyinghui/project/exRNA/result_yunfan/CRC/tissue/02.featurecount_matched_merge'
sample_list='/BioII/lulab_b/chenyinghui/project/exRNA/result_yunfan/CRC/tissue/matched_sample.list'
sep="'\t'"

script="$outDir/script"
if [ ! -d $script ]; then mkdir -p $script ; fi

for p in `cat $sample_list`
do
i=`echo $p|sed 's/-/_/g'`
if [ -f "$script/${p}.merge.featurecount.R" ]; then rm "$script/${p}.merge.featurecount.R" ; fi
echo -e """
$i <- read.table('$dataPath/${p}-N.featurecounts.txt',header=F, skip=2, quote='',check.names=F, sep=$sep)
$i <- as.data.frame($i)

# Add gene ID column
geneID <- ${i}[,'V1']
rownames($i) <- geneID
all_counts_data <- ${i}['V1']
colnames(all_counts_data) <- c('geneID')

# Add control sample's gene read counts
${i}_data <- ${i}['V7']
colnames(${i}_data) <- c('${i}_N')
all_counts_data <- cbind(all_counts_data, ${i}_data)

# Add Cancer sample's gene read counts
$i <- read.table('$dataPath/${p}-T.featurecounts.txt',header=F, row.names=1, skip=2, quote='',check.names=F, sep=$sep)
$i <- as.data.frame($i)
${i}_data <- ${i}['V7']
colnames(${i}_data) <- c('${i}_T')
all_counts_data <- cbind(all_counts_data, ${i}_data)

write.table(all_counts_data, file='$outDir/${p}.featurecounts.txt', quote=F, sep=$sep, row.names=F, col.names=T)
""" > "$script/${p}.merge.featurecount.R"
done

for p in `cat $sample_list`
do
Rscript "$script/${p}.merge.featurecount.R"
done
