qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path Anti-IL22-baseline-batch1-sequences \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-paired-end.qza

qiime demux summarize \
  --i-data demux-paired-end.qza \
  --o-visualization demux.qzv

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux-paired-end.qza \
  --p-trim-left-f 20 \
  --p-trim-left-r 17 \
  --p-trunc-len-f 300 \
  --p-trunc-len-r 285 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv

qiime tools export \
--input-path table.qza \
--output-path exported

biom convert \
--input-fp exported/feature-table.biom \
--output-fp feature-table.tsv \
--to-tsv
