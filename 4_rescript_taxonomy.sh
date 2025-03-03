qiime rescript get-silva-data \
  --p-version '138.1' \
  --p-target 'SSURef_NR99' \
  --o-silva-sequences silva-138.1-ssu-nr99-rna-seqs.qza \
  --o-silva-taxonomy silva-138.1-ssu-nr99-tax.qza \
  --p-ranks domain kingdom phylum class order suborder family subfamily genus

qiime rescript reverse-transcribe \
  --i-rna-sequences silva-138.1-ssu-nr99-rna-seqs.qza \
  --o-dna-sequences silva-138.1-ssu-nr99-seqs.qza

qiime rescript cull-seqs  \
  --i-sequences silva-138.1-ssu-nr99-seqs.qza \
  --o-clean-sequences silva-138.1-ssu-nr99-seqs-cleaned.qza

qiime rescript filter-seqs-length-by-taxon \
  --i-sequences silva-138.1-ssu-nr99-seqs-cleaned.qza \
  --i-taxonomy silva-138.1-ssu-nr99-tax.qza  \
  --p-labels Archaea Bacteria Eukaryota  \
  --p-min-lens 900 1200 1400  \
  --o-filtered-seqs silva-138.1-ssu-nr99-seqs-filt.qza  \
  --o-discarded-seqs silva-138.1-ssu-nr99-seqs-discard.qza

qiime rescript dereplicate  \
  --i-sequences silva-138.1-ssu-nr99-seqs-filt.qza \
  --i-taxa silva-138.1-ssu-nr99-tax.qza  \
  --p-mode 'uniq' \
  --o-dereplicated-sequences silva-138.1-ssu-nr99-seqs-derep-uniq.qza  \
  --o-dereplicated-taxa silva-138.1-ssu-nr99-tax-derep-uniq.qza

qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads  silva-138.1-ssu-nr99-seqs-derep-uniq.qza  \
  --i-reference-taxonomy silva-138.1-ssu-nr99-tax-derep-uniq.qza \
  --o-classifier silva-138.1-ssu-nr99-classifier.qza

qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path 18S_ASVs_lulu.fa \
  --output-path 18S_ASVs_lulu.qza

qiime feature-classifier classify-sklearn \
  --i-classifier silva-138.1-ssu-nr99-classifier.qza \
  --i-reads 18S_ASVs_lulu.qza \
  --o-classification taxonomy-fullssu.qza

qiime metadata tabulate \
  --m-input-file taxonomy-fullssu.qza \
  --o-visualization taxonomy-fullssu.qzv

qiime tools export \
  --input-path taxonomy-fullssu.qza \
  --output-path taxonomy-fullssu
