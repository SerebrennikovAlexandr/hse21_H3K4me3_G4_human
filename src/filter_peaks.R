source('lib.R')

###

NAME <- 'H3K4me3_H1.ENCFF178RTX.hg19'
#NAME <- 'H3K4me3_H1.ENCFF236DPM.hg19'

###

bed_df <- read.delim(paste0(DATA_DIR, NAME, '.bed'), as.is = TRUE, header = FALSE)
colnames(bed_df) <- c('chrom', 'start', 'end', 'name', 'score')
bed_df$len <- bed_df$end - bed_df$start

# Look at the peak's length
bed_df %>%
  arrange(-len) %>%
  head(50)

# Remove long peaks in H3K4me3_H1.ENCFF236DPM.hg19
#bed_df <- bed_df %>%
#  arrange(-len) %>%
#  filter(len < 5500)

# Remove long peaks in H3K4me3_H1.ENCFF178RTX.hg19
bed_df <- bed_df %>%
  arrange(-len) %>%
  filter(len < 4500)
  
ggplot(bed_df) +
  aes(x = len) +
  geom_histogram() +
  ggtitle(NAME, subtitle = sprintf('Number of peaks = %s', nrow(bed_df))) +
  theme_bw()
ggsave(paste0('filter_peaks.', NAME, '.filtered.hist.pdf'), path = OUT_DIR)

bed_df %>%
  select(-len) %>%
  write.table(file=paste0(DATA_DIR, NAME ,'.filtered.bed'),
            col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)

