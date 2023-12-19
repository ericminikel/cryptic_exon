options(stringsAsFactors=F)
library(tidyverse)
library(janitor)
setwd('~/d/sci/src/cryptic_exon')

source('../helper.R')



#####
# Figure 1A
####

txscripts = tribble(
  ~y, ~enst, ~start, ~stop, ~canonical,
  4, 'ENST00000379440.9*', 	4686456,	4686512, T,
  4, 'ENST00000379440.9*', 	4699211,	4701588, T,
  2, 'ENST00000430350.2',   4686350,	4686508, F,
  2, 'ENST00000430350.2', 	4699211,	4701588, F,
  1, 'ENST00000457586.2', 	4686543,	4686736, F,
  1, 'ENST00000457586.2',   4699211,	4701590, F,
  3, 'ENST00000424424.2',   4686459,	4686512, F,
  3, 'ENST00000424424.2', 	4699216,	4701590, F,
)

cds = tribble(
  ~start, ~stop,
  4699221,	4699982
)

# grch38_coding_offset = 4699605 - 4680251

# per https://en.wikipedia.org/wiki/BED_(file_format)
# note that start is inclusive, stop is non-inclusive so it is [start, stop)

# GRCh38 coords:

human_exon1_disp_start = 4686150
human_exon1_disp_stop = 4686736 # 4686512

human_exon1_atgs = human_exon1_disp_start + c(8, 24, 44, 70)

human_exon2_start = 4689138
human_exon2_stop  = 4689236

human_exon2_atgs = human_exon2_start + 39

human_exon3_start = 4699211
human_exon3_stop = 4701588

human_exon3_atgs = human_exon3_start + 10

all_atgs = c(human_exon1_atgs, human_exon2_atgs, human_exon3_atgs)

exons = tribble(
  ~exon, ~start, ~stop,
  1, human_exon1_disp_start, human_exon1_disp_stop,
  2, human_exon2_start, human_exon2_stop,
  3, human_exon3_start, human_exon3_stop
)

overhang = 10
bases_e1 = tibble(base = seq(human_exon1_disp_start-overhang, human_exon1_disp_stop+overhang-1, 1),
                  exon = 1)
bases_e2 = tibble(base = seq(human_exon2_start-overhang, human_exon2_stop+overhang-1, 1),
                  exon = 2)
bases_e3 = tibble(base = seq(human_exon3_start-overhang, human_exon3_stop+overhang-1, 1),
                  exon = 3)
bases = rbind(bases_e1, bases_e2, bases_e3)

dp = tibble(exon=integer(0), base=integer(0), depth=integer(0), tissue=character(0))

files = list.files('data/gtex_depth/', full.name = T)
for (file in files) {
  d = read.table(file, 
                 skip=1, 
                 col.names = c('chr','start','stop','depth'),
                 comment.char='#')
  d$tissue = gsub('\\.txt','',gsub('.*\\/','',file))
  
  crossing(bases, d) %>%
    filter(base >= start & base < stop) %>%
    select(exon, base, depth, tissue) -> this_dp
  
  dp = rbind(dp, this_dp)
}

dp %>%
  mutate(norm_to = case_when(exon %in% 1:2 ~ 1,
                             exon %in% 3 ~ 3)) %>%
  group_by(norm_to, tissue) %>%
  mutate(reldepth = depth / max(depth)) %>%
  ungroup() %>%
  group_by(exon, base) %>%
  summarize(.groups='keep',
            rel_mean = mean(reldepth),
            rel_l95 = lower(reldepth),
            rel_u95 = upper(reldepth),
            rel_min = min(reldepth),
            rel_max = max(reldepth)) %>%
  ungroup() -> dp_by_base




resx=300
png('display_items/figure-1a.png',width=13*resx,height=4*resx,res=resx)
depth_color = '#FA9A50'
start_color = '#00B0F0'

exons %>%
  mutate(len = stop - start) %>%
  mutate(width = len / sum(len)) -> widths
layout_matrix = matrix(1:8, byrow=T, nrow=2)
layout(layout_matrix, widths=c(.05, widths$width), heights=c(1,.25))
par(mar=c(6,0,3,0))
plot(NA, NA, xlim=0:1, ylim=0:1, axes=F, ann=F)
axis(side=4, line=0, at=0:4/4, labels=NA, las=2, tck=0.025)
axis(side=4, line=-3.5, at=0:4/4, labels=percent(0:4/4), las=2, lwd=0)
mtext(side=4, line=-4, text='depth (% max)')
par(mar=c(6,1,3,1))
for (i in 1:nrow(exons)) {
  xlims = c(exons$start[i], exons$stop[i]) + overhang*c(-1,1)
  ylims = c(0, 1)
  
  xbigs = seq(round(min(xlims), -3), round(max(xlims), -3), 1000)
  xats = seq(round(min(xlims), -2), round(max(xlims), -2), 100) 
  xmicros = seq(round(min(xlims), -1), round(max(xlims), -1), 10) 
  
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  axis(side=1, at=xmicros, tck=-0.0125, labels=NA)
  axis(side=1, at=xats, tck=-0.025, labels=NA)
  #axis(side=1, at=xats, labels=formatC(xats, big.mark=',', format='f', digits=0), tck=-0.025)
  axis(side=1, at=xbigs, labels=formatC(xbigs, big.mark=',', format='f', digits=0), tck=-0.05, las=2)
  axis(side=1, at=xlims, labels=formatC(xlims,  big.mark=',', format='f', digits=0), tck=-0.065, las=2)
  points(dp_by_base$base, dp_by_base$rel_mean, lwd=2, type='l', col=depth_color)
  polygon(x=c(dp_by_base$base, rev(dp_by_base$base)), y=c(dp_by_base$rel_min, rev(dp_by_base$rel_max)), col=alpha(depth_color,.15), border=NA)
  
  mtext(side=3, line=0, text=paste0('exon ',i), cex=0.8)
  
  points(x=all_atgs, y=rep(.9, length(all_atgs)), pch=25, bg=start_color, lwd=0)
}

ylims = range(txscripts$y) + c(-0.5, 0.5)
par(mar=c(0,1,0,1))
plot(NA, NA, xlim=0:1, ylim=ylims, axes=F, ann=F)
par(xpd=T)
txscripts %>%
  distinct(y, enst, canonical) -> ensts
mtext(side=4, line = -2, adj=0, at=ensts$y, text=ensts$enst, cex=0.6, las=2, font=ifelse(ensts$canonical,2,1))
par(xpd=F)
utrheight = 0.2 
cdsheight=0.4
ensts %>%
  crossing(cds) -> cds_each
for (i in 1:nrow(exons)) {
  xlims = c(exons$start[i], exons$stop[i]) + overhang*c(-1,1)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  rect(xleft=txscripts$start, xright=txscripts$stop, ybottom=txscripts$y-utrheight/2, ytop=txscripts$y+utrheight/2, border=NA, col=ifelse(txscripts$canonical, '#000000', '#999999'))
  rect(xleft=cds_each$start, xright=cds_each$stop, ybottom=cds_each$y-cdsheight/2, ytop=cds_each$y+cdsheight/2, border=NA, col=ifelse(cds_each$canonical, '#000000', '#999999'))
  
}


dp_by_base %>%
  inner_join(txscripts, by=c('base'='start')) %>%
  filter(exon==1)

dp_by_base %>%
  filter(base %in% 4686456:4686512) %>%
  summarize(mean(rel_mean)) %>%
  pull() -> mean_depth_canonical_exon1

dp_by_base %>%
  filter(base %in% 4686350:4686455) %>%
  summarize(mean(rel_mean)) %>%
  pull() -> mean_depth_extended_portion_exon1

dp_by_base %>%
  filter(base %in% 4686150:4686349) %>%
  summarize(mean(rel_mean)) %>%
  pull() -> mean_depth_claimed_portion_exon1

dp_by_base %>%
  filter(base %in% 4686150:max(human_exon1_atgs)) %>%
  summarize(mean(rel_mean)) %>%
  pull() -> mean_depth_uorfs_exon1



dev.off()









####
# Figure 1B
####


resx=300
png('display_items/figure-1b.png',width=6.5*resx,height=2*resx,res=resx)

# read in and break up reference sequences

ham_fa = read.table('data/sequences/mesAur1_PRNP_environs.fa',skip=1)
ham_raw = paste(ham_fa$V1, collapse='')
nchar(ham_raw)

ham_tss = 'GGCGTCCGAGC' # 'ggcgtccgag'
ham_end = 'ATGATGAGCACCC'

tss_location = gregexpr(ham_tss, ham_raw)[[1]][1]
end_location = gregexpr(ham_end, ham_raw)[[1]][1]

ham_sequence = substr(ham_raw, tss_location, end_location + nchar(ham_end) - 1)
nchar(ham_sequence)
# substr(ham_sequence, nchar(ham_sequence)-10, nchar(ham_sequence))

ham_fa_start = 1
ham_exon1_start = gregexpr(ham_tss, ham_sequence)[[1]][1]
ham_exon1_end = gregexpr('CGCGTCGGTG', ham_sequence)[[1]][1] + nchar('CGCGTCGGTG')

ham_exon2_start = gregexpr('GAAGGACTCCT', ham_sequence)[[1]][1]
ham_exon2_end = gregexpr('TCACAGCAGATC', ham_sequence)[[1]][1] + nchar('TCACAGCAGATC')

ham_exon3_start = gregexpr('ATCAGCCATC', ham_sequence)[[1]][1]
ham_exon3_end = gregexpr(ham_end, ham_sequence)[[1]][1] + nchar(ham_end)

ham_cds_start = gregexpr('ATGGCGAAC', ham_sequence)[[1]][1]
ham_cds_end = gregexpr('GTGGGATGA', ham_sequence)[[1]][1] + nchar('GTGGGATGA')

human_fa = read.table('data/sequences/human_prnp_hg37_20_4663808_4682964.fa',skip=1)
human_sequence = paste(human_fa$V1, collapse='')
nchar(human_sequence)

human_fa_start = 4686456
human_fa_end = 4701590

human_exon1_start = 4686456 # transcription start site
human_exon1_end = 4686512

human_exon2_start = 4689138
human_exon2_end = 4689236

human_exon3_start = 4699211
human_exon3_end = 4701588 # transcription end site

human_cds_start = 4699221
human_cds_end = 4699982

human_exon1 = substr(human_sequence, human_exon1_start - human_fa_start + 1, human_exon1_end - human_fa_start + 1)
nchar(human_exon1)

human_exon2 = substr(human_sequence, human_exon2_start - human_fa_start + 1, human_exon2_end - human_fa_start + 1)
nchar(human_exon2)

human_exon3 = substr(human_sequence, human_exon3_start - human_fa_start + 1, human_exon3_end - human_fa_start + 1)
nchar(human_exon3)

human_cds = substr(human_sequence, human_cds_start - human_fa_start + 1, human_cds_end - human_fa_start + 1)
nchar(human_cds)


# BAC engineering

human_up_to_exon2 = substr(human_sequence, 1, human_exon2_start - human_fa_start)
nchar(human_up_to_exon2)
substr(human_up_to_exon2, nchar(human_up_to_exon2)-6, nchar(human_up_to_exon2))

human_after_exon2 = substr(human_sequence, human_exon2_end - human_fa_start + 2, nchar(human_sequence))
substr(human_after_exon2, 1, 6)

human_sequence_check = paste(human_up_to_exon2, human_exon2, human_after_exon2, sep='')
nchar(human_sequence_check)
human_sequence_check == human_sequence

mouse_fa = read.table('data/sequences/mouse_prnp_mm10_2_131905248_131939011.fa',skip=1)
mouse_sequence = paste(mouse_fa$V1, collapse='')
nchar(mouse_sequence)

mouse_fa_start = 131905248
mouse_fa_end = 131939011
mouse_exon2_start = 131912195
mouse_exon2_end = 131912292

mouse_exon1_start = 131909928
mouse_exon1_end = 131910003

mouse_exon3_start = 131936420
mouse_exon3_end = 131938436
mouse_cds_start = 131936430
mouse_cds_end = 131937194

mouse_exon1 = substr(mouse_sequence, mouse_exon1_start - mouse_fa_start + 1, mouse_exon1_end - mouse_fa_start + 1)
nchar(mouse_exon1)

mouse_exon2 = substr(mouse_sequence, mouse_exon2_start - mouse_fa_start + 1, mouse_exon2_end - mouse_fa_start + 1)
nchar(mouse_exon2)

mouse_exon3 = substr(mouse_sequence, mouse_exon3_start - mouse_fa_start + 1, mouse_exon3_end - mouse_fa_start + 1)
nchar(mouse_exon3)

mouse_cds = substr(mouse_sequence, mouse_cds_start - mouse_fa_start + 1, mouse_cds_end - mouse_fa_start + 1)
nchar(mouse_cds)


mouse_up_to_exon2 = substr(mouse_sequence, 1, mouse_exon2_start - mouse_fa_start)
nchar(mouse_up_to_exon2)
substr(mouse_up_to_exon2, nchar(mouse_up_to_exon2)-6, nchar(mouse_up_to_exon2))

mouse_after_exon2 = substr(mouse_sequence, mouse_exon2_end - mouse_fa_start + 2, nchar(mouse_sequence))
substr(mouse_after_exon2, 1, 6)

mouse_sequence_check = paste(mouse_up_to_exon2, mouse_exon2, mouse_after_exon2, sep='')
nchar(mouse_sequence_check)
mouse_sequence_check == mouse_sequence



tss_starts = c(mouse_exon1_start, ham_exon1_start, human_exon1_start)
landmarks = tibble(species = c('mouse','hamster','human'),
                   exon1_start = c(mouse_exon1_start, ham_exon1_start, human_exon1_start) - tss_starts +1,
                   exon2_start = c(mouse_exon2_start, ham_exon2_start, human_exon2_start) - tss_starts +1,
                   exon3_start = c(mouse_exon3_start, ham_exon3_start, human_exon3_start) - tss_starts +1,
                   exon1_end   = c(mouse_exon1_end,   ham_exon1_end,   human_exon1_end)   - tss_starts +1,
                   exon2_end   = c(mouse_exon2_end,   ham_exon2_end,   human_exon2_end)   - tss_starts +1,
                   exon3_end   = c(mouse_exon3_end,   ham_exon3_end,   human_exon3_end)   - tss_starts +1,
                   cds_start   = c(mouse_cds_start,   ham_cds_start,   human_cds_start)   - tss_starts +1,
                   cds_end     = c(mouse_cds_end,     ham_cds_end,     human_cds_end)     - tss_starts +1,
                   exon2_fill = c('#000000', '#A3A3A3', '#FFFFFF'),
                   y = c(2,1,0))

intron_lwd = 1
utr_lwd = 10
cds_lwd = 20

default_fill = '#000000'
utr_height = .25
cds_height = .5

par(mar=c(0,4,0,1))
ylims = c(-0.5, 2.5)
xlims = c(0, max(landmarks$exon3_end))
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
mtext(side=2, at=landmarks$y, text=landmarks$species, las=2, line=.75)
segments(x0=landmarks$exon1_start, x1=landmarks$exon3_end, y0=landmarks$y, lwd=intron_lwd)
rect(xleft=landmarks$exon1_start, xright=landmarks$exon1_end, ybottom=landmarks$y-utr_height/2, ytop=landmarks$y+utr_height/2, col=default_fill, border=NA)
rect(xleft=landmarks$exon2_start, xright=landmarks$exon2_end, ybottom=landmarks$y-utr_height/2, ytop=landmarks$y+utr_height/2, col=landmarks$exon2_fill, border='#000000', lwd=.25)
rect(xleft=landmarks$exon3_start, xright=landmarks$exon3_end, ybottom=landmarks$y-utr_height/2, ytop=landmarks$y+utr_height/2, col=default_fill, border=NA)
rect(xleft=landmarks$cds_start, xright=landmarks$cds_end, ybottom=landmarks$y-cds_height/2, ytop=landmarks$y+cds_height/2, col=default_fill, border=NA)
scale_x = c(23000,28000)
scale_y = 1.25
arrows(x0=scale_x[1], x1=scale_x[2], y0=scale_y, code=3, angle=90, length=0.05)
text(x=mean(scale_x), y=scale_y, pos=3, labels='5 kb')

dev.off()




####
# Other analyses
####


all_kozak = read_tsv('data/annotations/noderer-2014-supplement-table-s2.txt',skip=1,col_types=cols()) %>%
  clean_names() %>%
  mutate(percentile_to_all = cume_dist(efficiency)) %>%
  mutate(relative_efficiency = efficiency / max(efficiency))

mane = read_tsv('data/mane/MANE_transcripts_v1.0.tsv.gz', col_types=cols())

mane %>%
  select(transcript_id, seq, five_prime_utr_length, three_prime_utr_length) %>%
  mutate(sequence = gsub('T','U',toupper(substr(seq, five_prime_utr_length - 5, five_prime_utr_length + 5)))) %>% # Generate the 11 nucleotide context for each transcript
  filter(nchar(sequence)==11) %>%  # Filter out <11 nucleotide contexts - some without no or short 5' UTRs
  select(-seq) %>%
  left_join(select(all_kozak, sequence, relative_efficiency), by=c('sequence')) %>%
  mutate(percentile_to_mane = cume_dist(relative_efficiency)) -> mane_annotated

orfs = read_tsv('data/annotations/prnp_orfs.tsv', col_types=cols()) %>%
  mutate(sequence = gsub('T','U',sequence)) %>%
  left_join(select(all_kozak, sequence, percentile_to_all, relative_efficiency), by='sequence') %>%
  left_join(select(mane_annotated, sequence, percentile_to_mane), by=c('sequence')) %>%
  filter(description %in% c('novel uORF','canonical ORF')) %>%
  group_by(order) %>%
  slice(1)

write_tsv(orfs, 'outputs/prnp_orfs_strength.tsv')

leg = tribble(
  ~id, ~desc, ~color,
  'all','all possible\nKozak sequences\n','#FFA824BF',
  'mane','canonical start codons\nof human genes','#62B1F6FF',
)



resx=300
png('display_items/figure-1d.png',width=4.5*resx,height=3*resx,res=resx)
par(mar=c(3,3,1,1))
xlims = c(0,1)
# ylims = c(0,4)
ylims = c(0, 10000)
plot(NA, NA, xlim = xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=0:4/4, labels=NA, cex.axis=0.8)
axis(side=1, at=0:4/4, lwd=0, labels=percent(0:4/4), cex.axis=0.8, line=-0.5)
mtext(side=1, line=1.6, text='relative translational efficiency')
axis(side=2, at=0:10*1000, las=2, labels=NA, tck=-0.025)
axis(side=2, at=0:2*5000, las=2, labels=NA, tck=-0.05)
axis(side=2, at=0:2*5000, las=2, labels=paste0(c(0,5,10),'K'), lwd=0, line=-0.5)
mtext(side=2, line=2, text='number of Kozak sequences')
#axis(side=2, at=0:3, las=2, labels=0:3)
breaks = seq(0,1,.05)
par(new=T)
hist(all_kozak$relative_efficiency, breaks=breaks, xlim = xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', border = NA, col=leg$color[leg$id=='all'])
par(new=T)
hist(mane_annotated$relative_efficiency, breaks=breaks, xlim = xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', border = NA, col=leg$color[leg$id=='mane'])
abline(v=orfs$relative_efficiency, lty=3, col='black')
mtext(side=3, line=0, at=orfs$relative_efficiency, text=orfs$description)
legend('topleft', leg$desc, col=leg$color, pch=15, cex=0.7, bty='n')
dev.off()