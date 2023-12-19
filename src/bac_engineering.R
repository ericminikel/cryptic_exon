options(stringsAsFactors=F)
setwd('~/d/sci/src/cryptic_exon')

write.fasta = function(textstring, path, title='', linelen = 50) {
  textvector = substring(textstring, seq(1, nchar(textstring), linelen), seq(linelen, nchar(textstring)+linelen, linelen))
  write(paste('>',title,sep=''), path)
  write(textvector, path, append=T)
}

# read in and break up reference sequences

human_fa = read.table('data/sequences/human_prnp_hg37_20_4663808_4682964.fa',skip=1)
human_sequence = paste(human_fa$V1, collapse='')
nchar(human_sequence)

human_fa_start = 4663808
human_fa_end = 4682964

human_exon1_start = 4666797 # transcription start site
human_exon1_end = 4667158

human_exon2_start = 4669784
human_exon2_end = 4669882

human_exon3_start = 4679857
human_exon3_end = 4682234 # transcription end site

human_cds_start = 4679867
human_cds_end = 4680628

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


# create frankenstein sequences

mouse_with_human_exon2 = paste(mouse_up_to_exon2, human_exon2, mouse_after_exon2, sep='')
substr(mouse_with_human_exon2, 6942, 7052)


human_up_to_exon2_minus6 = substr(human_up_to_exon2, 1, nchar(human_up_to_exon2)-6)
nchar(human_up_to_exon2)
nchar(human_up_to_exon2_minus6)
substr(human_up_to_exon2, nchar(human_up_to_exon2)-15, nchar(human_up_to_exon2))
substr(human_up_to_exon2, nchar(human_up_to_exon2_minus6)-9, nchar(human_up_to_exon2_minus6))
human_after_exon2_minus6 = substr(human_after_exon2, 7, nchar(human_after_exon2))
nchar(human_after_exon2)
nchar(human_after_exon2_minus6)
substr(human_after_exon2, 1, 15)
substr(human_after_exon2_minus6, 1, 9)
mouse_acceptor6 = substr(mouse_up_to_exon2, nchar(mouse_up_to_exon2)-6+1, nchar(mouse_up_to_exon2))
mouse_acceptor6
mouse_donor6 = substr(mouse_after_exon2, 1, 6)
mouse_donor6

human_with_mouse_splice6 = paste(human_up_to_exon2_minus6, mouse_acceptor6, human_exon2, mouse_donor6, human_after_exon2_minus6, sep='')
substr(human_with_mouse_splice6, human_exon2_start - human_fa_start + 1 - 10, human_exon2_end - human_fa_start + 1 + 10)



human_up_to_exon2_minus100 = substr(human_up_to_exon2, 1, nchar(human_up_to_exon2)-100)
nchar(human_up_to_exon2)
nchar(human_up_to_exon2_minus100)
substr(human_up_to_exon2, nchar(human_up_to_exon2)-15, nchar(human_up_to_exon2))
substr(human_up_to_exon2, nchar(human_up_to_exon2_minus100)-9, nchar(human_up_to_exon2_minus100))
human_after_exon2_minus100 = substr(human_after_exon2, 7, nchar(human_after_exon2))
nchar(human_after_exon2)
nchar(human_after_exon2_minus100)
substr(human_after_exon2, 1, 15)
substr(human_after_exon2_minus100, 1, 9)
mouse_acceptor100 = substr(mouse_up_to_exon2, nchar(mouse_up_to_exon2)-100+1, nchar(mouse_up_to_exon2))
mouse_acceptor100
mouse_donor100 = substr(mouse_after_exon2, 1, 100)
mouse_donor100

human_with_mouse_splice100 = paste(human_up_to_exon2_minus100, mouse_acceptor100, human_exon2, mouse_donor100, human_after_exon2_minus100, sep='')
substr(human_with_mouse_splice100, human_exon2_start - human_fa_start + 1 - 20, human_exon2_end - human_fa_start + 1 + 20)

# also create human cDNA with and without exon2

human_exon1_start = 4666797
human_exon1_end = 4667158 # this goes all the way up to the canonical GT donor
human_exon1 = substr(human_sequence, human_exon1_start - human_fa_start + 1, human_exon1_end - human_fa_start + 1)
nchar(human_exon1)

human_exon3_start = 4679857
human_exon3_end = 4682234

human_exon3 = substr(human_sequence, human_exon3_start - human_fa_start + 1, human_exon3_end - human_fa_start + 1)
human_exon3


human_cdna_exon1_3 = paste(human_exon1, human_exon3, sep='')
nchar(human_cdna_exon1_3)

human_cdna_exon1_2_3 = paste(human_exon1, human_exon2, human_exon3, sep='')
nchar(human_cdna_exon1_2_3)

nchar(human_cdna_exon1_2_3) - nchar(human_cdna_exon1_3)

# check
gregexpr(pattern='ATG',human_cdna_exon1_3)
gregexpr(pattern='ATG',human_cdna_exon1_2_3)


# Minigene engineering
# Start from TSS (assume we'll clone into vector w/ strong promoter)
# include 500bp on either side of each exon, wild-type human sequence

# how far from exon1 to exon2?
human_exon2_start - human_exon1_end # 2,626bp. so we can save significant length by cutting up this intron too

human_exon1_trail_500 = substr(human_sequence, human_exon1_start - human_fa_start + 1, human_exon1_end - human_fa_start + 1 + 500)
nchar(human_exon1_trail_500)

human_exon2_flank_500 = substr(human_sequence, human_exon2_start - human_fa_start + 1 - 500, human_exon2_end - human_fa_start + 1 + 500)
nchar(human_exon2_flank_500)

human_exon3_lead_500 = substr(human_sequence, human_exon3_start - human_fa_start + 1 - 500, human_exon3_end - human_fa_start + 1)
nchar(human_exon3_lead_500)

human_exon1_through_2_trail_500 = substr(human_sequence, human_exon1_start - human_fa_start + 1, human_exon2_end - human_fa_start + 1 + 500)
nchar(human_exon1_through_2_trail_500)

human_minigene_1 = paste(human_exon1_trail_500, human_exon2_flank_500, human_exon3_lead_500, sep='')
nchar(human_minigene_1)

human_minigene_2 = paste(human_exon1_through_2_trail_500, human_exon3_lead_500, sep='')
nchar(human_minigene_2)

# write out all files
write.fasta(human_sequence, 'outputs/human_sequence_check.fa', title='human_sequence_check')
write.fasta(mouse_sequence, 'outputs/mouse_sequence_check.fa', title='mouse_sequence_check')
write.fasta(mouse_with_human_exon2, 'outputs/mouse_with_human_exon2.fa', title='mouse_with_human_exon2')
write.fasta(human_with_mouse_splice6, 'outputs/human_with_mouse_splice6.fa', title='human_with_mouse_splice6')
write.fasta(human_with_mouse_splice100, 'outputs/human_with_mouse_splice100.fa', title='human_with_mouse_splice100')
write.fasta(human_cdna_exon1_3, 'outputs/human_cdna_exon1_3.fa', title='human_cdna_exon1_3')
write.fasta(human_cdna_exon1_2_3, 'outputs/human_cdna_exon1_2_3.fa', title='human_cdna_exon1_2_3')

write.fasta(human_minigene_1, 'outputs/human_minigene_1.fa', title='human_minigene_1')
write.fasta(human_minigene_2, 'outputs/human_minigene_2.fa', title='human_minigene_2')
write.fasta(human_cds, 'outputs/human_cds.fa', title='human_cds')

# make plots

ymin = -2.7
ymax = -0.5
xmin = human_fa_start
xmax = human_fa_end
landmarks = c(human_exon1_start, human_exon1_end, human_exon2_start, human_exon2_end, human_exon3_start, human_exon3_end)
landmark_labels = formatC(landmarks, format='fg', big.mark=',')
transcript = data.frame(starts = c(human_exon1_start),
                        ends = c(human_exon3_end))
exons = data.frame(starts = c(human_exon1_start, human_exon2_start, human_exon3_start),
                   ends = c(human_exon1_end, human_exon2_end, human_exon3_end))
cds = data.frame(starts = c(human_cds_start),
                 ends = c(human_cds_end))
reference_length = human_exon3_end - human_exon1_start
reference_length_kb = paste(formatC(reference_length/1000, format='f', digits=1), 'kb')
minigene_1_length = nchar(human_minigene_1)
minigene_1_length_kb = paste(formatC(minigene_1_length/1000, format='f', digits=1), 'kb')
minigene_2_length = nchar(human_minigene_2)
minigene_2_length_kb = paste(formatC(minigene_2_length/1000, format='f', digits=1), 'kb')

reference_y = -1.0
minigene_1_y = -1.5
minigene_2_y = -2.0

intron_lwd = 1
utr_lwd = 10
cds_lwd = 20


pdf('vector_graphics/minigenes_vs_reference_original.pdf',width=8,height=3)
par(mar=c(6,8,1,4))
plot(NA, NA, xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)
abline(h=ymin, lwd=3)
#abline(v=c(xmin, xmax), lwd=2)
arrows(x0=xmin+100, x1=xmin+1100, y0=ymin + 0.1, y1=ymin + 0.1, length=0.05, code=3, angle=90)
text(x=xmin+550, y=ymin+0.1, pos=3, labels='1 kb')
axis(side=1, at=landmarks, labels=c(landmark_labels[1:2],'','',landmark_labels[5:6]), lwd=0, lwd.ticks=1, cex.axis=0.8, las=2, srt=180)
axis(side=1, at=landmarks[3]-150, labels=landmark_labels[3], lwd=0, lwd.ticks=0, cex.axis=0.8, las=2, srt=180)
axis(side=1, at=landmarks[4]+150, labels=landmark_labels[4], lwd=0, lwd.ticks=0, cex.axis=0.8, las=2, srt=180)
mtext(side=1, line=1, text='hg37 coordinates')
# reference sequence
segments(y0=reference_y, x0=transcript$starts, x1=transcript$ends, lwd=intron_lwd, lend=1)
segments(y0=reference_y, x0=exons$starts, x1=exons$ends, lwd=utr_lwd, lend=1)
segments(y0=reference_y, x0=cds$starts, x1=cds$ends, lwd=cds_lwd, lend=1)
mtext(side=2, at=reference_y, line=0, text='reference', las=2)
mtext(side=4, at=reference_y, line=0.5, text=reference_length_kb, las=2)
# minigene 1 - 500 flanking each exon
segments(y0=minigene_1_y, x0=exons$starts + c(0,-500,-500), x1=exons$ends+c(500,500,0), lwd=intron_lwd, lend=1)
segments(y0=minigene_1_y, x0=exons$starts, x1=exons$ends, lwd=utr_lwd, lend=1)
segments(y0=minigene_1_y, x0=cds$starts, x1=cds$ends, lwd=cds_lwd, lend=1)
mtext(side=2, at=minigene_1_y, line=0, text='minigene 1', las=2)
mtext(side=4, at=minigene_1_y, line=0.5, text=minigene_1_length_kb, las=2)
# minigene 2 - retain intron 1, 500 flanking into intron 2
segments(y0=minigene_2_y, x0=exons$starts[1], x1=exons$ends[2] + 500, lwd=intron_lwd, lend=1)
segments(y0=minigene_2_y, x0=exons$starts[3] - 500, x1=exons$ends[3], lwd=intron_lwd, lend=1)
segments(y0=minigene_2_y, x0=exons$starts, x1=exons$ends, lwd=utr_lwd, lend=1)
segments(y0=minigene_2_y, x0=cds$starts, x1=cds$ends, lwd=cds_lwd, lend=1)
mtext(side=2, at=minigene_2_y, line=0, text='minigene 2', las=2)
mtext(side=4, at=minigene_2_y, line=0.5, text=minigene_2_length_kb, las=2)


dev.off()


# start for canonical 57-bp exon 1
transcript$starts[1] = 4667102
exons$starts[1] = 4667102
landmarks = c(4667102, human_exon1_end, human_exon2_start, human_exon2_end, human_exon3_start, human_exon3_end)
landmark_labels = formatC(landmarks, format='fg', big.mark=',')

pdf('vector_graphics/minigenes_vs_reference_canonical_exon1.pdf',width=8,height=3)
par(mar=c(6,8,1,4))
plot(NA, NA, xlim=c(xmin, xmax), ylim=c(ymin, ymax), xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)
abline(h=ymin, lwd=3)
#abline(v=c(xmin, xmax), lwd=2)
arrows(x0=xmin+100, x1=xmin+1100, y0=ymin + 0.1, y1=ymin + 0.1, length=0.05, code=3, angle=90)
text(x=xmin+550, y=ymin+0.1, pos=3, labels='1 kb')
axis(side=1, at=landmarks, labels=c(landmark_labels[1:2],'','',landmark_labels[5:6]), lwd=0, lwd.ticks=1, cex.axis=0.8, las=2, srt=180)
axis(side=1, at=landmarks[3]-150, labels=landmark_labels[3], lwd=0, lwd.ticks=0, cex.axis=0.8, las=2, srt=180)
axis(side=1, at=landmarks[4]+150, labels=landmark_labels[4], lwd=0, lwd.ticks=0, cex.axis=0.8, las=2, srt=180)
mtext(side=1, line=1, text='hg37 coordinates')
# reference sequence
segments(y0=reference_y, x0=transcript$starts, x1=transcript$ends, lwd=intron_lwd, lend=1)
segments(y0=reference_y, x0=exons$starts, x1=exons$ends, lwd=utr_lwd, lend=1)
segments(y0=reference_y, x0=cds$starts, x1=cds$ends, lwd=cds_lwd, lend=1)
mtext(side=2, at=reference_y, line=0, text='reference', las=2)
mtext(side=4, at=reference_y, line=0.5, text=reference_length_kb, las=2)
# minigene 1 - 500 flanking each exon
segments(y0=minigene_1_y, x0=exons$starts + c(0,-500,-500), x1=exons$ends+c(500,500,0), lwd=intron_lwd, lend=1)
segments(y0=minigene_1_y, x0=exons$starts, x1=exons$ends, lwd=utr_lwd, lend=1)
segments(y0=minigene_1_y, x0=cds$starts, x1=cds$ends, lwd=cds_lwd, lend=1)
mtext(side=2, at=minigene_1_y, line=0, text='minigene 1', las=2)
mtext(side=4, at=minigene_1_y, line=0.5, text=minigene_1_length_kb, las=2)
# minigene 2 - retain intron 1, 500 flanking into intron 2
segments(y0=minigene_2_y, x0=exons$starts[1], x1=exons$ends[2] + 500, lwd=intron_lwd, lend=1)
segments(y0=minigene_2_y, x0=exons$starts[3] - 500, x1=exons$ends[3], lwd=intron_lwd, lend=1)
segments(y0=minigene_2_y, x0=exons$starts, x1=exons$ends, lwd=utr_lwd, lend=1)
segments(y0=minigene_2_y, x0=cds$starts, x1=cds$ends, lwd=cds_lwd, lend=1)
mtext(side=2, at=minigene_2_y, line=0, text='minigene 2', las=2)
mtext(side=4, at=minigene_2_y, line=0.5, text=minigene_2_length_kb, las=2)


dev.off()
