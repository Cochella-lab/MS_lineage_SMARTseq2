################
#### Prepare tables
################

#### Tss filtered by DHS regions
ensgene = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/motif_analysis_encode/ensemble_mm10_TSS.txt', sep='\t', header=FALSE)
jj = which(ensgene[,3]=='-')
kk = which(ensgene[,3]=='+')
tss1 = cbind(ensgene[kk, c(2, 4, 4, 1)], rep(0, length(kk)), ensgene[kk, 3])
tss2 = cbind(ensgene[jj, c(2, 5, 5, 1)], rep(0, length(jj)), ensgene[jj, 3])
colnames(tss1) = c('chr', 'start', 'end', 'name', 'score', 'strand')
colnames(tss2) = c('chr', 'start', 'end', 'name', 'score', 'strand')
tss = rbind(tss1, tss2)
colnames(tss) = c('chr', 'start', 'end', 'name', 'score', 'strand')

tss = tss[which(tss[,2]>0), ]
write.table(tss, file='/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/motif_analysis_encode/TSS_all_mm10.bed', sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)
cmd1 = "bedtools sort -i TSS_all_mm10.bed > TSS_all_mm10_sorted.bed"
system(cmd1) ## filter TSS with DHS regions
cmd2 = "bedtools intersect -a TSS_all_mm10_sorted.bed -b encode_broadpeaks_mm10.bed > TSS_within_dhs_encode.bed"
system(cmd2) ## filter TSS with DHS regions

##########################
#### transcriptional outputs (rhythmic introns) of expressed genes from RNA-seq
##########################
source('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/f24_modified_1.0.r')

aa = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/WT_RF_Intron_Exon_RFP.txt',sep='\t', header=TRUE)

mrna = aa[, -c(3:26, 51:74)]
mrna = mrna[, c(2:14)]
intron = aa[, c(2:14)]
colnames(mrna) = c('gene', paste('mRNA.ZT', (c(0:11)*4+2), sep=''))
colnames(intron) = c('gene', paste('intron.ZT', (c(0:11)*4+2), sep=''))
#mrna = data.frame(mrna$gene, log2(as.matrix(mrna[,-1])))
#colnames(mrna) = c('gene', paste('mRNA.ZT', (c(0:11)*4), sep=''))
require(lattice)
require(ggplot2)
#pairs(mrna[,-1])
data.mrna = mrna[,-1]
tt = c(0:11)*4+2
res.mrna = t(apply(data.mrna, 1, f24_R2_alt2, t=tt))
res.mrna[,4] = t(apply(2^data.mrna, 1, f24_R2_alt2, t=tt))[,4]
colnames(res.mrna) = paste(colnames(res.mrna), '.mrna', sep='')
mrna = data.frame(mrna, res.mrna, stringsAsFactors=FALSE)

data = intron[,-1]
tt = c(0:11)*4+2
res = t(apply(data, 1, f24_R2_alt2, t=tt))
res[,4] = t(apply(2^data, 1, f24_R2_alt2, t=tt))[,4]
colnames(res) = paste(colnames(res), '.intron', sep='')
intron = data.frame(intron, res, stringsAsFactors=FALSE)

### select expressed genes
hist(mrna$mean.mrna, breaks=50);abline(v=0, col='red')
cutoff.expressed = 0
jj = which(mrna$mean.mrna>cutoff.expressed)
intron = intron[jj,]
intron$qv.intron = qvals(intron$pval.intron)
cutoff.qv = 0.1
sel = which(intron$qv.intron<cutoff.qv)
intron.sel = intron[sel,]

examples = c('Arntl', 'Clock', 'Per1', 'Per2', 'Per3', 'Cry1', 'Cry2', 'Rora', 'Rorc', 'Nr1d1', 'Nr1d2', 'Dbp', 'Hlf', 'Tef', 'Nfil3', 'Bhlhe40', 'Nampt', 'Parp1', 'Prkaa2', 'Prkag1', 'Gsk3a', 'Gsk3b', 'Csnk1d')
mm = match(examples, intron.sel$gene)

save(intron.sel, file='Rdata/Intron_sel_rhythmic.Rdata')

#####################################
#### find the associated DHS regions
#####################################
load(file='Rdata/Intron_sel_rhythmic.Rdata')
aa = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/WT_RF_Intron_Exon_RFP.txt',sep='\t', header=TRUE)
tss = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/motif_analysis_encode/TSS_within_dhs_encode.bed', sep='\t', header=FALSE, as.is=c(1))
dhs.all = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/motif_analysis_encode/encode_broadpeaks_mm10.bed', sep='\t', header=FALSE, as.is=c(1))
ensgene = read.table('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/motif_analysis_encode/ensemble_mm10_TSS.txt', sep='\t', header=FALSE)
save(intron.sel, aa, tss, dhs.all, ensgene, file='Rdata/Tables4preparation_elasticnet.Rdata')

window.size = 25000;
version = '_encode_25kb'

Prepare.table4elasticnet = TRUE
if(Prepare.table4elasticnet)
{
    load(file='Rdata/Tables4preparation_elasticnet.Rdata')
    annotate.dhs = function(x)
    {
        annot = paste(as.character(x[1]), ':', x[2], '-', x[3], sep='')
    }

    keep = rep(NA, nrow(intron.sel))
    for(n in 1:length(intron.sel$gene))
    {
        cat(n, '\n')
        gg = intron.sel$gene[n]
        g2 = aa[which(aa[,2]==gg), 1]
        transcripts = ensgene[match(g2, ensgene[,6]), 1]

        kk = match(transcripts, tss[,4])
        kk = kk[which(!is.na(kk))]
        if(length(kk)>0)
        {
            test = c()
            for(index in kk)
            {
                ii = which(dhs.all[,1]==tss[index, 1] & (dhs.all[, 2]> (tss[index, 2]-window.size) & dhs.all[, 3]< (tss[index, 2]+window.size)))
                if(length(ii)==1) test = c(test, annotate.dhs(dhs.all[ii,]))
                if(length(ii)>1) test = c(test, unlist(apply(dhs.all[ii,], 1, annotate.dhs)))
            }
            if(length(test)>0) keep[n] = paste(test, collapse=';', sep='')
        }
    }

    gene2dhs = data.frame(intron.sel$gene, keep, stringsAsFactors=FALSE)
    kk = which(!is.na(gene2dhs[,2])==TRUE)
    intron.sel = intron.sel[kk, ]
    gene2dhs = gene2dhs[kk,]

    cat('nb of intron ', nrow(intron.sel), '\n')

    dhs.list = gene2dhs[,2]
    dhs.list = unique(unlist(strsplit(as.character(dhs.list), ';')))

    cmd0 = "rm dhs_list.txt"
    system(cmd0)

    write.table(dhs.list, file='/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/motif_analysis_encode/dhs_list.txt', sep='\t', col.names=FALSE, row.names=FALSE, quote=FALSE)

    cmd1 = "LANGUAGE=en_US.UTF-8; LC_ALL=en_US.UTF-8; LANG=en_US.UTF-8; LC_TYPE=en_US.UTF-8"
    system(cmd1)

    cmd = "perl fimoParserToMatrixCount_selregion.pl -fimo encode_dhs_motif_scanning_all_curation_p_0.0001.txt > encode_dhs_motif_oc_selected_region.txt"
    system(cmd)

    #### DHS motif occurrence
    #Data_table<-read.table("encode_dhs_motif_occurence.txt",header=TRUE)

    ### pval by defaut
    Curation.v.expressed.tfs.pwms = TRUE
    {
        motif_occurence<-read.table("encode_dhs_motif_oc_selected_region.txt",check.names = FALSE,sep=" ",header=TRUE)  ### input motif_scan mat of DHS by curated PWMs
        mot_oc<-data.matrix(motif_occurence[,2:dim(motif_occurence)[2]])
        rownames(mot_oc)<-motif_occurence[,1]

        save(intron.sel, dhs.list, gene2dhs, mot_oc, file=paste('Rdata/DHS_matrix_oc_pwms_p_defaut_intron_sel', version, '.Rdata', sep=''))
    }

    ### contruct matrix of oc for selected genes
    mm = match(intron.sel$gene, gene2dhs[,1])
    gene2dhs = gene2dhs[mm, ]

    miss = c()
    mot_intron = matrix(NA, nrow = nrow(intron.sel), ncol = nrow(mot_oc))
    for(n in 1:nrow(intron.sel))
    #for(n in 1:10)
    {
        print(nrow(intron.sel)-n)

        dhs.sel = gene2dhs[n, 2]
        dhs.sel = unlist(strsplit(as.character(dhs.sel), ';'))

        kk = match(dhs.sel, colnames(mot_oc))
        kk = kk[which(!is.na(kk))]
        if(length(kk)>1) {
            sum_oc = apply(mot_oc[, kk], 1, sum)
            mot_intron[n, ] = sum_oc;
        }else{
            if(length(kk)==1) {
                mot_intron[n, ] = mot_oc[,kk]
            }else{
                miss = c(miss, n)
            }
        }

    }

    colnames(mot_intron) = rownames(mot_oc)

    if(length(miss)>0) {
        intron.sel = intron.sel[-miss, ]
        mot_intron = mot_intron[-miss, ]
    }

    #colnames(mot_pol2) = rownames(pol2.sel)
    Check.Matrix = FALSE
    if(Check.Matrix)
    {
        rss = apply(mot_intron, 1, sum)
        css = apply(mot_intron, 2, sum)

        xx = mot_pol2[which(rss>0),]
        res.pol2.sel = res.pol2.sel[which(rss>0),]
        xx = xx[,which(css>0)]
        mot_pol2 = xx
    }
    save(mot_intron, intron.sel, file=paste('Rdata/Elastic_CM_library_encode_dhs_intron', version, '.Rdata', sep=''))
}


######
#### Elastic-net
#version = '_encode_0.5kb'
#version = '_encode_1kb'
#version = '_encode_2.5kb'
version = '_encode_5kb'
#version = '_encode_10kb'
#version = '_encode_25kb'
alpha = 0.1
binary = 1;
binary.matrix = binary== 0

load = TRUE
if(load) load(file=paste('Rdata/Elastic_CM_library_encode_dhs_intron', version, '.Rdata', sep=''))
x=as.matrix(mot_intron)
#x = as.matrix(x>0)
y1=intron.sel$amp.intron*cos(2*pi*intron.sel$phase.intron/24)
y2=intron.sel$amp.intron*sin(2*pi*intron.sel$phase.intron/24)
y=cbind(y1,y2)

cat('nb of intron ', nrow(x), '\n')

Filter.intron = TRUE
if(Filter.intron)
{
    jj = which(intron.sel$relamp.intron>0.2)
    x = x[jj,]
    y = y[jj,]
}

require(glmnet)
filter.motifs = TRUE
if(filter.motifs)
{
    kk = match(c('DBP_p2', 'EP300_p2', 'CTCF_p2', 'HOX_A5_B5_p2', 'HOXA9_MEIS1_p2', "HOX_A4_D4_p2", "PPARD_f1", "RORG_f1", "CEBPG_si", "ZFHX3_f1", "MCR_f1" , "TEF_f1", "ZBT7A_f1", "IRF9_f1", "ATF1_si" , "FOXA3_f1", 'NR1D1_f1', "TBX3_f1", "NR2C1_si", "FOXD3_p2", "FOXL1_p2", "FOXN1_p2", "FOXQ1_p2","FOX_C1_C2_p2", "FOX_D1_D2_p2","FOX_F1_F2_J1_p2", "FOX_I1_J2_p2", "V_FOX_Q2", "NANOG_p2", "V_ETV3_01"), colnames(x))

    x = x[,-kk]
    rss = apply(x, 1, sum)
    jj = which(rss==0)
    if(length(jj)>0)
    {
        x = x[-jj,];
        y = y[-jj,];
    }
}

#x=as.matrix(mat_oc_ordered[sel,]) >0
intercept=0
standardize=TRUE ### standardize matrix of motif occurrence makes more sense because the absolute number of motif occurrence is not precise.
standardize.response=FALSE

if(binary.matrix){x = x >0; standardize=FALSE}

### use Cross-validation to select tuning paprameter
cv.fit=cv.glmnet(x,y,family='mgaussian',grouped=TRUE, alpha=alpha, nlambda=500, standardize=standardize, standardize.response=standardize.response, intercept=intercept)
plot(cv.fit)
#cv.fit$lambda

optimal = which(cv.fit$lambda==cv.fit$lambda.min)
#optimal = which(cv.fit$lambda==cv.fit$lambda.1se)
fit=glmnet(x,y,alpha=alpha, lambda=cv.fit$lambda,family='mgaussian', type.multinomial=c("grouped"), standardize=standardize, standardize.response=standardize.response, intercept=intercept)
#fit=glmnet(x,y,alpha=alpha, family='mgaussian', type.multinomial=c("grouped"), standardize=standardize, standardize.response=standardize.response, intercept=intercept)


colnames(x)[which(fit$beta[[1]][,optimal]!=0)]
#colnames(x)[which(fit$beta[[2]][,optimal]!=0)]
#coef(cv.fit)
#par(mfrow=c(1,1))
fit.pol2 = fit
optimal.pol2 = optimal

###
### Motifs inferred to be involved in circadian transcription networks
###
kk = c(which(fit.pol2$beta[[1]][,optimal.pol2]!=0),which(fit.pol2$beta[[2]][,optimal.pol2]!=0))
kk = unique(kk)
names.pol2 = colnames(x)[kk]
a.pol2 = fit.pol2$beta[[1]][,optimal.pol2][kk]
b.pol2 = fit.pol2$beta[[2]][,optimal.pol2][kk]
aa = a.pol2
bb = b.pol2
period = 24
keep.pol2 = c()

for(n in 1:length(names.pol2))
{
    phase=period/(2*pi)*atan2(bb[n], aa[n])
    if(phase<0) phase=phase+period
    if(phase>period) phase=phase-period
    amp = sqrt(aa[n]^2+bb[n]^2)
    keep.pol2 = rbind(keep.pol2, c(phase,amp))
}
rownames(keep.pol2) = names.pol2
colnames(keep.pol2) = c('phase','ampl')

infer = keep.pol2 # infered motifs by elastic-net

save(infer, file=paste('Rdata/Motifs_inferred_CV_elastic_net_alpha_0.1', version, '.Rdata', sep=''))
#save(infer, file=paste('Rdata/Motifs_inferred_CV_elastic_net_alpha_0.05', version, '_occurrence_matrix_exclude_Hocomoco_5kb2TSS.Rdata', sep=''))
#### compare with old version
load(file='/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/Elastic-net_analysis/Motifs_inferred_CV_elastic_net_alpha_0.1.Rdata')
infer0 = infer
load(file=paste('Rdata/Motifs_inferred_CV_elastic_net_alpha_0.1', version, '.Rdata', sep=''))
load(file='/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/motif_analysis_encode/Rdata/motifs_used_ref.Rdata')
motifs.used = c(motifs.used, 'E2F1_2_3_4_5_p2')

cat(paste('COMMON motifs :', length(intersect(rownames(infer), rownames(infer0))), sep=''), '\n' );intersect(rownames(infer), rownames(infer0));

cat(paste('NEW motifs :', length(setdiff(rownames(infer), rownames(infer0))), sep=''), '\n' );setdiff(rownames(infer), rownames(infer0));

cat(paste('LOST motifs :', length(setdiff(rownames(infer0), rownames(infer))), sep=''), '\n' ); setdiff(rownames(infer0), rownames(infer))

cat(paste('LOST USED motifs :', length(setdiff(motifs.used, rownames(infer))), sep=''), '\n' ); setdiff(motifs.used, rownames(infer))

##### Plot phases of inferred Motifs
motif.amp = 7.5;
phase.m = infer[,1]

o1 = order(phase.m)
phase.m = phase.m[o1]
motif.names = rownames(infer)[o1]

amp = motif.amp
motif.a = amp*cos(2*pi/24*phase.m)
motif.b = amp*sin(2*pi/24*phase.m)
CC = (motif.a -1i*motif.b) * exp(1i*pi/2)
motif.aa = Re(CC)
motif.bb = Im(CC)

amp = motif.amp + 0.25
motif.a = amp*cos(2*pi/24*phase.m)
motif.b = amp*sin(2*pi/24*phase.m)
CC = (motif.a -1i*motif.b) * exp(1i*pi/2)
motif.txt1.aa = Re(CC)
motif.txt1.bb = Im(CC)

amp = motif.amp-0.25
motif.a = amp*cos(2*pi/24*phase.m)
motif.b = amp*sin(2*pi/24*phase.m)
CC = (motif.a -1i*motif.b) * exp(1i*pi/2)
motif.txt2.aa = Re(CC)
motif.txt2.bb = Im(CC)

rmm=max(abs(c(motif.aa,motif.bb)))+2.0
rr=c(-rmm,rmm)
xlim = rr
ylim = rr

pdfname = paste('myplots/test/Motif_activities_elastic_net_introns', version, '_alpha_', alpha, '.pdf', sep='')
pdf(pdfname, width=6, height=6)
par(cex = 0.5, las = 1, mgp = c(1.6,0.5,0), mar = c(3,3,5,0.8)+0.1, tcl = -0.3)
plot(motif.aa, motif.bb, main=paste(version, '  alpha = ', alpha, ' binary =', binary, sep=''), type='n', xlim=xlim, ylim=ylim, axes=F, xlab='', ylab='', lwd=2, pch=23, col='black',bg='green',cex=2.0)

abline(v=0,h=0, col='darkgray',lwd=2.0)
phi=seq(0,2*pi,len=1000)
#lines(rm*cos(phi), rm*sin(phi), col='darkgray', lwd=2)
lines((motif.amp-0.5)*cos(phi), (motif.amp-0.5)*sin(phi), col='darkgray', lwd=2)
rainbow = rainbow(length(motif.names),s = 0.85, v = 0.85)
for(n in 1:length(motif.names))
{
    points(motif.aa[n], motif.bb[n], pch=21, cex=1.5, col='black',bg=rainbow[n])

    if(phase.m[n]<=6) srt = 90-phase.m[n]/24*360;
    if(phase.m[n]>6 & phase.m[n]<=12) srt = 90-phase.m[n]/24*360;
    if(phase.m[n]>12 & phase.m[n]<=18) srt = 270-phase.m[n]/24*360;
    if(phase.m[n]>=18) srt = 270-phase.m[n]/24*360;

    if(phase.m[n]<=12) {
        if(n%%2==1){
            adj = 0;
            pos=1;
            text(motif.txt1.aa[n], motif.txt1.bb[n], motif.names[n], cex=0.7, col='black', srt=srt, adj=adj)
        }else{
            adj = 1;
            text(motif.txt2.aa[n], motif.txt2.bb[n], motif.names[n], cex=0.7, col='black', srt=srt, adj=adj)
        }

    }else{
        if(n%%2==1)
        {
            adj = 1;
            pos=3.0
            text(motif.txt1.aa[n], motif.txt1.bb[n], motif.names[n], cex=0.7, col='black', srt=srt, adj=adj)
        }else{
            adj = 0;
            pos=1;
            text(motif.txt2.aa[n], motif.txt2.bb[n], motif.names[n], cex=0.7, col='black', srt=srt, adj=adj)
        }
    }
    #legend(x = 0.9*xlim[2],y = ylim[2]*(30-n)/30, legend = motif.names[n], cex=0.6, col='black', pt.bg = rainbow[n], pch=21, bty = 'n' )
}
dev.off()


Test.optimal.parameters
if(Test.optimal.parameters)
{
    pdfname = paste('myplots/test/Motif_activities_inferred_Pol2_elastic_net_version_alpha_binary_occurence_matrix_all.pdf', sep='')
    pdf(pdfname, width=6, height=6)

    for(version in c('_encode_1kb', '_encode_2.5kb', '_encode_5kb', '_encode_10kb', '_encode_25kb'))
    {
        for(binary in c(0, 1))
        {
            for(alpha in c(0.05, 0.1))
            {
                binary.matrix = binary==


            }
        }
    }
    dev.off()
}

#######
### Annotate motifs
#######
load(file=paste('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/motif_analysis_encode/Rdata/Motifs_inferred_CV_elastic_net_alpha_0.05_encode_5kb_occurrence_matrix_exclude_Hocomoco_5kb2TSS.Rdata'))

version = '_v3'
filename = paste('/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/Annotations_TFs/TFs_Motifs_Mapping_all_Database_Selection_weng_suissregulon_jaspar_hocomoco_rhythmic_transfac_final', version, '.txt', sep='')
mapping = read.table(file=filename, sep='\t', header=TRUE)
which(is.na(mapping$Motifs))

mapping.motif.tf = mapping

Motifs = c()
index = c()
for(n in 1:nrow(mapping))
{
    test = mapping$Motifs[n]
    test = unlist(strsplit(as.character(test), '[|]'))
    Motifs = c(Motifs, test)
    index = c(index, rep(n, length(test)))
}
xx = mapping[index,]
xx$Motifs = Motifs
motif.tf = data.frame(xx$TFs.pool, xx$Motifs, xx$Sources, xx$TFs.nuclear, xx$TFs.nuclear.init, xx$phase.nuclear, xx$qv.nuclear, stringsAsFactors=FALSE)
#colnames(motif.tf) = unlist(strsplit(colnames(motif.tf), '[.]'))[1]
colnames(motif.tf) = c('TFs.pool', 'Motifs', 'Sources', 'TFs.nuclear', 'TFs.nuclear.init', 'phase', 'qv')
motif.tf = motif.tf[order(motif.tf[, ncol(motif.tf)]), ]

index.motif = c()
tf.mapping = c()
phase.tf.mapping = c()
pval.tf.mapping = c()
qv.tf.mapping = c()

for(n in 1:nrow(infer))
{
    ### inferred motifs
    kk = which(mapping.motif.tf$Motifs==rownames(infer)[n])
    if(length(kk)==0)
    {
        print(rownames(infer)[n])
        index.motif = c(index.motif, n);
        tf.mapping = c(tf.mapping, NA)
        phase.tf.mapping = c(phase.tf.mapping, NA)
        pval.tf.mapping = c(pval.tf.mapping, NA)
        qv.tf.mapping = c(qv.tf.mapping, NA)
    }else{

        dd = mapping.motif.tf[kk,]
        index.motif = c(index.motif, rep(n, length(kk)));

        tf.mapping = c(tf.mapping, as.character(dd$TFs.pool))
        phase.tf.mapping = c(phase.tf.mapping, dd$phase.nuclear)
        pval.tf.mapping = c(pval.tf.mapping, dd$pval.nuclear)
        qv.tf.mapping = c(qv.tf.mapping, dd$qv.nuclear)
    }

}

res = data.frame(rownames(infer)[index.motif], rep(NA, length(index.motif)), infer[index.motif, ], tf.mapping, phase.tf.mapping, pval.tf.mapping, qv.tf.mapping, stringsAsFactors=FALSE)
colnames(res) = c('Motifs.inferred', 'Motif.names.modif', 'Phase.motifs.inferred', 'Amp.motifs.inferred', 'TFs.mapping', 'phases.tfs.mapping','pval.tfs.mapping','qv.tfs.mapping')

write.table(res, file='/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/motif_analysis_encode/TFs_motifs_acvity_mapping_dhs_encode_alpha_0.05_matrix_occurrence_5kb2TSS.txt', row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE)

######
### Modify motif names and also assigne TFs to some motifs
######
res = read.table(file='/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/motif_analysis_encode/TFs_motifs_acvity_mapping_dhs_encode_alpha_0.05_matrix_occurrence_5kb2TSS.txt', header=TRUE, sep='\t', as.is=c(1,2,5))
xx = read.table(file='/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Tables_DATA/Motifs_inferred_tfs_mapping_curate_alpha_0.1_manual_v3.txt', sep='\t', header=TRUE, as.is = c(1,2,5))
length(unique(xx[,1]))

res = res[, c(1:5)]

yy = res
motifs = unique(yy[,1])

add = c()
for(test in motifs)
{
    kk = which(xx[,1]==test)
    jj = which(res[,1]==test)

    if(length(kk)>0)
    {
        name.m = xx[kk, 2]
        name.m = name.m[which(!is.na(name.m))]
        if(length(name.m)==1) yy[jj, 2] = name.m;

        tfs.map = xx[kk, 5]
        if(length(setdiff(tfs.map, yy[jj,5])))
        {
            index.add = match(setdiff(tfs.map, yy[jj,5]), xx[kk, 5])
            for(iindex in index.add)
            {
                add = rbind(add, c((yy[jj[1], c(1:4)]), xx[kk[iindex], 5]))
            }
        }
    }
}

add = data.frame(matrix(unlist(add), nrow=21, byrow=FALSE))
colnames(add) = colnames(yy)
yy = rbind(yy, add)

yy = yy[order(yy[,1]),]
write.table(yy, file='/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Transcription_network/motif_analysis_encode/TFs_motifs_acvity_mapping_dhs_encode_alpha_0.05_matrix_occurrence_5kb2TSS_v1.txt', row.names=FALSE, col.names=TRUE, sep='\t', quote=FALSE)
