
##############################################
#########	Aspergillus nidulans		######
##############################################

#bowtie2 mapping
bowtie2 -p 6 --local  -x <bt2-idx> -1 <> -2 <> | samtools view -bS - | samtools sort  -O bam -o <>.bam


#Index and alignment stats
{
for i in `cat sample_ANidulans.list`
do
cd $i
bash /home/lakhanp/scripts/ChIPseq_scripts/template_ChIPseq_config.sh -c /home/lakhanp/database/reference_genomes.yaml -o A_nidulans --polII >> generalJob.sh
sed "s/SAMPLE_ID/$i/g" /home/lakhanp/scripts/ChIPseq_scripts/template_ChIPseq_process.sh >> generalJob.sh
cd ..
done
}


## Calculate GC bias
{
for i in `cat sample_ANidulans.list`
do
cd $i
printf "##Compute GC bias
computeGCBias -b %s_bt2.bam --effectiveGenomeSize 29850950 -g /home/lakhanp/database/A_nidulans_FGSC_A4/reference/A_nidulans_FGSC_A4_version_s10-m04-r03_chromosomes.2bit -l 60 --GCbiasFrequenciesFile %s_GCbiasFreq.txt --biasPlot %s_GCbias.png 
error_exit \$? \n\n" $i $i $i >> generalJob.sh
cd ..
done
}


#### Generate profile matrix
##{
##for i in `cat sample_ANidulans.list`
##do
##cd $i
##printf "##generate profile matris (-2kb == normalized(geneBody) to 2kb == +1kb)
##computeMatrix scale-regions -S %s_normalized.bw  -R /home/lakhanp/database/A_nidulans_FGSC_A4/annotation/A_nidulans_FGSC_A4_version_s10-m04-r03_CDS_Unique.bed -m 2000 -b 2000 -a 1000 --numberOfProcessors 2  --outFileName %s_normalized_profile.tab.gz
##error_exit \$? \n\n" $i $i >> generalJob.sh
##cd ..
##done
##}
####


## Print control information and peak type (narrow/broad) information
{
## IMP
## use the printf information from the excel file
for i in `cat sample_tf_macs2.list`
do cd $i
printf "peakType=\'narrow\'\n\n" >> generalJob.sh
cd ..
done
}


## MACS2 peak calling: use template
{
for i in `cat sample_tf_macs2.list`
do cd $i
sed "s/SAMPLE_ID/$i/g" /home/lakhanp/scripts/ChIPseq_scripts/template_ChIPseq_macs2.sh >> generalJob.sh
cd ..
done
}







## copy data to local
for i in `cat sample_ANidulans.list`
do 
cd $i
cp ${i}_normalized.bw ${i}_normalized_profile.tab.gz ${i}_polii_expr.tab.rel.mat ../localCopy/
cd ..
done

for i in `cat tf_AN.list`; do cd $i; 
cp macs2_*/${i}*{.narrowPeak,.broadPeak,.tab} ../localCopy/
cd ..
done


## Rscript E:\Chris_UM\Codes\Shuhui_SM_OE_project\kmeans_all.R > kmeasns_all.log


##############################################
## read count for polII data: deeptools
multiBamSummary BED-file \
--BED /home/lakhanp/database/A_nidulans_FGSC_A4/annotation/AN_genes_for_polII.bed \
--bamfiles /home/lakhanp/Analysis/23_Shuhui_SM_OE_project/mapping/*_polII_*/*bam \
--smartLabels -p 12 --transcript_id_designator transcript_id \
-out raw_count.polII.deeptools.npz --outRawCounts raw_count.polII.deeptools.tab


## read count for polII data: bedtools
bedtools multicov -bams /home/lakhanp/Analysis/23_Shuhui_SM_OE_project/mapping/*_polII_*/*bam \
-bed /home/lakhanp/database/A_nidulans_FGSC_A4/annotation/AN_genesForPolII.bed > raw_count.polII.bedtools.tab



###################################################
##### bigwig correlation using deeptools     ######
###################################################
###############
## generate the counts using deeptools multiBigwigSummary: gene body
multiBigwigSummary BED-file \
--bwfiles /home/lakhanp/Analysis/23_Shuhui_SM_OE_project/mapping/*/*_normalized.bw \
/home/lakhanp/Analysis/19_ChIPMix_process/CL2019_ChIPmix_46/mapping/AN*_xyl_*/*_normalized.bw \
/home/lakhanp/Analysis/19_ChIPMix_process/CL2019_ChIPmix_55/mapping/*/*_normalized.bw \
/home/lakhanp/Analysis/19_ChIPMix_process/CL2019_ChIPmix_56/mapping/{AN,WT,hepAdel}*/*_normalized.bw \
/home/lakhanp/Analysis/19_ChIPMix_process/CL2019_ChIPmix_60/mapping/*/*_normalized.bw \
--BED /home/lakhanp/database/A_nidulans_FGSC_A4/annotation/A_nidulans_FGSC_A4_version_s10-m04-r03_CDS_Unique.bed \
-out An_bigwig_summary_gene.npz --outRawCounts An_bigwig_summary_gene.tab -p 12


## plotCorrelation for all samples: pearson
plotCorrelation --corData An_bigwig_summary_gene.npz --corMethod pearson --skipZeros --plotTitle "pearson Correlation at gene body" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o An_bigwig_pearson_gene.png --outFileCorMatrix An_bigwig_pearson_gene.tab --plotHeight 100 --plotWidth 100

## plotCorrelation for all samples: spearman
plotCorrelation --corData An_bigwig_summary_gene.npz --corMethod spearman --skipZeros --plotTitle "spearman Correlation at gene body" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o An_bigwig_spearman_gene.png --outFileCorMatrix An_bigwig_spearman_gene.tab --plotHeight 100 --plotWidth 100

###############
## generate the counts using deeptools multiBigwigSummary: promoter
multiBigwigSummary BED-file \
--bwfiles /home/lakhanp/Analysis/23_Shuhui_SM_OE_project/mapping/*/*_normalized.bw \
/home/lakhanp/Analysis/19_ChIPMix_process/CL2019_ChIPmix_46/mapping/AN*_xyl_*/*_normalized.bw \
/home/lakhanp/Analysis/19_ChIPMix_process/CL2019_ChIPmix_55/mapping/*/*_normalized.bw \
/home/lakhanp/Analysis/19_ChIPMix_process/CL2019_ChIPmix_56/mapping/{AN,WT,hepAdel}*/*_normalized.bw \
/home/lakhanp/Analysis/19_ChIPMix_process/CL2019_ChIPmix_60/mapping/*/*_normalized.bw \
--BED /home/lakhanp/database/A_nidulans_FGSC_A4/annotation/A_nidulans_FGSC_A4_version_s10-m04-r03_CDS_upstream500.bed \
-out An_bigwig_summary_promoter.npz --outRawCounts An_bigwig_summary_promoter.tab -p 12


## plotCorrelation for all samples: pearson
plotCorrelation --corData An_bigwig_summary_promoter.npz --corMethod pearson --skipZeros --plotTitle "pearson Correlation at promoter region" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o An_bigwig_pearson_promoter.png --outFileCorMatrix An_bigwig_pearson_promoter.tab --plotHeight 100 --plotWidth 100

## plotCorrelation for all samples: spearman
plotCorrelation --corData An_bigwig_summary_promoter.npz --corMethod spearman --skipZeros --plotTitle "spearman Correlation at promoter region" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o An_bigwig_spearman_promoter.png --outFileCorMatrix An_bigwig_spearman_promoter.tab --plotHeight 100 --plotWidth 100


################
## generate the counts using deeptools multiBigwigSummary: 1kb bins
multiBigwigSummary bins \
--bwfiles /home/lakhanp/Analysis/23_Shuhui_SM_OE_project/mapping/*/*_normalized.bw \
/home/lakhanp/Analysis/19_ChIPMix_process/CL2019_ChIPmix_46/mapping/AN*_xyl_*/*_normalized.bw \
/home/lakhanp/Analysis/19_ChIPMix_process/CL2019_ChIPmix_55/mapping/*/*_normalized.bw \
/home/lakhanp/Analysis/19_ChIPMix_process/CL2019_ChIPmix_56/mapping/{AN,WT,hepAdel}*/*_normalized.bw \
/home/lakhanp/Analysis/19_ChIPMix_process/CL2019_ChIPmix_60/mapping/*/*_normalized.bw \
--binSize 1000 -out An_bigwig_summary_bin.npz --outRawCounts An_bigwig_summary_bin.tab -p 12

## plotCorrelation for all samples: pearson
plotCorrelation --corData An_bigwig_summary_bin.npz --corMethod pearson --skipZeros --plotTitle "pearson Correlation at 1kb bin regions" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o An_bigwig_pearson_bin.png --outFileCorMatrix An_bigwig_pearson_bin.tab --plotHeight 100 --plotWidth 100

## plotCorrelation for all samples: spearman
plotCorrelation --corData An_bigwig_summary_bin.npz --corMethod spearman --skipZeros --plotTitle "spearman Correlation at 1kb bin regions" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o An_bigwig_spearman_bin.png --outFileCorMatrix An_bigwig_spearman_bin.tab --plotHeight 100 --plotWidth 100

######################################
####   Motif enrichment analysis  ####
######################################

{
## meme-chip analysis: with control
while IFS=$'\t' read -r name fasta neg
do
	printf "## meme-chip: ${name}\n"
	printf "## meme-chip: ${name}
	meme-chip -order 1 -meme-minw 5 -meme-maxw 30 -meme-nmotifs 5 -meme-mod zoops \
	-db /home/lakhanp/tools/meme_motif_databases_12.19/JASPAR/JASPAR2018_CORE_fungi_non-redundant.meme \
	-desc ${name} -oc ${name} -neg ${neg} ${fasta} < /dev/null
	"
	printf "## done...\n\n"
done < memechip_de_conf.tab
}

{
## meme-chip analysis: without control
while IFS=$'\t' read -r name fasta
do
	printf "## meme-chip: ${name}\n"
	#printf "## meme-chip: ${line[0]}
	meme-chip -order 1 -meme-minw 5 -meme-maxw 30 -meme-nmotifs 5 -meme-mod anr -meme-p 12 \
	-db /home/lakhanp/tools/meme_motif_databases_12.19/JASPAR/JASPAR2018_CORE_fungi_non-redundant.meme \
	-desc ${name} -oc ${name} ${fasta} < /dev/null
	#"
	printf "## done...\n\n"
done < memechip_conf.tab
}



##fimo motif scanning in polII DEG promoter
fimo -oc fimo.AN0153_OE_vs_WT_DEG --thresh 0.001 --bfile /home/lakhanp/database/A_nidulans_FGSC_A4/reference/promoter_sequences/A_nidulans_promoters.500_TSS_100.meme_background.m2.model AN0153_motif.meme AN0153_OE_vs_WT.DEG_promoter.fasta


##########################
####   AflR analysis  ####
##########################
##fimo motif scanning in whole genome
fimo -oc fimo.AN7820_genome --thresh 0.001 --bfile /home/lakhanp/database/A_nidulans_FGSC_A4/reference/promoter_sequences/A_nidulans_promoters.500_TSS_100.meme_background.m2.model AN7820_motifs.meme /home/lakhanp/database/A_nidulans_FGSC_A4/reference/A_nidulans_FGSC_A4_version_s10-m04-r03_chromosomes.fasta






