#preprocessing.sh                             #command:
snakemake
snakemake all
# genomestarinit=$(date +%s.%N)
# snakemake -j 4 -R star_genome_generate
# genomestarend=$(date +%s.%N)
alignmentstarinit=$(date +%s.%N)
snakemake -j 4 -R star_alignment
alignmentstarend=$(date +%s.%N)
kallistoindexinit=$(date +%s.%N)
snakemake -j 4 -R kallisto_index
kallistoindexend=$(date +%s.%N)
kallistoquantinit=$(date +%s.%N)
snakemake -j 4 -R kallisto_pseudo
kallistoquantend=$(date +%s.%N)
samtoolsviewinit=$(date +%s.%N)
snakemake -j 4 -R samtools_view
samtoolsviewend=$(date +%s.%N)
samtoolsortinit=$(date +%s.%N)
snakemake -j 4 -R samtools_sort
samtoolsortend=$(date +%s.%N)
samtoolsindexinit=$(date +%s.%N)
snakemake -j 4 -R samtools_index
samtoolindexend=$(date +%s.%N)
samtoolpercentageinit=$(date +%s.%N)
snakemake -j 4 -R samtools_percentage
samtoolpercentageend=$(date +%s.%N)

# rtstargenome=$(python -c "print(${genomestarend} - ${genomestarinit})")
rtstaralignment=$(python -c "print(${alignmentstarend} - ${alignmentstarinit})")
rtstar=$(python -c "print(${alignmentstarend} - ${genomestarinit})")
rtkallistoindex=$(python -c "print(${kallistoindexend} - ${kallistoindexinit})")
rtkallistoquant=$(python -c "print(${kallistoquantend} - ${kallistoquantinit})")
rtkallisto=$(python -c "print(${kallistoquantend} - ${kallistoindexinit})")
rtstview=$(python -c "print(${samtoolsviewend} - ${samtoolsviewinit})")
rtstsort=$(python -c "print(${samtoolsortend} - ${samtoolsortinit})")
rtstindex=$(python -c "print(${samtoolindexend} - ${samtoolsindexinit})")
rtstpercentage=$(python -c "print(${samtoolpercentageend} - ${samtoolpercentageinit})")
rtst=$(python -c "print(${samtoolpercentageend} - ${samtoolsviewinit})")
rttotal=$(python -c "print(${samtoolpercentageend} - ${alignmentstarinit})")

echo "Runtime of genome generate of star = $rtstargenome"
echo "Runtime of alignment of star = $rtstaralignment"
echo "Runtime total of star = $rtstar"
echo "Runtime of index generation of kallisto = $rtkallistoindex"
echo "Runtime of pseudo kallisto = $rtkallistoquant"
echo "Runtime total of kallisto = $rtkallisto"
echo "Runtime of view of samtools = $rtstview"
echo "Runtime of sort of samtools = $rtstsort"
echo "Runtime of index of samtools = $rtstindex"
echo "Runtime of percentage of samtools = $rtstpercentage"
echo "Runtime total of samtools = $rtst"
echo "Runtime total = $rttotal"
