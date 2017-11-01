Input: vcf file must be annotated by SnpSift/SnpEff with information from dbNSFP.  Filtering of vcf file prior to input is recommended - i.e. mutation effect, quality filter, removal of artifact variants.

Parameters.txt file
Mark choices with an "x"
1. Under "Mutation_Effects" section, choose mutation effects from SnpEff to filter which variants WILL be included in the analysis.
2. Under "Ancestry" section, choose which annotated ancestries to include for minor allele frequency (MAF) cutoff.  
3. Under "Genotypes" section, choose type of mutation(s) you would like to consider: CompoundHeterozygotes, Homozygous, TwoMutations, OneMutation.  Input data must be phased to choose CompoundHeterozygotes.
4. Under "Search_Level", choose either "Gene" or "Transcript" for which locus to quantify variants.
5. Under "Minor_Allele_Frequency", choose MAF for which you would like to filter variants.  Annotated variants will be removed based on the ancestries chosen under the "Ancestry" section.

SubjectInfo.csv file
This file is optional.  If it is used, Individual IDs (first column) should correspond to the individual IDs in the vcd file.  The other 3 columns can be changed to adjust the parameters for which you want to consider in the counting (e.g. count samples by gender, phenotype, ancestry, etc.)

Other parameters:
The user may input a .bed file with chromosomal regions to narrow the queried variants in the vcf file.  

Output:

CountFile.txt (filename may be specified in the command) - the output is a text file with transcripts or genes as rows and sample IDs as columns.  For each transcript/gene, the total count of individuals harboring qualifying mutations is in the first column after the transcript/gene name, and the remaining columns show which sample harbors the mutations.

OutputFile.txt (filename may be specified in the command) - this output file lists the variants composing the counts in CountFile.txt, and which individuals have each variant.

Usage example:
python VarCount.py -i input.vcf -o OutputFile.txt -c CountFile.txt -p Parameter_file -b Optional_bed.bed -s Optional_Subject_Info.csv
