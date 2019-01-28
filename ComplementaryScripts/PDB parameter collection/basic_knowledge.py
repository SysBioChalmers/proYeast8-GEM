
"""
genotype
0/0 - the sample is homozygous reference
0/1 - the sample is heterozygous, carrying 1 copy of each of the REF and ALT alleles
1/1 - the sample is homozygous alternate
"""


'''
In order to extract homo or hetero variants, there are different tools like: GATK SelectVariants, SnpSift
'''

'''
VCFtools and VCFlib support filtering by genotype.

'''


'''
I thought it was standard that "./." meant there was no confident genotype to be called, and "0/0" in that case that homozygous ref was called. The best thing to do is just clarify this with your collaborators
'''


'''
AFAIK, "./." means the call at that position is missing. It could be missing for a variety of reasons. It might be missing because a call was made, but it didn't meet some threshold and was filtered out. It could also be missing because multiple VCFs are merged together, and samples that don't have a call at a position are listed as "./." Usually by default, when someone runs a variant caller on an individual sample, only the variant calls are emitted so there are no "0/0" reference calls in the VCF file, which is fine if you want to know the variants of a single sample. This becomes a problem when you want to compare SNP calls across samples, because you can't assume that an absent call in the VCF means it was a reference call (because it could have also been a position where the caller couldn't make an accurate call). One option is to force the caller to emit all calls, even reference calls. This will generate very large files. Another option is to call SNPs simultaneously on the samples and output a multi-sample VCF (which can be done with samtools mpileup or GATK). For this solution, if one of the samples has a variant, then the calls for all the other samples will be emitted too, even if it's reference or low quality. I prefer this way, because the VCF file is still small, and it is then an easy task to find positions that have different calls across the samples.

Hope that is of some help. We've wrestled with how best to handle this quite a bit, so definitely interesting to hear how others manage it.

Justin
'''


'''
If you read the VCF format specification pdf (http://samtools.github.io/hts-specs/VCFv4.1.pdf) you'll have the answer. Summarizing, GT represent the genotype, encoded as allele values separated by / or |. If the allele value is 0, means that it is equal to the reference allele (what is in REF field), if 1 mean that is equal to alternative (first allele listed in ALT), and if 2 it is equal to the second allele listed in ALT (if it exists).

So a SNP tagged with GT = 1/1 represent a SNP homozygous for the ALT allele (1/0 heterozygous, 0/0 homozygous for the reference).

For determining if the call is homozygous or heterozygous I suggest you to read this thread: How To Distinguish Heterozygotes And Homozygotes From Variants In Vcf Format?. It can be done according to different criteria.
'''


'''
GATK
If for any of these reasons you find that you cannot perform variant recalibration on your data (after having tried the workarounds that we recommend, where applicable), you will need to use hard-filtering instead. This consists of setting flat thresholds for specific annotations and applying them to all variants equally. See the methods articles and FAQs for more details on how to do this.
The 1/2 genotype simply means the genotype contains the first alternate allele and the second alternate allele. The reference allele is 0, the first alternate allele is 1, the second alternate allele is 2...and so on.
'''