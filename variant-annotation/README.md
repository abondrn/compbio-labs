
Link to repo: https://github.com/abondrn/compbio-labs/tree/master/variant-annotation

# Variant Annotation

Prototype of a variant annotation tool. Data provided as a part of a Tempus Labs challenge. Effect and frequency lookups done using the Broad Institute [ExAC API](http://exac.hms.harvard.edu/).

## Setup

1. Install the required packages, for example by running `pip install -r requirements.txt`. Note that this requires Python 3.5 or newer as I utilize type annotations.

2. Run the script, for example `python3 tempus_va.py Challenge_data.vcf -o Challenge_annotated.tsv` was used to generate the data in this repo.

```
usage: tempus_va.py [-h] [-c CONSEQUENCES] [-o TSV] VCF

Prototype of a variant annotation tool

positional arguments:
  VCF              variant call filename in vcf format

optional arguments:
  -h, --help       show this help message and exit
  -c CONSEQUENCES  a TSV containing information to rank variant consequences
  -o TSV           output filename in TSV format
```

## Output Format

Each variant is annotated with the following pieces of information:

1. Locus of the variant: chromosome `chr`, position `pos`, reference allele `ref`, alterate allele `alt`.

1. Type of variation, which can be encoded as:

   - [SNPs](https://en.wikipedia.org/wiki/Single-nucleotide_polymorphism) are point mutations in which one base is substituted for another, which are broken down into transition (Purine⇒Purine or Pyrimidine⇒Pyrimidine, notated as `ts`) and transversion (Purine⇒Pyrimidine and vice versa, `tv`) mutations
   - Analogously, there are also `mnp` (Multiple Nucleotide Polymorphism) mutations for substitutions than span multiple bases
   - Indels collectively refer to insertions (`ins`) and deletions (`del`) in which the variant results in bases added or removed
   - The remainder are grouped under `complex`. A better overview is given [here](http://www2.csudh.edu/nsturm/CHEMXL153/DNAMutationRepair.htm).

1. `coverage_depth`: Depth of sequence coverage at the site of variation.

1. `variant_reads`: Number of reads supporting the variant.

1. `effect`: The most delererious effect (missense, silent, intergenic, etc.). A list of ranked effects can be found [here](https://github.com/abondrn/compbio-labs/blob/master/variant-annotation/consequences.tsv). For variants annotated as `missense_variant`, `SIFT` and `PolyPhen` scores may be available which predict pathogenic loss-of-function in protein-coding regions.

1. `allele_frequency`: Allele frequency of variant in the population.

1. Percentage of reads supporting the variant (`variant_percent`) versus those supporting reference reads (`reference_percent`).

## Details

 - Some variants in the VCF input have multiple alternate alleles. In the corresponding output table, this variant is split up into multiple annotated alleles so that there is one row per allele.
   - In the case where one of these alleles is a SNP (or an indel), there is superfluous context appended for confirmity with the other alternate alleles. In this case, normalization (detailed below) becomes especially important.
   - When calculating percentage of reads supporting the variant, reads supporting other alternate variants confound the interpretation somewhat.
 - There are multiple ways to represent the same allele, which can present problems at query time. To address this, I follow the normalization algorithm detailed [here](https://genome.sph.umich.edu/wiki/Variant_Normalization) and use the normalized representation of the allele to query ExAC, resulting in greatly increased chances of a hit.

#### Effect Ranking

When there is a variant hit in ExAC, it is common to have multiple effect matches; how then do you choose the one with the most significant consequence?

 - Variant effects in ExAC come with Sequence Ontology consequence terms, which describe how the change disrupts the gene expression process. Reading the documentation for the [Ensembl genome browser](https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html) and  [GEMINI framework schema](https://gemini.readthedocs.io/en/latest/content/database_schema.html#details-of-the-impact-and-impact-severity-columns), I arrived at a grouping and ordering of the consequences by impact, which I compiled into [consequences.tsv](https://github.com/abondrn/compbio-labs/blob/master/variant-annotation/consequences.tsv). This has the ability being easily updated with new terms, and the client is even able to supply alternate consequences with the `-c` option.
 - Based on the Ensembl IMPACT levels, the ratings able to be given to consequences are (in order of severity): `HIGH, MODERATE, MODIFIER, LOW`. While it is hard to pin down the severity of a `MODIFIER` effect, it can present a larger tail risk and may be more worth investigating than one of known `LOW` impact, and so is given priority when ranking. Unknown variants are given `LOW` impact.
 - In addition to impact score, the average predicted loss-of-function scores and relative position in the table to perform multikey sort.

## Future Work

 - Rather than one bulk request, it might be better to batch the variant lookups by asynchronously fetching them in groups. Since currently processing times are long, this might provide a more consistent running time, and allow better estimation of time left which could come in handy for large workloads.

 - The selection of the most deleterious effect could be improved. See [REVEL](https://illumina.github.io/NirvanaDocumentation/3.14/data-sources/revel) for how to better evaluate missense variants by ensembling scores by several systems, but that is just one approach.

 - It might be useful to think about how the client plans to ingest the data. The GEMINI project is interesting in this respect as it not only annotates the inputted variants, but then allows you to manipulate the output data using SQL, thus enabling downstream analysis.

 - Like all variant annotation tools, the significance of each variant is evaluated independently. It might be useful to combine the read counts of each variant, as well as their individual pathogenicity, to come to a holistic prediction, such as the presence or absense of a particular phenotype.