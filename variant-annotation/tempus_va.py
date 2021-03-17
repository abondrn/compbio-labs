
#!/usr/bin/env python3

import argparse
import logging
from statistics import mean
from math import nan
from typing import Any, List, Dict, Tuple
from collections import defaultdict
from collections.abc import Sequence

import pandas as pd
import requests
from cyvcf2 import VCF


BULK_VARIANT_URL = 'http://exac.hms.harvard.edu/rest/bulk/variant'


def ExAC_POST_bulk_variants(variant_ids: List[str]) -> Dict[str, Any]:
    '''
    Given a list of variant IDs, retrieve the variant information from the Broad ExAC API.

    Alleles are specified by the ExAC API as a  string.

    arguments:
    variant_ids: A list of variant IDs, in the form of
    CHROMOSOME-POSITION-REFERENCE-VARIANT, to be passed to the ExAC service

    return:
    JSON object from API call, as a Python nested dictionary
    '''
    logging.info('Fetching bulk variants')
    res = requests.post(BULK_VARIANT_URL, json=variant_ids)
    res.raise_for_status()
    logging.info('Done fetching')
    return res.json()


IMPACTS = {
    'HIGH': 0,
    'MODERATE': 1,
    'MODIFIER': 2,
    'LOW': 3,
}


def parse_consequences(tsv_file: str) -> Dict[str, Tuple[int, int]]:
    '''
    Parses a file containing consequence information into a dictionary used to rank the severity of the effects.

    arguments:
    tsv_file: A filename containing consequence information in the form of a tab-seperated values

    returns:
    A dictionary of consequences indexed by the SO term, and containing the enumeration of the impact category (whether HIGH, MODERATE, etc) and relative position (for resolving ties)
    '''
    df = pd.read_csv(tsv_file, sep='\t')
    key = defaultdict(lambda c: 1)
    for i, r in df.iterrows():
        key[r['SO term']] = (IMPACTS[r['IMPACT']], i)
    return key


def most_deleterious_effect(record: Dict[str, Any], consequences: Dict[str, int]) -> Dict[str, Any]:
    '''
    Uses information associated with a variant from the ExAC service to elevate the effect determined to be the most deleterious

    arguments:
    record: An ExAC variant record as serialized JSON
    consequences: A consequence ranking dictionary, such as one that might be returned from `parse_consequences`

    returns:
    The JSON associated with the effect deemed to be the most clinically significant
    '''
    def key(effect):
        consequence_score = consequences[effect['major_consequence']]
        missense_scores = []
        for program in ['SIFT', 'PolyPhen']:
            if effect[program] != '':
                between_parens = effect[program].partition('(')[2].rstrip(')')
                score = float(between_parens)
                if program == 'PolyPhen':
                    score = 1-score
                missense_scores.append(score)
        mean_missense_score = mean(missense_scores) if len(missense_scores) != 0 else 1
        return (consequence_score[0], mean_missense_score, consequence_score[1])

    if 'consequence' not in record or record['consequence'] is None:
        return None
    effects = [c3 for c1 in record['consequence'].values() for c2 in c1.values() for c3 in c2]
    if len(effects) == 0:
        return None
    return min(effects, key=key)


def normalize_variant(pos: int, ref: str, alt: str) -> Tuple[int, str, str]:
    '''Perform variant normalization'''
    limit = min(len(ref), len(alt))
    # right trim
    k = 0
    while limit+k > 1 and ref[k-1] == alt[k-1]:
        k -= 1
    # left trim
    j = 0
    while j+1 < limit+k and ref[j] == alt[j]:
        j += 1
    
    return pos+j, ref[j:len(ref)+k], alt[j:len(alt)+k]


PURINE = 'AG'
PYRIMIDINE = 'CT'


def tabulate_vcf(vcf_file: str) -> pd.DataFrame:
    vcf_rows = []
    for var in VCF(vcf_file):
        # If there is more than one non-reference allele, there will be
        # separate read information that must be extracted for each allele
        types = var.INFO.get('TYPE').split(',')
        if var.INFO.get('NUMALT') == 1:
            variant_reads = [var.INFO.get('AO')]
        else:
            variant_reads = var.INFO.get('AO')

        for i in range(var.INFO.get('NUMALT')):
            pos, ref, alt = normalize_variant(var.POS, var.REF, var.ALT[i])

            type = types[i]
            if type == 'snp':
                try:
                    assert len(ref) == len(alt) == 1
                except:
                    print(i, j, k, var)
                    raise
                # SNP was grouped with ins or del, resulting in superfluous read lengths
                ref, alt = ref[0], alt[0]
                if (ref in PURINE and alt in PURINE) or (ref in PYRIMIDINE and alt in PYRIMIDINE):
                    type = 'ts'
                else:
                    type = 'tv'
        
            variant_id = '%s-%d-%s-%s' % (var.CHROM, pos, ref, alt)
            vcf_rows.append({
                'id': variant_id,
                'chr': var.CHROM,
                'pos': var.POS,
                'type': type,
                'ref': var.REF,
                'alt': var.ALT[i],
                'coverage_depth': var.INFO.get('DP'),
                'reference_reads': var.INFO.get('RO'),
                'variant_reads': variant_reads[i],
            })
    return pd.DataFrame.from_records(vcf_rows)


def tabulate_ExAC_variants(variant_ids: List[str], cons_file: str, missing: str) -> pd.DataFrame:
    consequences = parse_consequences(cons_file)
    variant_data = ExAC_POST_bulk_variants(variant_ids)

    variant_rows = []
    for variant_id in variant_ids:
        record = variant_data[variant_id]
        effect = most_deleterious_effect(record, consequences)
        if 'variant' not in record or 'allele_freq' not in record['variant']:
            allele_freq = nan
        else:
            allele_freq = float(record['variant']['allele_freq'])
        variant_rows.append({
            'effect': effect['major_consequence'] if effect else missing,
            'SIFT': effect['SIFT'] or '.' if effect else missing,
            'PolyPhen': effect['PolyPhen'] or '.' if effect else missing,
            'allele_frequency': allele_freq,
        })
    return pd.DataFrame.from_records(variant_rows)


def variant_annotate(vcf_file: str, cons_file: str, missing: str = '.') -> pd.DataFrame:
    '''
    Generates a table of annotated variants drawn from a VCF file.

    arguments:
    vcf_file: The filename of a VCF file containing the variants to-be-annotated
    conf_file: The filename of TSV file ranking potential effects to be encountered in ExAC
    missing: To be used in place of undetermined values associated with unannotated variants in the ExAC service

    returns:
    A Pandas DataFrame containing the annotated variants
    '''

    vcf_df = tabulate_vcf(vcf_file)
    vcf_df['variant_percent'] = 100*vcf_df['variant_reads'] / vcf_df['coverage_depth']
    vcf_df['reference_percent'] = 100*vcf_df['reference_reads'] / vcf_df['coverage_depth']

    exac_df = tabulate_ExAC_variants(list(vcf_df.id), cons_file, missing)
    
    df = pd.concat([
        vcf_df,
        exac_df
    ], axis='columns')
    return df


def main():
    '''Outputs variant annotation table based on the supplied command-line-arguments'''
    parser = argparse.ArgumentParser(description='Prototype of a variant annotation tool')

    parser.add_argument('input', metavar='VCF', type=str, help='variant call filename in VCF format')
    parser.add_argument('-c', metavar='CONSEQUENCES', type=str, default='consequences.tsv', help='a TSV containing information to rank variant consequences')
    parser.add_argument('-o', metavar='TSV', required=False, type=str, help='output filename in TSV format')
    parser.add_argument('-p', metavar='PRECISION', required=False, type=int, help='round figures to that many digits past the decimal point')
    args = parser.parse_args()

    annotated = variant_annotate(args.input, args.c)
    kwargs = dict(sep='\t', index=False)
    if args.p:
        assert args.p >= 0, 'Precision must be nonnegative'
        kwargs['float_format'] = f'%.{args.p}f'
    if args.o:
        annotated.to_csv(path_or_buf=args.o, **kwargs)
    else:
        print(annotated.to_csv(**kwargs))


if __name__ == '__main__':
    main()
