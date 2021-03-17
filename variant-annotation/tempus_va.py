import argparse
import logging
import statistics
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

    arguments:
    variant_ids: A list of variant IDs, in the form of CHR-POS-REF-ALT, to be passed to the ExAC service

    return:
    JSON object from API call, as a Python nested dictionary
    '''
    res = requests.post(BULK_VARIANT_URL, json=variant_ids)
    res.raise_for_status()
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
        mean_missense_score = statistics.mean(missense_scores) if len(missense_scores) != 0 else 1
        return (consequence_score[0], mean_missense_score, consequence_score[1])

    if 'consequence' not in record or record['consequence'] is None:
        return None
    effects = [c3 for c1 in record['consequence'].values() for c2 in c1.values() for c3 in c2]
    if len(effects) == 0:
        return None
    return min(effects, key=key)


def variant_annotate(vcf_file: str, cons_file: str, missing: str = '.') -> pd.DataFrame:
    '''
    Performs variant annotation.

    arguments:
    vcf_file: The filename of a VCF file containing the variants to-be-annotated
    conf_file: The filename of TSV file ranking potential effects to be encountered in ExAC
    missing: To be used in place of undetermined values associated with unannotated variants in the ExAC service

    returns:
    A Pandas DataFrame containing the annotated variants
    '''
    variant_ids = []
    vcf_rows = []
    for var in VCF(vcf_file):
        # If there is more than one non-reference allele, there will be
        # separate read information that must be extracted for each allele
        row = {
            'chr': var.CHROM,
            'pos': var.POS,
            'ref': var.REF,
            'coverage_depth': var.INFO.get('DP'),
        }

        types = var.INFO.get('TYPE').split(',')
        if var.INFO.get('NUMALT') == 1:
            variant_reads = [var.INFO.get('AO')]
        else:
            variant_reads = var.INFO.get('AO')

        for i in range(var.INFO.get('NUMALT')):
            variant_ids.append('%s-%d-%s-%s' % (var.CHROM, var.POS, var.REF, var.ALT[i]))
            row.update({
                'alt': var.ALT[i],
                'type': types[i],
                'variant_reads': variant_reads[i],
            })
            vcf_rows.append(row.copy())
    
    consequences = parse_consequences(cons_file)
    data = ExAC_POST_bulk_variants(variant_ids)

    variant_rows = []
    for var in variant_ids:
        record = data[var]
        effect = most_deleterious_effect(record, consequences)
        if 'variant' not in record or 'allele_freq' not in record['variant']:
            allele_freq = missing
        else:
            allele_freq = record['variant']['allele_freq']
        variant_rows.append({
            'allele_frequency': allele_freq,
            'effect': effect['major_consequence'] if effect else missing,
            'SIFT': effect['SIFT'] if effect else missing,
            'PolyPhen': effect['PolyPhen'] if effect else missing,
        })
    
    df = pd.concat([
        pd.DataFrame.from_records(vcf_rows),
        pd.DataFrame.from_records(variant_rows)
    ], axis='columns')
    df['variant_percent'] = 100*df['variant_reads'] / df['coverage_depth']
    return df


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('input', metavar='VCF', type=str, help='variant call filename in vcf format')
    parser.add_argument('-c', metavar='CONSEQUENCES', type=str, default='consequences.tsv', help='a TSV containing information to rank variant consequences')
    parser.add_argument('--out', metavar='TSV', required=False, type=str, help='output filename in tsv format')
    args = parser.parse_args()

    annotated = variant_annotate(args.input, args.c)
    if args.out:
        annotated.to_csv(path_or_buf=args.out, sep='\t', index=False)
    else:
        print(annotated.to_csv(sep='\t', index=False))


if __name__ == '__main__':
    main()
