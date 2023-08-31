#!/usr/bin/env python
"""
Modified autoPVS1 wrapper script, as original did not work but overall scoring logic preserved.
Overall goal is to classify "Null Variants".
A good definition of what a PVS1 Null variant is: https://www.baylorgenetics.com/variant-classification/
"""
# -*- coding:utf-8 -*-
# author: Jiguang Peng
# created: 2019/6/27 19:13
# Modified: 2022/10/05
# Mod_author: Miguel Brown

import sys
import argparse
import pysam
from read_data_mod import transcripts_hg38
from pvs1 import PVS1
from utils import get_transcript, vep_consequence_trans, VCFRecord
__version__ = 'v2.0.0'



def tie_breaker(record, candidate_csqs, tx_candidates, pick, csq_fields, summary):
    """
    Takes in pysam record, candidate_csqs and tx_candidates from get_candidates(), csq_fields from pysam description, and summary dict defined in main()
    Tries to rank candidates, if still > 1, go with PICK if possible, transcript length as last resort
    Returns best ranked candidates and updated summary dict

    """
    # from https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences
    csq_severity_order = ["transcript_ablation","splice_acceptor_variant","splice_donor_variant","stop_gained","frameshift_variant","stop_lost","start_lost","transcript_amplification","feature_elongation","feature_truncation","inframe_insertion","inframe_deletion","missense_variant","protein_altering_variant","splice_donor_5th_base_variant","splice_region_variant","splice_donor_region_variant","splice_polypyrimidine_tract_variant","incomplete_terminal_codon_variant","start_retained_variant","stop_retained_variant","synonymous_variant","coding_sequence_variant","mature_miRNA_variant","5_prime_UTR_variant","3_prime_UTR_variant","non_coding_transcript_exon_variant","intron_variant","NMD_transcript_variant","non_coding_transcript_variant","coding_transcript_variant","upstream_gene_variant","downstream_gene_variant","TFBS_ablation","TFBS_amplification","TF_binding_site_variant","regulatory_region_ablation","regulatory_region_amplification","regulatory_region_variant","intergenic_variant","sequence_variant"]
    best_rank = 1000
    champ = 0
    csq_idx = csq_fields.index("Consequence")
    score_dict = {}
    # for c in range(len(candidate_csqs)):
    for candidate_index, consequence_string in enumerate(candidate_csqs):
        # Some variants have "compound csqs" separated by '&', subdivide and go with best first
        csq_values = consequence_string[csq_idx].split("&")
        for csq_value in csq_values:
            score = csq_severity_order.index(csq_value)
            if score <= best_rank:
                # even if found, keep going in case a later rank is even better
                best_rank = score
                champ = candidate_index
                if score not in score_dict:
                    score_dict[score] = []
                score_dict[score].append(candidate_index)
    if len(score_dict[best_rank]) > 1:
        print("Multiple candidates found for {} {}, looking for PICK and lengths".format(record.contig, record.pos), file=sys.stderr)
        # also track tx lengths
        longest = 0
        tx_len_dict = {}
        pick_break = 0
        for tied in score_dict[best_rank]:
            tx_len = tx_candidates[tied].tx_length
            if tx_len > longest:
                longest = tx_len
                tx_len_dict[longest] = tied
            if candidate_csqs[tied] == pick[1]:
                summary['pick'] +=1
                tie_break_return =  pick[0], pick[1], summary
                pick_break = 1
                break
        # If not found, warn that first rank hit used
        if not pick_break:
            print("WARN: No matching PICK, default to longest transcript", file=sys.stderr)
            summary['length'] +=1
            tie_break_return =  tx_candidates[tx_len_dict[longest]], candidate_csqs[tx_len_dict[longest]], summary
    else:
        summary['rank'] += 1
        tie_break_return = tx_candidates[champ], candidate_csqs[champ], summary
    return tie_break_return


def get_candidates(record_csq, csq_fields):
    """
    Takes in a pysam csq record and the csq_fields list parsed from the record description.
    Returns a list of picked csqs, matching transcript objects, and the PICK csq
    1. See if the tx candidate is canonical
    2. If canonical, check to see if it exists in the tool transcript reference
    3. Record the PICK value for rule 5 from pick_transcript
    """
    canon_idx = csq_fields.index("CANONICAL")
    pick_index = csq_fields.index("PICK")
    refseq_index = csq_fields.index('RefSeq')
    pick_value = "1"
    canon_value = "YES"
    
    candidate_csqs = []
    tx_candidates = []

    for csq in record_csq:
        info_list = csq.split("|")
        refseq_id = info_list[refseq_index]
        if info_list[pick_index] == pick_value:
            pick = (get_transcript(refseq_id, transcripts_hg38), info_list)
        if info_list[canon_idx] == canon_value and refseq_id != "":
            check = get_transcript(refseq_id, transcripts_hg38)
            if check is not None:
                candidate_csqs.append(info_list)
                tx_candidates.append(check)
    return candidate_csqs, tx_candidates, pick


def pick_transcript(record, csq_fields, summary):
    """
    Function to process a record and choose a representative transcript
    Does so by doing the following:
    1. See if the tx candidate is canonical
    2. If canonical, check to see if it exists in the tool transcript reference
    3. If multiple hits exists, use rank list to chose the worst
    4. If still multiple (meaning ranks are tied), use longest transcript length as defined by the reference
    5. If no hits from step 1 and 2, just revert to using PICK
    """
    candidate_csqs, tx_candidates, pick = get_candidates(record.info['CSQ'], csq_fields)
    if len(candidate_csqs) == 1:
        summary['canonical'] += 1
        return_values =  tx_candidates[0], candidate_csqs[0], summary

    # tie break with rank - if multiple share rank, go with PICK, if not go with tx length
    if len(candidate_csqs) > 1:
        return_values = tie_breaker(record, candidate_csqs, tx_candidates, pick, csq_fields, summary)
    # go with PICK
    if len(candidate_csqs) == 0:
        summary['pick'] +=1
        return_values =  pick[0], pick[1], summary
    return return_values


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vep_vcf', help='VCF classify')
    parser.add_argument('--genome_version', help='Currently only hg38 or GRCh38 supported')
    parser.add_argument('--version', action='version',
                        version='%(prog)s {version}'.format(version=__version__))
    args = parser.parse_args()
    lof_type = ['frameshift', 'nonsense', 'splice-5', 'splice-3', 'init-loss']
    genome_version = args.genome_version
    in_vcf = pysam.VariantFile(args.vep_vcf, threads=8)

    csq_fields = in_vcf.header.info['CSQ'].description.replace("Consequence annotations from Ensembl VEP. Format: ", "").split("|")
    print ("vcf_id",'SYMBOL','Feature','trans_name','consequence', 'strength_raw', 'strength','criterion', sep="\t")
    summary = { "canonical": 0, "pick": 0, "rank": 0, "length": 0 }
    for record in in_vcf.fetch():
        # Get representative transcript to score
        transcript, first_picked, summary = pick_transcript(record, csq_fields, summary)
        info = {}
        for i in range(len(first_picked)):
            info[csq_fields[i]] = first_picked[i]
            # adjust this field to match behavior of original auto pvs1 script
            if csq_fields[i] == 'HGVSp':
                info[csq_fields[i]] = first_picked[i].replace('%3D', '=')

        vcfrecord = VCFRecord(record.contig.replace('chr', ''), str(record.pos), record.ref, record.alts[0])

        consequence = vep_consequence_trans(info['Consequence'])
        vcf_id = "-".join([vcfrecord.chrom, str(vcfrecord.pos), vcfrecord.ref, vcfrecord.alt])
        if consequence in lof_type and transcript:
            lof_pvs1 = PVS1(vcfrecord, consequence, info['HGVSc'], info['HGVSp'], transcript, genome_version)
            trans_name = lof_pvs1.transcript.full_name
            print(vcf_id,
                  info['SYMBOL'],
                  info['Feature'],
                  trans_name,
                  lof_pvs1.consequence,
                  lof_pvs1.strength_raw.name,
                  lof_pvs1.strength.name,
                  lof_pvs1.criterion,
                  sep="\t")
        elif consequence in lof_type:
            print(vcf_id,
                  info['SYMBOL'],
                  info['Feature'],
                  'not_canonical',
                  consequence,
                  'Unmet',
                  'Unmet',
                  'na',
                  sep="\t")
        else:
            print(vcf_id,
                  info['SYMBOL'],
                  info['Feature'],
                  'not_lof',
                  consequence,
                  'Unmet',
                  'Unmet',
                  'na',
                  sep="\t")
    print("Summary of transcript pick categories: Canonical: {}, Rank: {}, Pick: {}, Length: {}".format(summary['canonical'], summary['rank'], summary['pick'], summary['length']), file=sys.stderr)

if __name__ == '__main__':
    main()
