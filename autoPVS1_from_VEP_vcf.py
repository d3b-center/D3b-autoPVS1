#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: Jiguang Peng
# created: 2019/6/27 19:13
# Modified: 2022/10/05
# Mod_author: Miguel Brown

import sys
import argparse
from collections import namedtuple
import pysam
from read_data_mod import transcripts_hg38
from pvs1 import PVS1
from utils import get_transcript, vep_consequence_trans, VCFRecord
__version__ = 'v0.3.0'


def pick_transcript(record, csq_fields, summary, csq_severity_order):
    # Rules - is canonical, is in transcript reference
    canon_idx = csq_fields.index("CANONICAL")
    p_idx = csq_fields.index("PICK")
    csq_idx = csq_fields.index("Consequence")
    pick_value = "1"
    canon_value = "YES"
    r_idx = csq_fields.index('RefSeq')
    picked_csqs = []
    tx_candidates = []

    for csq in record.info['CSQ']:
        info_list = csq.split("|")
        refseq_id = info_list[r_idx]
        if info_list[p_idx] == pick_value:
            pick = (get_transcript(refseq_id, transcripts_hg38), info_list)
        if info_list[canon_idx] == canon_value:
            if refseq_id != "":
                # pdb.set_trace()
                check = get_transcript(refseq_id, transcripts_hg38)
                if check is not None:
                    picked_csqs.append(info_list)
                    tx_candidates.append(check)
    # tie break with rank - if multiple share rank, go first best found
    if len(picked_csqs) > 1:
        best_rank = 1000
        champ = 0
        summary['rank'] += 1
        for c in range(len(picked_csqs)):
            # Some variants have "compound csqs" seprated by '&', subdivide and go with worst first
            csq_values = picked_csqs[c][csq_idx].split("&")
            for csq_value in csq_values:
                score = csq_severity_order.index(csq_value)
                if score < best_rank:
                    # even if found, keep going in case a later rank is even worse
                    best_rank = score
                    champ = c
        return tx_candidates[champ], picked_csqs[champ], summary
    # go with PICK
    if len(picked_csqs) == 0:
        summary['pick'] +=1
        return pick[0], pick[1], summary
    else:
        summary['canonical'] += 1
        return tx_candidates[0], picked_csqs[0], summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vep_vcf', help='VCF classify')
    parser.add_argument('--genome_version', help='Currently only hg38 or GRCh38 supported')
    parser.add_argument('--version', action='version',
                        version='%(prog)s {version}'.format(version=__version__))
    args = parser.parse_args()
    lof_type = ['frameshift', 'nonsense', 'splice-5', 'splice-3', 'init-loss']
    csq_severity_order = ["transcript_ablation","splice_acceptor_variant","splice_donor_variant","stop_gained","frameshift_variant","stop_lost","start_lost","transcript_amplification","feature_elongation","feature_truncation","inframe_insertion","inframe_deletion","missense_variant","protein_altering_variant","splice_donor_5th_base_variant","splice_region_variant","splice_donor_region_variant","splice_polypyrimidine_tract_variant","incomplete_terminal_codon_variant","start_retained_variant","stop_retained_variant","synonymous_variant","coding_sequence_variant","mature_miRNA_variant","5_prime_UTR_variant","3_prime_UTR_variant","non_coding_transcript_exon_variant","intron_variant","NMD_transcript_variant","non_coding_transcript_variant","coding_transcript_variant","upstream_gene_variant","downstream_gene_variant","TFBS_ablation","TFBS_amplification","TF_binding_site_variant","regulatory_region_ablation","regulatory_region_amplification","regulatory_region_variant","intergenic_variant","sequence_variant"]

    genome_version = args.genome_version
    in_vcf = pysam.VariantFile(args.vep_vcf, threads=8)

    # Use VEP PICK field to choose a representative transcript
    csq_fields = in_vcf.header.info['CSQ'].description.replace("Consequence annotations from Ensembl VEP. Format: ", "").split("|")

    print ("vcf_id",'SYMBOL','Feature','trans_name','consequence', 'strength_raw', 'strength','criterion', sep="\t")
    summary = { "canonical": 0, "pick": 0, "rank": 0 }
    for record in in_vcf.fetch():
        # Parse CSQ for PICKed transcript
        # Populate dict with CSQ key-value pairs
        
        transcript, first_picked, summary = pick_transcript(record, csq_fields, summary, csq_severity_order)
        info = {}
        for i in range(len(first_picked)):
            info[csq_fields[i]] = first_picked[i]
            # adjust this filed to match behavior of original auto pvs1 script
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
    print("Summary of transcript pick categories: Canonical: {}, Rank: {}, Pick: {}".format(summary['canonical'], summary['rank'], summary['pick']), file=sys.stderr)

if __name__ == '__main__':
    main()
