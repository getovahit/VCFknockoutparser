import vcf
import subprocess
import json
import os
from typing import List, Dict, Tuple
import unittest

def split_vcf_by_chromosome(vcf_file: str) -> Dict[str, str]:
    """
    Split the input VCF file into separate files by chromosome.
    """
    chromosome_files = {}
    reader = vcf.Reader(filename=vcf_file)
    
    for record in reader:
        if record.CHROM not in chromosome_files:
            chrom_file = f"{record.CHROM}_temp.vcf"
            chromosome_files[record.CHROM] = chrom_file
            writer = vcf.Writer(open(chrom_file, 'w'), reader)
        else:
            writer = vcf.Writer(open(chromosome_files[record.CHROM], 'a'), reader)
        writer.write_record(record)
    
    return chromosome_files

def run_vep(vcf_file: str) -> List[Dict]:
    """
    Run VEP on the input VCF file and return the results for all alleles.
    """
    vep_command = [
        "vep",
        "--input_file", vcf_file,
        "--format", "vcf",
        "--output_file", "STDOUT",
        "--json",
        "--cache",
        "--everything",
        "--allele_number"
    ]
    
    try:
        result = subprocess.run(vep_command, capture_output=True, text=True, check=True)
        vep_output = result.stdout.strip().split("\n")
        return [json.loads(line) for line in vep_output if line.strip()]
    except subprocess.CalledProcessError as e:
        print(f"Error running VEP: {e}")
        return []

def is_knockout(variant_data: Dict) -> bool:
    """
    Determine if a variant is likely to cause a gene knockout based on VEP annotations.
    """
    knockout_consequences = {
        'stop_gained',
        'frameshift_variant',
        'splice_donor_variant',
        'splice_acceptor_variant',
        'start_lost',
        'transcript_ablation'
    }
    
    for transcript_consequence in variant_data.get('transcript_consequences', []):
        if set(transcript_consequence.get('consequence_terms', [])) & knockout_consequences:
            return True
        
        if 'deletion' in transcript_consequence.get('consequence_terms', []) and transcript_consequence.get('cds_change_length', 0) > 100:
            return True
        
        if 'splice_region_variant' in transcript_consequence.get('consequence_terms', []):
            exon = transcript_consequence.get('exon', '').split('/')
            if len(exon) == 2 and (exon[0] == '1' or exon[0] == exon[1]):
                return True
    
    return False

def parse_vep_results(vep_results: List[Dict], vcf_file: str) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """
    Parse VEP results and identify potential gene knockouts for all alleles.
    """
    all_knockouts = {}
    homozygous_knockouts = {}

    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    
    for variant_data in vep_results:
        if is_knockout(variant_data):
            variant_id = variant_data['input']
            for transcript_consequence in variant_data.get('transcript_consequences', []):
                gene_symbol = transcript_consequence.get('gene_symbol', 'Unknown')
                if gene_symbol != 'Unknown':
                    if gene_symbol not in all_knockouts:
                        all_knockouts[gene_symbol] = []
                    all_knockouts[gene_symbol].append(variant_id)
                    
                    vcf_record = next(record for record in vcf_reader if f"{record.CHROM}:{record.POS}" in variant_id)
                    if vcf_record.num_het == 0 and vcf_record.num_hom_alt > 0:
                        if gene_symbol not in homozygous_knockouts:
                            homozygous_knockouts[gene_symbol] = []
                        homozygous_knockouts[gene_symbol].append(variant_id)
    
    return all_knockouts, homozygous_knockouts

def process_chromosome(chrom_file: str) -> Tuple[Dict[str, List[str]], Dict[str, List[str]]]:
    """
    Process a single chromosome file and return knockout results.
    """
    vep_results = run_vep(chrom_file)
    return parse_vep_results(vep_results, chrom_file)

def main(vcf_file: str):
    chromosome_files = split_vcf_by_chromosome(vcf_file)
    
    all_knockouts = {}
    homozygous_knockouts = {}
    
    for chrom, chrom_file in chromosome_files.items():
        print(f"Processing chromosome {chrom}...")
        chrom_all_knockouts, chrom_homozygous_knockouts = process_chromosome(chrom_file)
        
        all_knockouts.update(chrom_all_knockouts)
        homozygous_knockouts.update(chrom_homozygous_knockouts)
        
        os.remove(chrom_file)  # Clean up temporary files
    
    print(f"Total potential gene knockouts found: {len(all_knockouts)}")
    for gene, variants in sorted(all_knockouts.items()):
        print(f"{gene}: {', '.join(variants)}")
    
    print(f"\nHomozygous gene knockouts found: {len(homozygous_knockouts)}")
    for gene, variants in sorted(homozygous_knockouts.items()):
        print(f"{gene}: {', '.join(variants)}")

class TestKnockoutFunctions(unittest.TestCase):
    def test_is_knockout(self):
        # Test case for a knockout variant
        knockout_variant = {
            'transcript_consequences': [
                {'consequence_terms': ['stop_gained']}
            ]
        }
        self.assertTrue(is_knockout(knockout_variant))
        
        # Test case for a non-knockout variant
        non_knockout_variant = {
            'transcript_consequences': [
                {'consequence_terms': ['synonymous_variant']}
            ]
        }
        self.assertFalse(is_knockout(non_knockout_variant))
    
    def test_parse_vep_results(self):
        # Mock VEP results and VCF file for testing
        mock_vep_results = [
            {
                'input': 'chr1:1000A>T',
                'transcript_consequences': [
                    {'consequence_terms': ['stop_gained'], 'gene_symbol': 'GENE1'}
                ]
            },
            {
                'input': 'chr1:2000G>C',
                'transcript_consequences': [
                    {'consequence_terms': ['synonymous_variant'], 'gene_symbol': 'GENE2'}
                ]
            }
        ]
        
        # Create a temporary VCF file for testing
        temp_vcf = 'temp_test.vcf'
        with open(temp_vcf, 'w') as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
            f.write("chr1\t1000\t.\tA\tT\t.\tPASS\t.\tGT\t1/1\n")
            f.write("chr1\t2000\t.\tG\tC\t.\tPASS\t.\tGT\t0/1\n")
        
        all_knockouts, homozygous_knockouts = parse_vep_results(mock_vep_results, temp_vcf)
        
        self.assertEqual(len(all_knockouts), 1)
        self.assertIn('GENE1', all_knockouts)
        self.assertEqual(len(homozygous_knockouts), 1)
        self.assertIn('GENE1', homozygous_knockouts)
        
        os.remove(temp_vcf)  # Clean up

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        unittest.main(argv=['first-arg-is-ignored'], exit=False)
    else:
        if len(sys.argv) != 2:
            print("Usage: python script.py <vcf_file>")
            sys.exit(1)
        main(sys.argv[1])
