import pysam
import sys
import argparse
import pandas as pd
from collections import defaultdict

def split_bam_by_cells(input_bam, cell_table_file, output_prefix=""):
    """
    Split BAM file into multiple output BAMs based on cell assignments.
    
    Args:
        input_bam: Path to input BAM file
        cell_table_file: Path to TSV/CSV file with two columns: cell_barcode, output_bam_name
        output_prefix: Optional prefix for output BAM files
    """
    # Read cell assignment table
    cell_table = pd.read_csv(cell_table_file, sep='\t', header=None)
    cell_table.columns = ['cell_barcode', 'output_bam']
    cell_table['output_bam'] = cell_table['output_bam'] + '.bam'
    print(cell_table.head())
    
    # Create dictionary mapping cell barcodes to output BAM names
    cell_to_bam = dict(zip(cell_table['cell_barcode'], cell_table['output_bam']))
    
    # Get unique output BAM names
    output_bams = cell_table['output_bam'].unique()
    
    print(f"Loaded {len(cell_to_bam)} cell barcodes")
    print(f"Splitting into {len(output_bams)} output BAM files: {list(output_bams)}")
    
    # Open input BAM and create output BAM file handles
    with pysam.AlignmentFile(input_bam, "rb") as inbam:
        # Create dictionary of output BAM file handles
        outbam_handles = {}
        for bam_name in output_bams:
            output_path = f"{output_prefix}{bam_name}" if output_prefix else bam_name
            outbam_handles[bam_name] = pysam.AlignmentFile(output_path, "wb", header=inbam.header)
        
        try:
            total_reads = 0
            numKeptReads=0
            kept_reads = defaultdict(int)
            
            for read in inbam:
                total_reads += 1
                
                # Get CB tag if it exists
                if read.has_tag('CB'):
                    cb = read.get_tag('CB')
                    
                    # Write read to appropriate output BAM if cell barcode is in our table
                    if cb in cell_to_bam:
                        output_bam_name = cell_to_bam[cb]
                        outbam_handles[output_bam_name].write(read)
                        kept_reads[output_bam_name] += 1
                        numKeptReads += 1

                # Progress update every 1,000,000 reads
                if total_reads % 1000000 == 0:
                    print(f"Processed {total_reads} reads")
                    print(f"\nTotal kept reads: {numKeptReads} out of {total_reads} total reads")
        
        finally:
            # Close all output BAM files
            for handle in outbam_handles.values():
                handle.close()
    
    print(f"\nDone! Processed {total_reads} total reads")
    for bam_name in output_bams:
        print(f"  {bam_name}: {kept_reads[bam_name]} reads")

def main():
    parser = argparse.ArgumentParser(
        description='Split BAM file into multiple BAMs based on cell barcode assignments'
    )
    parser.add_argument('input_bam', help='Input BAM file')
    parser.add_argument('cell_table', help='Table file (TSV/CSV) with columns: cell_barcode, output_bam_name')
    parser.add_argument('--prefix', default='', help='Optional prefix for output BAM files')
    
    args = parser.parse_args()
    
    split_bam_by_cells(args.input_bam, args.cell_table, args.prefix)

if __name__ == '__main__':
    main()
