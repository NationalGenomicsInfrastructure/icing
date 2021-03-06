import click
from Bio import SeqIO


@click.command(context_settings = dict( help_option_names = ['-h', '--help'] ))
@click.option('--fasta','-f', type=str, help='multi FASTA input file to separate')
@click.option('--prefix','-p', type=str, help='prefix to add to the ID and filename',default='')
def separateFasta(fasta,prefix):
    """Utility to separate a multi-FASTA to separate FASTA files"""
    
    for seq_record in SeqIO.parse(fasta, "fasta"):
        # we are changing the IDs from "HLA:HLA..." to "HLA...". This makes IGV and other tools much happier
        # would be nice to know what lead to this idiocity in naming conventions
        seq_record.id = seq_record.id.replace("HLA:","")
        seq_record.id = prefix + seq_record.id
        print("writing " + seq_record.id)
        SeqIO.write(seq_record,seq_record.id + ".fasta", "fasta")

if __name__ == '__main__':
    separateFasta()
    
