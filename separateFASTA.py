import click
from Bio import SeqIO


@click.command(context_settings = dict( help_option_names = ['-h', '--help'] ))
@click.option('--fasta','-f', type=str, help='multi FASTA input file to separate')
def separateFasta(fasta):
    """Utility to separate a multi-FASTA to separate FASTA files"""
    
    for seq_record in SeqIO.parse(fasta, "fasta"):
        print("writing " + seq_record.id)
        SeqIO.write(seq_record,seq_record.id + ".fasta", "fasta")

if __name__ == '__main__':
    separateFasta()
    
