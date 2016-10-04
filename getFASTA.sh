#!/bin/bash
wget http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=imgthla\&id=${1}\&format=fasta\&style=raw\&Retrieve=Retrieve -O ${1}.fasta
