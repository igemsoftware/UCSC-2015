#Test files for FOCUS using a beta-glucosidase from Haloarcula Hispanica

```beta-gluc_gene.txt``` is the FASTA nucleotide sequence for a beta-glucosidase from H. hispanica
```Hispanica_CodonBias.txt``` is the GCG codon bias table. 
```SS_beta-glucosidase.txt``` is the secondary structure prediction from PSSpred in I-TASSER format. This file is optional when running ```FOCUS.py```

To use these test files, they must be in the same directory as ```FOCUS.py``` (See user guide).

# Usage
```python3 FOCUS.py Hispanica_CodonBias.txt beta-gluc_gene.txt SS_beta-glucosidase.txt```
