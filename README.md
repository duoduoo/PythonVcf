# Python_vcf
This folder contains python code related to manipulating vcf files.

* **`extract_vcf.py`:** Extract a region from an uncompressed vcf file.

```python
python extract_vcf.py [inputvcf] [chr_number] [start] [end] [output]
```

* **`cal_ld_1vcf.py`:** Calculate linkage disequilibrium for two sets of SNPs, using vcf file from 1000 Genome Phase 3.
SNP list id in "rsXXXXX" format.

```python
python cal_ld_1vcf.py [SNP_vcf] [SNP1_IDs] [SNP2_IDs] [output] [error]
```
