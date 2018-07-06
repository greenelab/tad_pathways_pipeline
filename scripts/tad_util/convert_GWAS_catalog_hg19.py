"""
(C) 2016 Gregory Way
convert_GWAS_catalog_hg19.py

Description:
The GWAS catalog is downloaded in hg38 coordinates. Use pyliftover to convert
coordinates to hg19.

Usage:
Called by ANALYSIS.sh

Output:
The GWAS Catalog but with hg19 coordinates
"""

from pyliftover import LiftOver
lo = LiftOver('data/hg38ToHg19.over.chain.gz')

new_coordinates = []
with open('data/gwas_catalog_v1.0.1.tsv', 'r') as gwas_fh:
    header = next(gwas_fh)
    for snp in gwas_fh:
        snp = snp.split('\t')

        # Get Coordinates
        chrom = 'chr' + str(snp[11])

        if len(snp[12]) > 0:
            coord = int(snp[12]) + 1

            # Convert from hg38 to hg19
            converted = lo.convert_coordinate(chrom, coord)
            if converted is not None:
                if len(converted) > 0:
                    new_coord = converted[0][1]
                    new_chrom = converted[0][0]
                else:
                    new_coord = 'Not Mapped'
                    new_chrom = 'Not Mapped'
            else:
                new_coord = 'Not Mapped'
                new_chrom = 'Not Mapped'
        else:
            new_coord = 'NA'
            new_chrom = 'NA'

        # Reassign back to original line
        snp[11] = str(new_chrom)
        snp[12] = str(new_coord)

        new_coordinates.append(snp)

with open('data/gwas_catalog_hg19.tsv', 'w') as gwas_fh:
    gwas_fh.write(header)
    for snp in new_coordinates:
        gwas_fh.write('\t'.join(snp))
