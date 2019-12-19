"""
Gregory Way 2019
TAD Pathways
tad_pathways/genes/genelist.py

Description:
stores function to build genelists of all genes in TADs harboring significant SNPs

Usage:
Import Only

    from genes.genelist import TADGeneList
"""

import pandas as pd


class TADGeneList:
    """
    Class to build a genelist for TAD associated genes
    """

    def __init__(self, snp, gene_index_file, remove_hla_tad=False, hla_tad="1198"):
        """
        Arguments:
        snp - the file name, or the dataframe of the input snp list.
              The SNP DataFrame has the columns: [RSid, database source, chromosome,
              genomic position, minor allele, and group]
        gene_index_file - file location of TAD index file
        remove_hla_tad - boolean if HLA TAD should be removed
        hla_tad - string of tad identifier storing the HLA TAD
        """

        try:
            self.snp = pd.read_table(snp, sep="\t")
        except FileNotFoundError:
            self.snp = snp
        self.remove_hla_tad = remove_hla_tad
        self.hla_tad = str(hla_tad)
        self.gene_index_file = gene_index_file
        self.tad_genes_df = pd.read_table(self.gene_index_file, index_col=0)

    def assign_snp_to_tad(self, snp_signal, tad_boundary_df):
        """
        Take an input snp signal and TAD boundary coordinates to output the
        coordinates of TAD where the SNP signal resides.

        Arguments:
        snp_signal -  Pandas Series with SNP descriptive attributes
        tad_boundary_df - a Pandas DataFrame of all TADs and genomic locations

        Output:
        a DataFrame storing all genes residing in the signal TAD
        """
        chrom = str(snp_signal.chrom[3:])
        snp_loc = int(snp_signal.position)

        tad = tad_boundary_df.loc[
            (tad_boundary_df.chromosome == chrom)
            & (tad_boundary_df.TAD_start <= snp_loc)
            & (tad_boundary_df.TAD_end > snp_loc)
        ]

        return tad.reset_index(drop=True)

    def compile(self):
        """
        Given a snp dataframe and a TAD gene index file, compile all genes in the TAD
        and identify which genes are closest to the input signals
        """

        # Initialize empty objects to store info
        self.tad_genes_df = pd.DataFrame(
            columns=self.tad_genes_df.columns.tolist() + ["group"]
        )
        nearest_gene_list = []
        for group in self.snp.group.unique():
            snp_sub_df = self.snp.query("group == @group")
            for snp, snp_series in snp_sub_df.iterrows():
                # Find TAD that the SNP resides in and all genes that are in the TAD
                tad_assign_df = self.assign_snp_to_tad(snp_series, self.tad_genes_df)
                tad_id = tad_assign_df.loc[:, "TAD_id"]

                # Check to see if we should continue considering TAD genes
                if len(tad_id) == 0:
                    continue

                if self.remove_hla_tad:
                    if tad_id[0] == self.hla_tad:
                        continue

                # Only collect genes if the SNP actually resides in a TAD
                if tad_assign_df.shape[0] == 0:
                    continue

                # Subset to protein coding genes for nearest gene designation
                protein_coding = tad_assign_df.query('gene_type == "protein_coding"')

                # Find nearest gene
                # Closest to the start of the gene
                near_start_df = pd.DataFrame(
                    (protein_coding.start - snp_series["position"]).abs().sort_values()
                )

                # Closest to the end of the gene
                near_end_df = pd.DataFrame(
                    (protein_coding.stop - snp_series["position"]).abs().sort_values()
                )

                dist_df = near_start_df.join(near_end_df)

                # Find the minimum distance and the index of a protein coding gene
                min_dist = dist_df.min().min()
                near_gene_idx = dist_df[dist_df == min_dist].dropna(thresh=1).index

                # If gene falls in TAD without protein coding gene return empty
                if len(near_gene_idx) == 0:
                    nearest_gene = ""
                else:
                    near_gene_idx = near_gene_idx.tolist()[0]
                    nearest_gene = tad_assign_df.iloc[near_gene_idx,].gene_name

                # Append to nearest_gene_df
                nearest_gene_return = pd.DataFrame(
                    [nearest_gene, snp_series.snp, group]
                ).T
                nearest_gene_list.append(nearest_gene_return)

                # Assign new columns to each TAD assignment for the RSid and group
                tad_assign_df = tad_assign_df.assign(
                    custom_snp=snp_series.snp, group=group
                )
                self.tad_genes_df = self.tad_genes_df.append(
                    tad_assign_df, ignore_index=True, sort=True
                )

        # Output results
        nearest_columns = ["MAPPED_GENE", "snp", "group"]
        self.nearest_gene_df = pd.concat(nearest_gene_list)
        self.nearest_gene_df.columns = nearest_columns
        self.tad_genes_df.columns = [
            "TADEnd",
            "TADidx",
            "TADStart",
            "chrom",
            "custom_snp",
            "db",
            "gene_name",
            "gene_type",
            "group",
            "start",
            "stop",
            "strand",
            "type",
        ]
