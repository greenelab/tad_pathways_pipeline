"""
Gregory Way 2019
tad_pathways/tad.py

Description:
Functions for setting up TAD_Pathways configuration

Usage: Import only

    from tad import TADPathways
"""

import os
import subprocess


class TadPathways:
    """
    Class to perform a tad_pathways analysis
    """

    def __init__(
        self,
        snp_list_name,
        snp_list_location,
        base_dir,
        remove_hla=False,
        all_pathways=False,
    ):
        """
        Arguments:
        snp_list_name - the name of the input snp list
        snp_list_location - the file path pointing to the location of the snp list
        base_dir - the root directory of where to save files
        remove_hla - boolean if HLA TAD should be removed
        all_pathways - boolean if consider all pathways, instead of just top pathway
        """
        self.snp_list_name = snp_list_name
        self.snp_list_location = snp_list_location
        self.base_dir = base_dir
        self.remove_hla = remove_hla
        self.all_pathways = all_pathways
        self._create_filenames()

    def _create_filenames(self):
        """
        Create several variables that correspond to the specific snp list
        """

        if not os.path.exists(self.base_dir):
            os.makedirs(self.base_dir)

        # File to compile RSid lookup results
        self.snp_location_file = os.path.join(
            self.base_dir, "{}_location.tsv".format(self.snp_list_name)
        )

        # Files for TAD based genesets - all genes that fall in candidate variant TADs
        self.genelist_file = os.path.join(
            self.base_dir, "{}_tad_genelist.tsv".format(self.snp_list_name)
        )

        # Files from WebGestalt output
        self.pathway_p_values_file = os.path.join(
            self.base_dir, "{}_pvals.tsv".format(self.snp_list_name)
        )

        # Files used to summarize results
        self.nearest_gene_file = os.path.join(
            self.base_dir, "{}_tad_genelist_nearest_gene.tsv".format(self.snp_list_name)
        )

        self.gene_evidence_file = os.path.join(
            self.base_dir, "{}_gene_evidence.csv".format(self.snp_list_name)
        )

        self.evidence_summary_file = os.path.join(
            self.base_dir, "{}_gene_evidence_summary.tsv".format(self.snp_list_name)
        )

        self.all_pathways_evidence = os.path.join(
            self.base_dir,
            "{}_all-sig-pathways_gene_evidence.csv".format(self.snp_list_name),
        )

        self.all_pathways_summary = os.path.join(
            self.base_dir,
            "{}_all-sig-pathways_gene_evidence_summary.tsv".format(self.snp_list_name),
        )

    def build_snp_list(self):
        command_list = [
            "Rscript",
            "--vanilla",
            "scripts/build_snp_list.R",
            "--snp_file",
            self.snp_list_location,
            "--output_file",
            self.snp_location_file,
        ]
        subprocess.call(command_list)

    def build_custom_tad_genelist(self):
        command_list = [
            "python",
            "scripts/build_custom_tad_genelist.py",
            "--snp_data_file",
            self.snp_location_file,
            "--output_file",
            self.genelist_file,
        ]

        if self.remove_hla:
            command_list += ["--remove_hla_tad"]

        subprocess.call(command_list)

    def webgestalt_run(self):
        command_list = [
            "Rscript",
            "--vanilla",
            "scripts/webgestalt_run.R",
            "--tad_genelist_file",
            self.genelist_file,
            "--output_name",
            self.snp_list_name,
            "--output_directory",
            self.base_dir,
        ]
        subprocess.call(command_list)

    def check_webgestalt(self):
        output_pval_file = os.path.join(
            self.base_dir, "{}_pvals.tsv".format(self.snp_list_name)
        )
        return os.path.exists(output_pval_file)

    def get_evidence(self):
        # The first command builds the evidence
        command_list = [
            "python",
            "scripts/construct_evidence.py",
            "--trait",
            self.snp_list_name,
            "--gwas",
            self.nearest_gene_file,
            "--pathway_file",
            self.pathway_p_values_file,
            "--results_directory",
            self.base_dir,
            "--gestalt_directory",
            self.base_dir,
        ]
        if self.all_pathways:
            command_list += ["--pathways", "all"]
            evidence_file = self.all_pathways_evidence
            output_file = self.all_pathways_summary
        else:
            evidence_file = self.gene_evidence_file
            output_file = self.evidence_summary_file

        subprocess.call(command_list)

        # The second command summarizes this evidence
        command_list = [
            "python",
            "scripts/summarize_evidence.py",
            "--evidence",
            evidence_file,
            "--snps",
            self.genelist_file,
            "--output_file",
            output_file,
        ]
        subprocess.call(command_list)

    def compile_results(self):
        self.build_snp_list()
        self.build_custom_tad_genelist()
        self.webgestalt_run()
        webgestalt_exists = self.check_webgestalt()
        if webgestalt_exists:
            self.get_evidence()
