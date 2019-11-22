"""
Gregory Way 2019
tad_pathways/config.py

Description:
Functions for setting up TAD_Pathways configuration

Usage: Import only

    from config import config_yaml
"""

import yaml


def config_yaml(config_file):
    """
    Load the tad_pathways parameters and output file paths

    Arguments:
    config_file - the input configuration file to set arguments

    Output:
    A dictionary storing tad_pathway parameters and defaults
    """

    # Load configuration
    with open(config_file, "r") as stream:
        config = yaml.load(stream)

    return config
