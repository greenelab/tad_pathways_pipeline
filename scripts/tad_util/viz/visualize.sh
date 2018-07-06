#!/bin/bash

python scripts/viz/visualize_TAD_locations.py -t 'hESC' -g 'hg19'
python scripts/viz/visualize_TAD_locations.py -t 'IMR90' -g 'hg19'
python scripts/viz/visualize_TAD_locations.py -t 'mESC' -g 'mm9'
python scripts/viz/visualize_TAD_locations.py -t 'cortex' -g 'mm9'

python scripts/viz/gc_content_distribution.py -g 'hg'
python scripts/viz/gc_content_distribution.py -g 'mm'

