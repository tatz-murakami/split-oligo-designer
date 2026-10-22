# split-oligo-designer

## Overview

This repository contains the code to design the oligonucleotides for mFISH3D.

The code is heavily relying on [OligoMiner](http://dx.doi.org/10.1073/pnas.1714530115) tool. 



## Software prerequisites
Tested on Ubuntu 20.04 LTS with the following version of software.

- [BLAST](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (2.13.0)
- Python 3.8

Make sure to set the path to the BLAST.

## Installing dependencies

1. Make sure you have [conda](https://docs.conda.io/en/latest/miniconda.html) installed.
2. Create a new environement (python=3.8), and install the dependencies.

```
$ git clone https://github.com/tatz-murakami/split-oligo-designer.git
$ cd split-oligo-designer
$ conda env create -f environment.yml
$ conda activate oligo

$ pip install . # install oligodesigner
```

## Database and template sequence

### Make database
1. Download and make the database for blast. The database for the cDNA  can be found at [Ensembl](http://www.ensembl.org/info/data/ftp/index.html.)
2. Generate the blast database according to the instruction of [NCBI](https://www.ncbi.nlm.nih.gov/books/NBK569841/).   
	```
	makeblastdb -in my_transcriptome.fa -parse_seqids -dbtype nucl -out my_transcriptome_db
	```
	You can find an example in data/human_transcriptome`
	
### Download fasta template
Download the fasta template from [NCBI](https://www.ncbi.nlm.nih.gov/) or else.
	You can find an example in data/slc17a7/slc17a7_hs.fasta`

## Example usage
The jupyter notebook is included `OligoDesign.ipynb`.

## Parameters
The example of the parameters are shown below.
```python
mFISH3D_param = {
    'fasta':'/home/tmurakami/src/split-oligo-designer/data/slc17a7/slc17a7_hs.fasta',
    'database':'/home/tmurakami/src/split-oligo-designer/data/human_transcriptome/human_transcriptome_db',
    'minimum_offtarget_gap':100,
    'hcr_seqs':{
            'seq_even_l':'GAGGAGGGCAGCAAACGGaa',
            'seq_odd_r':'atGAAGAGTCTTCCTTTACG',
            'seq_even_r':'',
            'seq_odd_l':''
        },
    'self_remove': True, # set this True if your template sequence appears in database.
    'export_blast': False # optional. set this True to save the BLAST result for debugging.
}
```

`minimum_offtarget_gap`
If the gap between two non-specific binding is more than minimum_offtarget_gap, the pair is not regarded to cause a off-target signal. Recommended value: 100.

`hcr_seqs`
The sequences of HCR fragments. The example design can be found in `./oligodesigner/parameters.py`

`self_remove` 
Set self_remove True if you want to remove the template sequence from
off-target analysis. Otherwise, the gene of your interest could be regarded as an off-target product.
Turn this to False when your template is not found in database (e.g. GFP).

`export_blast`
Optional, default False. Set export_blast True to save the BLAST result as `<fasta basename>_blast.csv` next to the input fasta file.
The file has the 12 standard BLAST columns plus `self_hit`, which is True for hits on the template gene. Rows with `self_hit` True are excluded from the off-target analysis. All values are False when self_remove is False.
This file is for debugging and is not needed for the oligo design.



The code requires the parameters for OligoMiner. The example is below.
```python
oligominer_param = {
    'l':20,
    'L':20,
    'gcPercent':25,
    'GCPercent':75,
    'tm':20,
    'TM':100,
    'X':'AAAAAA,TTTTTT,CCCCCC,GGGGGG',
    'sal':390,
    'form':30,
    'sp':1,
    'concA':25,
    'concB':25,
    'headerVal':None,
    'bedVal':False,
    'OverlapModeVal':False,
    'verbocity':False,
    'reportVal':True,
    'debugVal':False,
    'metaVal':False,
    'outNameVal':None,
    'nn_table':'DNA_NN3'
}
```


## Selection
You will obtain multiple output files including output files from OligoMiner. The final sequences are found in (your_fasta_name)_oligosets.csv.  
The file contains the binding sequences with HCR reaction sites. You may need to further select the oligonucleotides. Since the choice of the number of oligos, binding positions, and combinations are highly dependent on the scientific question, the size of mRNA, and the research budget, I intentionally did not automate this selection process.  
My recommended workflow for the selections of oligo is as follows:
1. If the number of oligos in the oligosets.csv is less than 24 (12 pairs), select all oligos.  
2. If the number is more than 24, you may want to narrow down the binding positions. In my experience, I have never seen a situation where you need more than 48 oligos. To narrow down the binding positions, survey the past FISH literature that provides the oligonucleotide sequence of your interest. For mouse and human brain, Allen Brain ISH database (https://mouse.brain-map.org/ or https://human.brain-map.org/ish/search) is a good place to start. Use the binding positions where the past literature has used, and exclude the oligos which do not bind to the positions. If there are no past literature available, skip this step.
3. Select the oligos with the small intervals. "interval_after" indicates the distance between the oligo and the next oligo. You should select 24 to 48 oligos while keeping the intervals as small as possible. 
4. If you could not get 24 oligos after the step 2 and 3, include oligos you have disregarded. The step 2 and 3 are not the absolute criteria. You can relax the selection criteria until you get 24 oligos.

## Citation

> Murakami and Heintz. Multiplexed and scalable cellular phenotyping toward the standardized three-dimensional human neuroanatomy. bioRxiv (2022) https://doi.org/10.1101/2022.11.23.517711 

> Beliveau, B.J., Kishi, J.Y., Sasaki, H.M. et al. OligoMiner provides a rapid, flexible environment for the design of genome-scale oligonucleotide in situ hybridization probes. PNAS (2018) https://doi.org/10.1073/pnas.1714530115
	


## License

We provide this open source software without any warranty under the [MIT license](https://opensource.org/licenses/MIT).
