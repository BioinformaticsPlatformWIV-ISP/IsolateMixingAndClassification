This repository contains the code used in the following publication:

> Van Uffelen, A., Posadas, A., Fraiture, M.-A., Roosens, N. H. C., De Keersmaecker, S. C. J., Marchal, K., & Vanneste, K. (2025).  
> **Detection of *Bacillus* production strains and contaminants in food enzyme products.**  
> *Food Chemistry: Molecular Sciences*, 11, 100309.  
> https://doi.org/10.1016/j.fochms.2025.100309

This codebase is intended for reproducibility and transparency of the analyses described in the publication.  
It is not a packaged or user-friendly bioinformatics tool.

---

## Requirements

The code was tested with:

- **Python 3.10**

Later Python 3 versions may also work but have not been explicitly tested.
You can install the dependencies using:

```bash
pip install -r requirements.txt
```
---

## Downloading Database Entries

Database entries are available from the following sources:
- https://doi.org/10.5281/zenodo.15388579

  Entries with `BTDB` in their name can be downloaded from:
  - https://btyperdb.btyper.app/
    
  Other entries originate from NCBI and can be downloaded using **NCBI Datasets**:
  - https://www.ncbi.nlm.nih.gov/datasets/docs/v2/  
  Use the `GCF` and `GCA` genome accession identifiers to retrieve the corresponding assemblies.

A complete database archive is also available at:

- https://doi.org/10.5281/zenodo.17435876
