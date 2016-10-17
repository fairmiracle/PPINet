# PPINet
Contructing protein-protein interaction network from multiple data sources:
- Gene expression profiles: to filter gene set in network
- Ensembl: to match gene symbols and protein id
- STRING: one of the most popular PPI databases, edges are donated by Ensembl Protein ID
- BioGRID: one of the most popular PPI databases, edges are donated by OFFICIAL SYMBOL

The basic idea is to determine a significantly expressed gene set based on expression profiles, and match the gene symbols with associated protein ids from popular PPI databases, as the following figure shows:
<img src="https://github.com/fairmiracle/PPINet/blob/master/preprocess.png" align="center" width="600">

If you are using this package, please cite the following paper:

Dong Li et. al. Active module identification in intracellular networks using a memetic algorithm with a new binary decoding scheme.
