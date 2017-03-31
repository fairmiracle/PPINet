# PPINet
## Contructing protein-protein interaction network from multiple data sources:

- Gene expression profiles: to filter gene set in network
- Ensembl: to match gene symbols and protein id
- STRING: one of the most popular PPI databases, edges are donated by Ensembl Protein ID
- BioGRID: one of the most popular PPI databases, edges are donated by OFFICIAL SYMBOL

The basic idea is to determine a significantly expressed gene set based on expression profiles, and match the gene symbols with associated protein ids from popular PPI databases, as the following figure shows:
<img src="https://github.com/fairmiracle/PPINet/blob/master/preprocess.png" align="center" width="600">

## Contructing Multilayer (dynamic) protein-protein interaction network

A Multilayer PPI network consists of several layers, each layer is a PPI network. There are two ways to construct Multilayer PPI:

- The structures of each layer are exactly the same as constructed above, but the vertecies weights are dertermined by gene expression level
- The structures of each layer are also determined by vertecies weights (e.g. gene expression level)

## Reference

If you are using this package, please cite the following paper:

Active module identification in intracellular networks using a memetic algorithm with a new binary decoding scheme. *BMC Genomics* 2017 18(Suppl 2):209. [Link](http://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-3495-y)
