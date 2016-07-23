#' PPI network from STRING
#' 
#' Constructing protenin-protein interaction network from STRING database, when 
#' STRING contains gene names from gene expression platform, like Saccharomyces 
#' cerevisiae
#' 
#' @param pvalues A N-dimensional vector of gene expression level, processed by limma
#' @param GeneNames A N-dimensional vector of gene symbols, pre-filtering by expression profiles
#' @param STRINGfilename The file downloaded from STRING, like 4932.protein.links.v10.txt
#' @param STRINGTheta The threshold for PPI edges selections when using STRING 
#' @param savename The file to save filtered gene symbols and network edgelist
#' 
#' @return net two or three column edge list for \code{\link{NetSimplify}} or 
#' \code{\link{PPIplot}} and p-values vector.
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @keywords STRING PPINet
#' 
#' 
PPIFromString <- function(pvalues, GeneNames, STRINGfilename, STRINGTheta = 0, savename){
    ##PPI pairs from STRING, protein name start with ENSPxxxx
    cat('Reading edge list from STRING...\n')
    PPI=read.delim(STRINGfilename,sep=' ')
    selectedrow = which(PPI[,3] >= STRINGTheta)
    PPI = PPI[selectedrow,]
    septag = strsplit(STRINGfilename,'.',fixed = TRUE)[[1]][1]
    linkssource = gsub(paste(septag,'.',sep=''),'',PPI[,1])
    linkstarget = gsub(paste(septag,'.',sep=''),'',PPI[,2])
    selectedp = union(linkssource,linkstarget)

    ##Matching gene names for further gene filtering
    cat('Filtering genes...\n')
    selectedgenes = intersect(selectedp,GeneNames)  
    furtherselectedrows = match(selectedgenes,GeneNames)
    GeneNames = GeneNames[furtherselectedrows]	
    pvalues = pvalues[furtherselectedrows]

    ##Network construction
    cat('Matching gene names and constructing...\n')
    net = matrix(0,nrow = dim(PPI)[1], ncol=3)
    j = 1
    for (i in 1:dim(PPI)[1]) {
    if (is.element(linkssource[i], GeneNames) && is.element(linkstarget[i], GeneNames)){
        net[j,1] = match(linkssource[i],GeneNames)
        net[j,2] = match(linkstarget[i],GeneNames)
        net[j,3] = PPI[i,3]
        j = j+1
    }
    }
    net = net[1:(j-1),]
    cat('Saving files...\n')
    #save(net,GeneNames,pvalues,file = paste(savename,'.RData',sep=''))
    write.table(pvalues, paste(savename,'Pvalues.dat',sep=''), row.names = FALSE, 
                col.names = FALSE, sep="\t")
    write.table(GeneNames, paste(savename,'GeneNames.dat',sep=''), 
                row.names = FALSE, col.names = FALSE,sep="\t",quote = FALSE)
    names(pvalues) = GeneNames
    return (list(net = net,p.values = pvalues))
}

#' PPI network from STRING and Ensembl
#' 
#' Constructing protenin-protein interaction network from STRING database, when 
#' protein name is needed for mapping with gene name from Ensembl
#' 
#' @param pvalues A N-dimensional vector of gene expression level, processed by limma
#' @param GeneNames A N-dimensional vector of gene symbols, pre-filtering by expression profiles
#' @param STRINGfilename The file downloaded from STRING, like 9606.protein.links.v10.txt
#' @param STRINGTheta The threshold for PPI edges selections when using STRING 
#' @param Ensemblfilename The file downloaded from ensembl (Customise your download), 
#' like mart_export.txt, first column is Associated Gene Name like MT-ND1, 
#' second column is Ensembl Protein ID like ENSP00000354687
#' @param savename The file to save filtered gene symbols and network edgelist
#' 
#' @return net two or three column edge list for \code{\link{NetSimplify}} or 
#' \code{\link{PPIplot}} and p-values vector. 
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @keywords STRING Ensembl PPINet
#' 
#' 
PPIFromStringEnsembl <- function(pvalues, GeneNames, STRINGfilename, STRINGTheta = 0, 
                                 Ensemblfilename, savename){
    ##PPI pairs from STRING, protein name start with ENSPxxxx
    cat('Reading edge list from STRING...\n')
    PPI=read.delim(STRINGfilename,sep=' ')
    selectedrow = which(PPI[,3] >= STRINGTheta)
    PPI = PPI[selectedrow,]
    septag = strsplit(STRINGfilename,'.',fixed = TRUE)[[1]][1]
    linkssource = gsub(paste(septag,'.',sep=''),'',PPI[,1])
    linkstarget = gsub(paste(septag,'.',sep=''),'',PPI[,2])
    selectedp = union(linkssource,linkstarget)
    
    ## Mapping from Ensembl.Protein.ID ENSPxxxx to Associated.Gene.Name C20orf144
    cat('Mapping gene names from Ensembl...\n')
    maps = read.delim(Ensemblfilename, sep='\t')
    maps = maps[!(is.na(maps[,1])|maps[,2]==""), ]
    selectedp = intersect(selectedp,maps[,2])
    selectedmapsrow = match(selectedp,maps[,2])
    maps = maps[selectedmapsrow,]
    interg = intersect(GeneNames,maps[,1])

    selectedmapsrow = match(interg,maps[,1])
    maps = maps[selectedmapsrow,]
    selectedp = maps[,2]	
    
    ##Matching gene names for further gene filtering
    cat('Filtering genes...\n')
    selectedgenes = match(interg,GeneNames)  
    GeneNames = GeneNames[selectedgenes]    
    pvalues = pvalues[selectedgenes]
    
    ##Network construction
    cat('Matching gene names and constructing...\n')
    net = matrix(0,nrow = dim(PPI)[1], ncol=3)
    j = 1
    for (i in 1:dim(PPI)[1]) {
        if (is.element(linkssource[i], selectedp) && is.element(linkstarget[i], selectedp)){
            ipg1 = maps[which(maps[,2] == linkssource[i]),1]
            net[j,1] = match(ipg1,GeneNames)
            ipg2 = maps[which(maps[,2] == linkstarget[i]),1]
            net[j,2] = match(ipg2,GeneNames)
            net[j,3] = PPI[i,3]
            j = j+1
        }
    }
    net = net[1:(j-1),]
    cat('Saving files...\n')
    #save(net,GeneNames,pvalues,file = paste(savename,'.RData',sep=''))
    write.table(pvalues, paste(savename,'Pvalues.dat',sep=''), row.names = FALSE, 
                col.names = FALSE, sep="\t")
    write.table(GeneNames, paste(savename,'GeneNames.dat',sep=''), 
                row.names = FALSE, col.names = FALSE,sep="\t",quote = FALSE)
    names(pvalues) = GeneNames
    return (list(net = net,p.values = pvalues))
}

#' PPI network from BioGRID
#' 
#' Constructing protenin-protein interaction network from BioGRID database
#' 
#' @param pvalues A N-dimensional vector of gene expression level, processed by limma
#' @param GeneNames A N-dimensional vector of gene symbols, pre-filtering by expression profiles
#' @param BioGRIDfilename The file downloaded from BioGRID, like 
#' BIOGRID-ORGANISM-Homo_sapiens-3.4.138.tab.txt
#' @param savename The file to save filtered gene symbols and network edgelist
#' @param column1 Source column in BioGRID file, for Saccharomyces_cerevisiae 
#' column1=1 and for Homo_sapiens column1=3
#' @param column2 Target column in BioGRID file, for Saccharomyces_cerevisiae 
#' column2=2 and for Homo_sapiens column3=4
#' 
#' @return net two or three column edge list for \code{\link{NetSimplify}} or 
#' \code{\link{PPIplot}} and p-values vector. 
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @keywords BioGRID PPINet
#' 
#' 
PPIFromBioGRID <- function(pvalues, GeneNames, BioGRIDfilename, savename, 
                           column1,column2){
    ##PPI pairs from STRING, protein name start with ENSPxxxx
    cat('Reading edge list from STRING...\n')
    PPI=read.delim(BioGRIDfilename,skip=35)
    links = PPI[,c(column1,column2)]
    linkssource = PPI[,column1]
    linkstarget = PPI[,column2]
    selectedp = union(links[,1],links[,2])
    
    selectedgenes = intersect(selectedp,GeneNames)
    
    ##Matching gene names for further gene filtering
    cat('Filtering genes...\n')
    selectedgenes = match(selectedgenes,GeneNames)  
    GeneNames = GeneNames[selectedgenes]    
    pvalues = pvalues[selectedgenes]
    
    ##Network construction
    cat('Matching gene names and constructing...\n')
    net = matrix(0,nrow = dim(PPI)[1], ncol=2)
    j = 1
#     edges <- links[,1] %in% GeneNames & links[,2] %in% GeneNames
#     net <- links[edges,]
    for (i in 1:dim(PPI)[1]) {
        if (is.element(linkssource[i], GeneNames) && is.element(linkstarget[i], GeneNames)){
            net[j,1] = match(linkssource[i],GeneNames)
            net[j,2] = match(linkstarget[i],GeneNames)
            j = j + 1
        }
    }
    net = net[1:(j-1),]
    cat('Saving files...\n')
    #save(net,GeneNames,pvalues,file = paste(savename,'.RData',sep=''))
    write.table(pvalues, paste(savename,'Pvalues.dat',sep=''), row.names = FALSE, 
                col.names = FALSE, sep="\t")
    write.table(GeneNames, paste(savename,'GeneNames.dat',sep=''), 
                row.names = FALSE, col.names = FALSE,sep="\t",quote = FALSE)
    names(pvalues) = GeneNames
    return (list(net = net,p.values = pvalues))
}

#' PPI network simplification
#' 
#' Simplify protenin-protein interaction network, for module identification
#' 
#' @param net Edgelist from \code{\link{PPIFromString}}, 
#' \code{\link{PPIFromStringEnsembl}} and \code{\link{PPIFromBioGRID}}
#' @param p.values The p-values vector
#' @param savename Name for saving PPI edge list as plain file and also saving 
#' nodes in largest connected component and figure as eps format by default
#' 
#' @return a list containing module size, best score, module as a list of nodes
#' 
#' @seealso \code{\link{PPIFromString}}, \code{\link{PPIFromStringEnsembl}} 
#' and \code{\link{PPIFromBioGRID}}
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @keywords STRING CONSTRUCTION
#' 
#' @examples
#' require(genalg)
#' Sp <- NetSimplify(net,0)
#' GA_result <- GA_search_connected(lambda=0.5,Sp$node_score,Sp$edge_score,
#' Sp$Edgelist,num_iter=1000, muCh=0.05, zToR=10, minsize=10)
#' 
NetSimplify <- function(net, p.values, STRINGTheta = 0){
    if (dim(net)[2] > 2 & STRINGTheta > 0)
        net = net[which(net[,3] > STRINGTheta),]
    GeneNames = names(p.values)
    el = cbind(GeneNames[net[,1]],GeneNames[net[,2]])
    #el = apply(net[,1:2], 2, as.character)
    require(igraph)
    g <- graph.edgelist(el, directed = FALSE)
    sg = simplify(g)        
    comps <- decompose.graph(sg,min.vertices=3)
    
    comps_edgelist = get.edgelist(comps[[1]], names=TRUE)
    Escores = numeric(length = dim(comps_edgelist)[1])
    if (dim(net)[2] == 3){
        edges <- comps_edgelist[,1] %in% GeneNames & comps_edgelist[,2] %in% GeneNames
        Escores = net[edges,3]
        Escores = Escores/max(Escores)
    }
    compsids = match(V(comps[[1]])$name,GeneNames)
    Nscores = p.values[compsids]
    names(Nscores) = GeneNames[compsids]
    return (list(node_score = Nscores, edge_score = Escores, 
                 Edgelist = comps_edgelist))
}

#' PPI network visulization
#' 
#' Plot protenin-protein interaction network by igraph
#' 
#' @param net Edgelist from \code{\link{PPIFromString}}, 
#' \code{\link{PPIFromStringEnsembl}} and \code{\link{PPIFromBioGRID}}
#' @param STRINGTheta The threshold for PPI edges selections when using STRING 
#' @param savename Name for saving PPI edge list as plain file and also saving 
#' nodes in largest connected component and figure as eps format by default
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @keywords STRING CONSTRUCTION
#' 
#' 
PPIplot <- function(net, STRINGTheta = 0, savename){
    if (dim(net)[2] > 2 & STRINGTheta > 0)
        net = net[which(net[,3] > STRINGTheta),1:2]
    el = apply(net[,1:2], 2, as.character)
    require(igraph)
    g <- graph.edgelist(el, directed = FALSE)
    sg = simplify(g)
    comps <- decompose.graph(sg,min.vertices=3)
    comps_edgelist = get.edgelist(sg, names=TRUE)
    write.table(comps_edgelist, file = paste(savename,'Edgelist.dat',sep=''), 
                row.names = FALSE, col.names = FALSE,sep="\t",quote = FALSE)
    write.table(V(comps[[1]])$name, file = paste(savename,'Compnent.dat',sep=''), 
                row.names = FALSE, col.names = FALSE, sep="\t", quote = FALSE)
    postscript(paste(savename,'.eps',sep=''), fonts=c("serif", "Palatino"))
    plot(comps[[1]],edge.width=0.1,vertex.label = NA,vertex.size=0.5)
    dev.off()
}