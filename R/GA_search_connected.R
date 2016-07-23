#' Genetic algorithm for module identification on weighted network
#' 
#' Using genetic algorithm (GA) to identify active modules on weighted network. 
#' The connectedness of resulted module is ensured by connected components 
#' finding of binary encoded set of genes. Rewritten from COSINE package.
#' 
#' @param lambda weight parameter in the objective
#' @param node_score A N-length vector storing the nodes
#' @param edge_score A M-length vector storing the edges
#' @param EdgeList A two-columns matrix with each row a pair of weighted edges
#' @param num_iter Number of iterations in GA
#' @param muCh Mutation rate in GA
#' @param zToR Zero to one ratio
#' @param minsize The minimal size of module
#' 
#' @return a list containing module size, best score, module as a list of nodes 
#' and a list of GA functions.
#' 
#' @references Ma, Haisu, et al. "COSINE: COndition-SpecIfic sub-NEtwork 
#' identification using a global optimization method." 
#' Bioinformatics 27.9 (2011): 1290-1298.
#' @references Dong Li et al. Memetic algorithm for finding active connected subnetworks 
#' in intracellular networks. 2016.
#' 
#' @author Dong Li, \email{dxl466@cs.bham.ac.uk}
#' @keywords GA, module identification
#' 
#' @examples
#' library(COSINE)
#' data(scaled_node_score)
#' data(scaled_edge_score)
#' data(PPI)
#' GA_result <- GA_search_connected(lambda=0.5,scaled_node_score,
#' scaled_edge_score,PPI,num_iter=100, muCh=0.05, zToR=10, minsize=50)
#' ## visualized by igraph
#' selected = c()
#' for (i in 1:dim(PPI)[1]){
#'     if (PPI[i,1] %in% GA_result2$Subnet || PPI[i,2] %in% GA_result2$Subnet)
#'         selected = c(selected,i)
#' }
#' library(igraph)
#' g <- graph.data.frame(PPI[selected,], directed=FALSE)
#' V(g)$color <- "blue"
#' V(g)$color[match(GA_result2$Subnet,V(g)$name)] <- "red"
#' layout <- layout.reingold.tilford(g, circular=T)
#' plot(g,layout=layout, vertex.size=5, vertex.label.cex=5,vertex.label.dist=0.5,edge.width=5)

#' @export
#' 
GA_search_connected <- function(lambda,node_score,edge_score,EdgeList,
                                num_iter=1000, muCh=0.05, zToR=10, minsize=10){
    
    ## define the objective scoring function for condition specific subnetwork 
    all_genes<<-names(node_score)
    
    subset_extract<-function(sub){    
        genes = all_genes[sub]
        edges <- EdgeList[,1] %in% genes & EdgeList[,2] %in% genes
        
        g <- graph.data.frame(EdgeList[edges,], directed=FALSE)
        comps <- decompose.graph(g)
        
        m<-sum(edges)
        total_score<- -10000
        if(m>0){
            for (i in 1:length(length(comps))) {
                subgenes = names(V(comps[[i]]))
                subedges = EdgeList[,1] %in% subgenes & EdgeList[,2] %in% subgenes
                node_score<-sum(node_score[match(subgenes,all_genes)])/sqrt(length(V(comps[[i]])))
                edge_score <-sum(edge_score[subedges])/sqrt(length(subedges))
                
                max_score<- lambda*edge_score + (1-lambda)*node_score
                
                if (max_score > total_score){
                    total_score = max_score
                    targetgenes = subgenes
                }
                
            }
            
        }
        return (targetgenes)
    }
    
    subset_score<-function(sub){    
        genes<-all_genes[sub==1]
        n<-length(genes)    
        if(n<minsize){return(10000)}
        else{
            edges<- EdgeList[,1] %in% genes & EdgeList[,2] %in% genes
            m<-sum(edges)
            total_score<- -10000
            if(m>2){
                #print(m)
                g <- graph.data.frame(EdgeList[edges,], directed=FALSE)
                comps <- decompose.graph(g)
                
                for (i in 1:length(length(comps))) {
                    subgenes = names(V(comps[[i]]))
                    subedges = EdgeList[,1] %in% subgenes & EdgeList[,2] %in% subgenes
                    node_score <-sum(node_score[match(subgenes,all_genes)])/sqrt(length(V(comps[[i]])))
                    edge_score <-sum(edge_score[subedges])/sqrt(length(subedges))
                    max_score <- lambda*edge_score + (1-lambda)*node_score
                    
                    #if (node_score > total_score)
                    #  total_score = node_score
                    if (max_score > total_score)
                        total_score = max_score
                }
                
            }
            return (-total_score)
        }
    }
    
    monitor <- function(obj) {
        minEval = min(obj$evaluations);
        filter = obj$evaluations == minEval;
        #   print(table(filter))
        bestObjectCount = sum(rep(1, obj$popSize)[filter]);
        if (bestObjectCount > 1) {
            bestSolution = obj$population[filter,][1,];
        } else {
            bestSolution = obj$population[filter,];
        }
        
        outputBest = paste(obj$iter, " #selected=", sum(bestSolution),
                           " Best (Score=", -minEval, "):\n", sep="");
        print(outputBest)
    }
    
    
    ##Start search
    
    gene_num <- length(node_score)
    print(paste("Working on lambda=",lambda)) 
    GA_result <- rbga.bin(size=gene_num,evalFunc=subset_score,iters=num_iter,mutationChance=muCh,monitorFunc=monitor,zeroToOneRatio=zToR)      
    a <- which.min(GA_result$evaluations)
    final <- GA_result$population[a,]
    b <- which(final==1)
    
    optimal_subnet = subset_extract(b)
    num_gene_selected <- length(optimal_subnet)
    
    best_score <- (-1)*min(GA_result$evaluations)
    print(paste("Finished lambda=",lambda))
    
    return(list( Subnet_size = num_gene_selected, Best_Scores = best_score, Subnet = optimal_subnet, GA_obj = GA_result))
}