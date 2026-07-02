#' @title as.ultrametric
#' @param tree tree object
#' @param ... additional parameters
#' @return treedata or phylo object
#' @export 
as.ultrametric <- function(tree, ...){
    UseMethod("as.ultrametric")
}

#' @method as.ultrametric phylo
#' @export
## reference 
## https://github.com/PuttickMacroevolution/MCMCtreeR/blob/master/R/readMCMCTree.R
as.ultrametric.phylo <- function(tree, ...){
    edge <- tree$edge
    outer <- edge[, 2]
    inner <- edge[, 1]
    ntip <- Ntip(tree)
    root <- rootnode(tree)
    children <- split(seq_along(outer), inner)
    distance <- numeric(max(edge))
    queue <- integer(max(edge))
    queue[1] <- root
    n_queue <- 1L
    i <- 1L
    while (i <= n_queue) {
        edge_index <- children[[as.character(queue[i])]]
        n_child <- length(edge_index)
        if (n_child > 0) {
            child <- outer[edge_index]
            distance[child] <- distance[queue[i]] + tree$edge.length[edge_index]
            idx <- n_queue + seq_len(n_child)
            queue[idx] <- child
            n_queue <- n_queue + n_child
        }
        i <- i + 1L
    }
    tipindx <- which(outer <= ntip)
    totalPath <- distance[outer[tipindx]]
    addLength <- max(totalPath) - totalPath
    tree$edge.length[tipindx] <- tree$edge.length[tipindx] + addLength
    return (tree)
}

#' @method as.ultrametric treedata
#' @export
as.ultrametric.treedata <- function(tree, ...){
    tree@phylo <- as.ultrametric(tree=tree@phylo,...)
    return (tree)
}

#' @method as.ultrametric tbl_tree
#' @export
as.ultrametric.tbl_tree <- function(tree, ...){
    tree <- as.treedata(tree)
    tree <- as.ultrametric(tree, ...)
    return(tree)
}
