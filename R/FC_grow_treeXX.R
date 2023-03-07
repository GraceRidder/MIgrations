##### tree function
######

grow.treeX <-function(speciesMatrix, pine, abcd) {

 # pine <- ape::makeNodeLabel(allotree, method ="number")

  tree <- pine
  tampa <- unlist(speciesMatrix)
  temp <- unique(tampa)

  ##counts how many speciating species exist across the site
  trip2 <- list()

  if (length(temp) > 0) {
    for (f in 1:length(temp)) {

       trip <- c()

      for (e in 1:2) {
        trip <- append(trip, abcd[[e]])
        abcd <- abcd[-c(1:2)]
      }

      trip2[[f]] <- trip

      names=as.factor(trip)
      tree=ape::as.phylo(~names)
      tree$edge.length=rep(1,2)
      subtreetobind= tree
      subtreetobind <- ape::makeNodeLabel(subtreetobind, method ="number")
      subtreetobind$node.label <- temp[f]
      pine =ape::bind.tree(pine,subtreetobind,where=which(pine$tip.label==temp[f]))

    }
  }

  results = list(tree = pine, trip2 = trip2, abcd = abcd)

  return(results)

}




