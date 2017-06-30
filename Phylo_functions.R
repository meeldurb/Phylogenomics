#######################
### phylo functions ###
#######################

require(ape, quietly = T, warn.conflicts = F)
require(seqinr, quietly = T, warn.conflicts = F)
require(phangorn, quietly = T, warn.conflicts = F)

#---------------------------------------fix.names---------------------------------------------#
# keeps species identifier in beginning of tip.labels separated with '_'

fix.names=function(tree, sep='_', rm.tag=c('\\.1|\\.2')){
  t = tree
  if(!is.null(rm.tag)) { t$tip.label  <- gsub(rm.tag, '', t$tip.label) }
  t$tip.label = sapply(strsplit(as.character(t$tip.label), sep), '[', 1)
  return(t)
}

fix.names.vector=function(names, sep='_', rm.tag=c('\\.1|\\.2')){
  sapply(strsplit(as.character(names), sep), '[', 1)
}


#-----------------------------------------------------------------------------------------------#


#-------------------------------------check.boot-------------------------------------------------#
## --> filter trees based on minimum bootstrap values in trees: 

check.boot = function(tree, min.bs=50){
  t = tree
  if(TRUE %in% (((as.numeric(t$node.label)))<min.bs)) return(NA)
  else return(t)
}
#-----------------------------------------------------------------------------------------------#



#-------------------------------------is.na.tree-------------------------------------------------#

is.na.tree= function(trees, rm.na=T){
  if(rm.na) { trees[!sapply(trees, function(i) is.na(i[1]), USE.NAMES = FALSE)] }
  else { sapply(trees, function(i) is.na(i[1]), USE.NAMES = FALSE) }
}


#-----------------------------------------------------------------------------------------------#




#------------------------------------root.tree----------------------------------------------------#

root.tree=function(tree, out.gr= c('Dr', 'Ga'), fix.names=T, verbose=F){
  t = tree
  if(fix.names){ t = fix.names(t, sep='_') }
  og = out.gr[which(out.gr %in% t$tip.label==T)[1]]
  if(verbose) cat('rooting with', og)
  if(sum(out.gr %in% og)<1) { return(NA) }
  res  <- try(root(t, outgroup = og, resolve=T))
  if(class(res)=='try-error') { return(NA) }
  else { res$tip.label  <- tree$tip.label ; return(res) }
}


#------------------------------------is.mono----------------------------------------------------#

is.mono=function(tree, mono.tips=c('Ssa','Om', 'Ssa','Om'), fix.names=T, rm.tag=NULL){
    if(length(tree)==1) { return(NA) }
    if(fix.names) {
                    t = fix.names(tree, sep='_', rm.tag=rm.tag)
                    return(is.monophyletic(phy=t, tips=mono.tips))
    }
  
    else { return(is.monophyletic(phy=trees, tips=mono.tips)) }

}


#------------------------------------count.salmo----------------------------------------------------#

count.salmo=function(tree, salmo.tip.tags=c('Ssa', 'Om'), for.each=F){
  if(for.each) { 
    res = c()
      for(i in 1:length(salmo.tip.tags)){
      res[i] = sum(fix.names(tree)$tip.label %in% salmo.tip.tags[i])
      }
    names(res)  <- salmo.tip.tags
    return(res)
  }
  else { return(sum(fix.names(tree)$tip.label %in% salmo.tip.tags)) }
}


#-------------------------------------topo_count------------------------------------------------#
# counts different topologies in a multiphylo-list. Require rooted trees
# input: a list of **rooted** phylo's.
# output: a list of pairs (tree, count) where 'tree' is a phylo object and 'count' is the number of occurrences.

topo_count = function(tree.list, sort=T, fix.names=T) { 
  tl = tree.list
  if(fix.names)  { tl = lapply(tl, fix.names) }
  
  res = list()
  for(tr in tl) {
    seen_before = F
    for(i in seq_along(res))
      if(all.equal(res[[i]]$tree, tr, use.edge.length = F)) {
        res[[i]]$count = res[[i]]$count + 1
        seen_before = T
        break
      }
    if(!seen_before)  res = c(res, list(list(tree=tr, count=1)))
  }
  if(sort) res = res[order(sapply(res, '[[', 2))]
  res   
}


#example: rooting the trees. Don't think edge lenghts matter here, but removing them as well.
#rooted_nolength = lapply(trees.filt, function(tr) {tr$edge.length=NULL; root(tr, 'Os', resolve=T)})
#count = topo_count(rooted_nolength)
#sapply(count, '[[', 2)



#-------------------------- extract k-values -----------------------------------------------------#

getKmatix=function(k.obj, k.value=c('ks','ka'), na.cuttoff=2){
  id = ifelse(k.value=='ks',2, 1 )
  sapply(k.obj, function(i) { m = as.matrix(dist(i[[id]], upper=T, diag=T))
                              colnames(m)  <-  fix.names.vector(colnames(m))
                              rownames(m)  <-  colnames(m)
                              m[m>na.cuttoff]  <- NA
                              m }
                )
}


pairwise.k = function(k.list, seq1, seq2){
  sapply(k.list, function(i) {
    if(sum(colnames(i) %in% c(seq1, seq2))==2) return(i[seq1, seq2])
    else return(NA)
  })
}

#--------------------------------------------------------------------------------------------------#





