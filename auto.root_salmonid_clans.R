
library(ape);library(phangorn)


auto.root = function(input.tree, outgr=c('Mmus', 'Hsap', 'Drer', 'Gacu', 'Olat')){
  
  if(class(input.tree)=='list') tree = input.tree[[1]]
  tree = unroot(input.tree)
  put.root = unlist(sapply(outgr, function(i) grep(i, tree$tip.label)))
  
  # if NO root sequences - defined by outgr parameter - midroot tree...(need phangorn library)
  if(length(put.root)==0){
    return(list(rooted.clans=midpoint(tree), root.info='midpoint', fail.root=F))
  }
  
  # fix names
  #names(put.root) <- sub('[:0-10000:]', '', names(put.root))
  
  # if mammal outgroups present
  
  if(sum(unique(names(put.root)) %in% c('Mmus', 'Hsap'))==2){
    cat('Minimum two mammal outgroups\n')
    mam.root = grep('Mmus|Hsap', tree$tip.label)
    
    
    if(is.monophyletic(tree, tips = mam.root, reroot = T)) { 
      cat('Monophyletic mammal outgroup\n')
      root.node = getMRCA(tree, mam.root)
      tr = try(root(tree, node = root.node, resolve.root = T), silent=T)
      if(class(tr)!='try-error') { return(list(rooted.clans=tr, root.info='Monophyletic_mammal_root', fail.root=F) ) }
      if(class(tr)=='try-error') { 
        cat('Problems with mammal node rooting\n....rooting with human\n')
        tr = try(root(tree, outgroup=grep('Hsap', tree$tip.label)[1], resolve.root = T), silent=T)
        if(class(tr)!='try-error') return(list(rooted.clans=tr, root.info='Human_root - mammals not monophyletic', fail.root=F))
        if(class(tr)=='try-error') return(list(rooted.clans=input.tree, root.info='Human_root - mammals not monophyletic', fail.root=T))
      }
    }
    
    if(!is.monophyletic(tree, mam.root, reroot = F)) { 
      cat('NOT Monophyletic mammal outgroup\n')
      root.node = getMRCA(tree, mam.root)
      tr = try(root(tree, node = root.node, resolve.root = T), silent = T)
      if(class(tr)!='try-error') { 
        cat('Rooting with MRCA-node of all outgroups\n')
        return(list(rooted.clans=tr, root.info='MRCA-node of all outgroups', fail.root=F))
      }
      if(class(tr)=='try-error') { 
        cat('Problems with mammal node rooting\n....rooting with a random human\n')
        tr = try(root(tree, outgroup=grep('Hsap', tree$tip.label)[1], resolve.root = T), silent=T)
        if(class(tr)!='try-error') return(list(rooted.clans=tr, root.info='Random_human', fail.root=F))
        if(class(tr)=='try-error') return(list(rooted.clans=input.tree, root.info='Random_human', fail.root=T))
      }
    }
    
  }
  if(sum(unique(names(put.root)) %in% c('Mmus', 'Hsap'))==1){
    cat('One mammal outgroup\n')
    tr = try(root(tree, outgroup= put.root[1], resolve.root = T), silent=T)
    if(class(tr)!='try-error') return(list(rooted.clans=tr, root.info='Single_mammal_outgroup', fail.root=F))
    if(class(tr)=='try-error') return(list(rooted.clans=input.tree, root.info='Single_mammal_outgroup', fail.root=T))
  }
  if(sum(unique(names(put.root)) %in% c('Mmus', 'Hsap'))==0){
    cat('No mammal outgroup\n')
    tr = try(ingroupMRCA.rooting(tree), silent = T)
    if(class(tr)!= 'try-error') return(list(rooted.clans=tr, root.info='ingroup_MRCA', fail.root=F))
    if(class(tr)=='try-error') {
      tr = try(root(tree, outgroup = put.root[1], resolve.root = T), silent = T)
      if(class(tr)!= 'try-error') return(list(rooted.clans=tr, root.info='Most_distant_teleost', fail.root=F))
      if(class(tr)=='try-error') return(list(rooted.clans=input.tree, root.info='Most_distant_teleost', fail.root=T))
    }
  }
}






tr = "((((((((Omyk2|CIGENEomyV6.41072:0.00055,Omyk|GSONMT00002897001:0.03238)0.966:0.00292,Ssal|XP_013993288.1:0.00055):0.00159,Okis|gb|GDQG01036265.1||m.63677:0.00934)0.823:0.00152,Tthy|Tthy_00032335-RD:0.00413)0.505:0.00907,(((Omyk2|CIGENEomyV6.34600:0.00054,Omyk|GSONMT00074661001:0.00055)0.984:0.0156,Ssal|XP_014052884.1:0.0029)0.966:0.00343,Tthy|Tthy_00027454-RA:0.00951)0.969:0.00438)0.999:0.00155,Eluc|XP_010902685.1:0.01043)0.816:0.0661,Mmus|ENSORLP00000000831.1:0.02215)1.000:0,Hsap|ENSGACP00000026666.1:0.01773)Root;"

tree = read.tree(text=tr)

par(mfrow=c(2,1))
plot(tree)

tree = unroot(tree)
plot(tree, type='unrooted')

plot(auto.root(tree)$rooted.clans)

