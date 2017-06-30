
library(ape); library(phangorn)

# kjapp funksjon for å plotte med farger
plot_mdv = function(tre, tips, ...) 
    plot(tre, tip.color = ifelse(tre$tip.label %in% tips, 'red', 'black'), ...)

#input: "tre" = urota tre; "ut" = character med utgruppe-artsnavn; "summaryFUN" = funksjon (ikke implementert)

clanFinder = function(tre, ut, minSize=3, nameSep="|", output=c("tree", "tips"), summaryFUN=NULL) { 
    clans = getClans(tre) # phangorn
    tipnames = colnames(clans)
    specs = sapply(strsplit(tipnames, nameSep, fixed=T), "[", 1) # species names (default: split on pipe)

    clanList_all = apply(clans, 1, function(r) which(as.logical(r)))
    clanList_in = clanList_all[vapply(clanList_all, function(cln) length(cln) >= minSize && !any(specs[cln] %in% ut), logical(1))]
    clanList_ut = clanList_all[vapply(clanList_all, function(cln) length(cln) >= 2 && all(specs[cln] %in% ut), logical(1))]
    
    # test if a clan (given as a subset of 1:ncol(clans)) is the union of an in-clan and an out-clan 
    isGoodClan = function(cln) {
        if(length(cln) < minSize+1) return(F)
        sp = specs[cln]             # associated species
        utgruppe = cln[sp %in% ut]  # out group clan members
        ingruppe = cln[!sp %in% ut] # in group clan members
        if(length(ingruppe) < minSize || length(utgruppe) == 0) # at least 'minsize' in and 1 out
            return(F)
        in_isClan = list(ingruppe) %in% clanList_in  # test if 'in group' members form a clan
        ut_isClan = length(utgruppe) == 1 || list(utgruppe) %in% clanList_ut # same for 'out' members (if more than 1)
        in_isClan && ut_isClan
    }
    
    tip_seq = seq_along(tre$tip.label)
    if (isGoodClan(tip_seq)) # if the whole tree has an in/out split, return it without looking further
        returnClans = list(tip_seq)
    else
        returnClans = clanList_all[vapply(clanList_all, isGoodClan, logical(1))]
    
    switch(match.arg(output),
        tree = lapply(returnClans, function(cln) drop.tip(tre, setdiff(tip_seq, cln))),
        tips = lapply(returnClans, function(cln) tipnames[cln]))
}

## Eksempel fra Lars, der bare hele treet skal returneres
#myTree <- read.tree(text = "(((((((Bd_R:1,BrDi:1):1,(Hv_R:1,HoVu:1):1):1,(MeNu1:1,MeNu2:1):2):1,StLa:4):1,NaSt:5):1,Os_R:6):1,(Zm_R:2,Sb_R:2):5):0;")
 
#clans = clanFinder(myTree, ut=c("Os_R", "Sb_R", "Zm_R"), output="tree")
#clans
#all.equal(clans[[1]], myTree)

# Gammelt eksempel fra Simen
#load('OGtrees.Rdata')
#exa = og.trees[["grp100029:"]]

#exa.clans = clanFinder(exa, ut=c("Gacu", "Drer","Mmus", "Hsap"), output="tree")

#par(mfrow=c(1,2))
#plot_mdv(exa, exa.clans[[1]]$tip.label, no.mar=F, cex=0.5)
#plot(exa.clans[[1]])

