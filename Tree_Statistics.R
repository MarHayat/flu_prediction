library(igraph)

WIENER = TRUE		# TRUE if the Wiener index of each tree is calculated (only makes sense if DIST = TRUE and FULL = TRUE)
DOUBLE = TRUE	
### This function determines the first and last position of occurrence of a "promise" vector,
### where the promise is that it has consecutive elements, each occurs twice and min is first
findFirstPositions = function(dupVector) {
  Min = dupVector[1]
  N = length(dupVector)/2
  Max = N + Min - 1
  reps = rep(FALSE, Max)
  first = rep(NA, N)
  second = rep(NA, N)
  pos1 = 1
  pos2 = 1
  for (ind in 1:(2 * N)) {
    cur = dupVector[ind]
    if (reps[cur]) {
      second[pos2] = ind
      pos2 = pos2 + 1
    }
    else {
      first[pos1] = ind
      reps[cur] = TRUE
      pos1 = pos1 + 1
    }
  }
  output = list(first, second)
  output
}

### This function computes the sizes of the left and right subtrees rooted at internal nodes
computeLRSizes = function(tree) {
  return(computeLRValues(tree, sum))
}

### This function computes the depths of the left and right subtrees rooted at internal nodes
computeLRDepths = function(tree) {
  return(computeLRValues(tree, max))
}

### This function factory recursively computes values for subtrees rooted at internal nodes  
computeLRValues = function(tree, FUN) {
  n = Ntip(tree)
  N = 2 * n - 1
  Tab = matrix(NA, n - 1, 2)
  edges = tree$edge
  for (ind in (N - 1):1) {
    curRow = edges[ind,] - n
    pos = Tab[curRow[1], 1]
    Tab[curRow[1], 2 - is.na(pos)] = 1 + ifelse(curRow[2] <= 0, 0, FUN(Tab[curRow[2],]))
  }
  Tab
}

### This function computes the number of  leaves for subtrees rooted at internal nodes   
computeLRLeaves = function(tree) {
  n = Ntip(tree)
  N = 2 * n - 1
  Tab = matrix(NA, n - 1, 2)
  edges = tree$edge
  for (ind in (N - 1):1) {
    curRow = edges[ind,] - n
    pos = Tab[curRow[1], 1]
    Tab[curRow[1], 2 - is.na(pos)] = ifelse(curRow[2] <= 0, 1, sum(Tab[curRow[2],]))
  }
  Tab
}
### This function computes the number of  pitchforks for subtrees rooted at internal nodes  
computeLRpitchforks = function(tree) {
  n = Ntip(tree)
  N = 2 * n - 1
  Tab = matrix(NA, n - 1, 2)
  Pit = matrix(0, n - 1, 2)
  edges = tree$edge
  for (ind in (N - 1):1) {
    curRow = edges[ind,] - n
    pos = Tab[curRow[1], 1]
    Tab[curRow[1], 2 - is.na(pos)] = 1+ifelse(curRow[2] <= 0, 0, sum(Tab[curRow[2],]))
    if(Tab[curRow[1], 2 - is.na(pos)]== 5){Pit[curRow[1], 2 - is.na(pos)]=1}
    else if (Tab[curRow[1], 2 - is.na(pos)]< 5){Pit[curRow[1], 2 - is.na(pos)]=0}
    else {Pit[curRow[1], 2 - is.na(pos)]=sum(Pit[curRow[2],])}
    
  }
  Pit
}

### This function uses a BFS to compute the depth of each node (internal or leaf) in a tree
computeDepths = function(tree) {
  n = Ntip(tree)
  myGraph = as.igraph(tree)
  depths = bfs(myGraph, root = 1, neimode = "out", dist = TRUE)$dist
  depths = c(tail(depths, n), head(depths, n - 1))
  depths
}

### This function computes the LP statistic of a tree
computeLP = function(tree) {
  n = Ntip(tree)
  LRTip = computeLRLeaves(tree)
  LRPit=computeLRpitchforks(tree)
  values = (( LRTip[,1] -  LRTip[,2]))^2 + (LRPit[,1] - LRPit[,2])^2 
  output = sum(values)
  output
}

### This function computes the I2 statistic of a tree
computeI2 = function(tree) {
  n = Ntip(tree)
  LRMat = computeLRSizes(tree)
  values = abs(LRMat[,1] - LRMat[,2]) / (rowSums(LRMat) - 2)
  output = sum(values[is.finite(values)])
  output
}

### This function computes the B1 statistic of a tree
computeB1 = function(tree) {
  n = Ntip(tree)
  depths = node.depth(tree, method = 2)
  output = sum(1/(depths[-(1:(n+1))] - 1))
  output
}

### This function computes the B2 statistic of a tree
computeB2 = function(tree) {
  n = Ntip(tree)
  depths = node.depth.edgelength(tree)[1:n]
  output = sum(depths/2^depths)
  output
}

# statB2 = function(Nis, ntip){
#   B2 <- 0 
#   for (i in 1:ntip){
#     B2 <- B2 + (Nis[i] / (2^ Nis[i]))
#   }
#   return (B2)
# }


#Branching speed feature
BS=function(clade){
  tips=length(clade$tip.label)
  L=mean(node.depth.edgelength(clade)[1:tips])
  return(tips/L)
}



#mean of pairwise distance between all the tips
tips_pairwise_distance=function(clade){
  Dt=node.depth.edgelength(clade)[1:length(clade$tip.label)]
  D=outer(Dt,Dt,'-')
  D[lower.tri(D)] <- 0
  return(as.numeric(mean(abs(D))))
}
#max of pairwise distance between all the tips
tips_pairwise_distance_max=function(clade){
  Dt=node.depth.edgelength(clade)[1:length(clade$tip.label)]
  D=outer(Dt,Dt,'-')
  D[lower.tri(D)] <- 0
  return(as.numeric(max(abs(D))))
}

### This function computes the diameter of a phylo tree without multifurcations
computeDiameter = function(tree) {
  Tab = computeLRDepths(tree)
  diam = max(rowSums(Tab))
  diam
}

### This function computes the Wiener index of a phylo tree without multifurcations
computeWienerIndex = function(tree, double = DOUBLE) {
  q = rowSums(computeLRSizes(tree)) + 1
  n = length(q)
  N = 2 * n + 1
  stopifnot(q[1] == N)
  W = (1 + double) * (sum(q * (N - q)) + (n + 1) * (N - 1))
  W
}

### This function computes the betweenness centrality of a phylo tree without multifurcations
computeBetweenness = function(tree) {
  Tab = computeLRSizes(tree)
  n = nrow(Tab)
  rSums = rowSums(Tab)
  Centralities = c(rep(0, n + 1), Tab[,1] * Tab[,2] + rSums * (2 * n - rSums))
  Centralities
}

### This function computes the closeness centrality of a phylo tree without multifurcations
computeCloseness = function(tree) {
  return(1/computeFarness(tree))
}

### This function computes the farness of each node of a phylo tree without multifurcations
computeFarness = function(tree) {
  sizes = rowSums(computeLRSizes(tree))
  n = Ntip(tree)
  N = 2 * n - 1
  Farness = rep(NA, N)
  Farness[n + 1] = sum(sizes)
  edges = tree$edge
  for (ind in 1:(N - 1)) {
    curRow = edges[ind,]
    kid = curRow[2]
    Farness[kid] = N + Farness[curRow[1]] - 2 * (1 + ifelse(kid <= n, 0, sizes[kid - n]))
  }
  Farness
}


getstattest <- function(tree) {
  if (class(tree)!="phylo") tree=as(tree,"phylo")
  #    
  # lttbelow gives times, nodeids, and number lineages just below each node
  # for each node, i want to know whether one of its descendants is the NEXT one to branch (after it)
  lttb<-lttbelow(tree)
  IND<-which(lttb$nodeids > length(tree$tip.label))
  intonly=lttb$nodeids[IND] # internal nodes in chron order
  nlins=lttb$lttbelow[IND] # number of lineages below each internal node, ch order
  nexttobranch=c(intonly[-1],0) # next one to branch after each intl node, ch order
  declist <- vapply(intonly, function(x) {tree$edge[which(tree$edge[,1]==x),2]},FUN.VALUE=c(1,1)) # each internal node's 2 descendants
  colnames(declist)=as.character(intonly)     
  # For each node in intonly, I want to know: is one of its descendants next? 
  
  sresult=vapply(1:length(intonly),function(x) {is.element(nexttobranch[x],declist[,as.character(intonly[x])])},FUN.VALUE=0.5)
  
  probs=2/nlins; probs[length(probs)]=0 # the last one can't have a descendant that branches next
  W=(sum(sresult-probs))/sqrt((sum(probs*(1-probs))))
  L=length(probs); # L=round(reducefrac*L)
  return(list(details=data.frame(nodeids=intonly,s=sresult,nlins=nlins,probs=probs,topdepths=lttb$topdepths[IND],times=lttb$times[IND]),W=W))
}

# extending getstattest: I could easily ask whether ONE of the node's immediate desc is in the next m to branch. easy to code. 
# the probability works out like this: 
# there are k lineages desc from i. we throw br events randomly at lineages
# pick you fave lineage. P(doesn't get an event first time) = 1-1/k. 
# P(doesn't next time) = (1-1/(k+1))
# P(doesn't next time) = (1-1/(k+2)) 
# multiply these out; they telescope: P(doesn't in m attempts) = k-1 / (k+m -1)
# P(both desc from i don't) = (k-1 / (k+m-1)) ^2. 
# P(at least one is in the next m somewhere) = 1- above expression

# assumes no death events before the next m branching events! If there is a death event, need to adjust denominator accordingly; telescoping not as good. 
# here the P (no event) = Prod (j=0 to m-1) (1-1/k_j) where k_j is the number of lineages at the jth branching event after node i. 
# what I need to compute this: desc of i. whether they are in the next m. and, for the next m-1 events, the num lins (for the probability). 
# and for the scaling, need to work out every mth of these, so only good for big trees. 

# going to add two more options to explore: (1) above: is an immediate desc in the next m branching events? 
# and (2) are k of the next m branching events descended from node i? 

################################################
#### desc in next m NOTE THIS NEEDS TO BE CHECKED 
################################################
#' is a descendant of the node among the next m to branch? 
#' @param tree object of class phylo
#' @param m integer; we ask whether each node's desc is among the next m 
#' @return A list. details: a data frame with success (was descn in the next m), 
#' probabilities (how probable was that), nodeids (as in tree); mW: W values for all 
#' offsets (eg 1, 3, 5, .. for m=2, and 2, 4, 6 ...). W: mean of mW. 
#' @examples
#' descinm(rtree(12))
descinm <- function(tree,m=2) {
  if (!inherits(tree, "phylo")) 
    stop("Tree should be an object of class \"phylo\".")
  lttb  <- lttbelow(tree)
  IND <- lttb$nodeids > length(tree$tip.label) # LOG
  intonly=lttb$nodeids[IND] # internal nodes in chron order
  nlins=lttb$lttbelow[IND]
  L=length(intonly)
  # init outputs
  success=0*intonly # success: was descendant in next m branching events? 
  probs=0*intonly # probability of success
  # to get next m to branch, do intonly[ thisone:(thisone+m)]. nlines[same]
  for (j in 1:(L-1)) {
    dj=tree$edge[which(tree$edge[,1]==intonly[j]),2]
    nextm=intonly[ (j+1):min(j+m,L)] # next m to branch
    nlj=nlins[j:min(j+m-1,L)] # nlins now, next, next .. m-1 after
    success[j]=any(is.element(dj,nextm))
    probs[j]=prod(1 - 1/nlj)^2
  }
  mW=0*(1:m)
  for (Offset in 0:(m-1)) {
    Filter=seq(from=1+Offset,by=m,to=L-1) # ignore last entry, success 0, prob 0.
    mW[Offset+1]=(sum(success[Filter]-probs[Filter]))/sqrt(sum(probs[Filter]*(1-probs[Filter]))) 
  }
  return(list(details=data.frame(success=success,probs=probs,nodeids=intonly),mW=mW,W=mean(mW)))
}

################################################
#### are k of the next m to branch among the descendants of node i?
################################################
#' Are k of the next m to branch among the descendants of node i?
#' @param tree object of class phylo
#' @param k integer: are k of next m descs descended from node i?
#' @param m integer: are k of next m descs descended from node i?
#' @return A list.  details: a data frame with success (k desc of i in the next m), 
#' probabilities (how probable was that), nodeids (as in tree); kmW: W values for all 
#' offsets (eg 1, 3, 5, .. for m=2, and 2, 4, 6 ...). W: mean of kmW. 
#' @examples
#' kinm(rtree(12))
kinm <- function(tree,k=2,m=3) {
  if (!inherits(tree, "phylo")) 
    stop("Tree should be an object of class \"phylo\".")
  num.tips=length(tree$tip.label)
  lttb<-lttbelow(tree)
  IND<-lttb$nodeids > length(tree$tip.label) # LOG
  intonly=lttb$nodeids[IND] # internal nodes in chron order
  nlins=lttb$lttbelow[IND] # total lineages in tree below each node
  L=length(intonly) # total number of internal nodes
  mck=choose(m,k)
  success=0*intonly # success: were k of next m among i's descendants?
  probs=0*intonly # probability of success
  for ( j in 1:(L-m)) {
    dj=phytools::getDescendants(tree,intonly[j])
    Nj=(sum(dj <= num.tips) -2)  # number of br events desc from j
    #     (L-j) the  num chronologically after j 
    # 		dj=dj[dj>num.tips] # internal nodes only. may NOT be needed. 
    nextm=intonly[ (j+1):min(j+m,L)] # next m to branch
    #   		nlj=nlins[j:min(j+m-1,L)] # nlins now, next, next .. m-1 after
    success[j]=sum(is.element(dj,nextm))>=k
    probs[j]= choose(Nj,k)*choose( L-j-Nj, m-k)/choose(L-j,m) # note old version was wrong: mck*pj^k*(1-pj)^(m-k)
    # number of ways of choosing k of Nj nodes descending from this one, times # ways to choose m-k of the others, divided by all the 
    # ways to chooes m of the L-j remaining nodes
  }
  mW=0*(1:m)
  for (Offset in 0:(m-1)) {
    Filter=seq(from=1+Offset,by=m,to=L-1) # ignore last entry, success 0, prob 0.
    mW[Offset+1]=(sum(success[Filter]-probs[Filter]))/sqrt(sum(probs[Filter]*(1-probs[Filter])))
  } 
  return(list(details=data.frame(success=success,probs=probs,nodeids=intonly),kmW=mW,W=mean(mW)))
}




################################################
#### lttbelow
################################################

#' Compute times of nodes and tips and arrange in chronological order along with LTT
#' @param tree Object of class phylo, or convertible with as(tree,"phylo")
#' @return data frame with times, topdepths (topological depth from root, in number of edges, 
#' lttbelow: the number of lineages in the tree just below the current node, and nodeids: the ID of the node. 
#' 1: number tips are the tips. Order corresponds to the order in the phylo tree. 
#' @examples
#' lttbelow(rtree(10))
lttbelow <- function(tree) {
  if (class(tree)!="phylo") tree=as(tree,"phylo")
  N=2*tree$Nnode+1 # total number of nodes. Each will get a time. 
  Ntips=(N+1)/2 # number of tips
  times=vector(mode="numeric",length=N)
  times[Ntips+1]=0; # roots has time 0. this line's not needed (already 0). for clarity
  topdepths=times # initialize to all 0s. 
  # set both times and topological depths for all nodes by traversing edge list
  for (k in 1:(N-1)) {
    times[tree$edge[k,2]]=tree$edge.length[k]+ times[tree$edge[k,1]]
    topdepths[tree$edge[k,2]]=1+topdepths[tree$edge[k,1]]
  }
  # times now has all node's times. and it is in the order 1,2,3,.. n-1 where first are tips, then internals. the Ntips+1'st one is 0 because it is the root. order(times) is the index of them ordered increasingly. times[order(times)] is jjust like sort(times)
  nodeids=order(times)
  sorttimes=times[nodeids]
  
  # now I want a vector that I create by : start with 1. at each time, in order, add 1 to the previous entry if the point was an internal node. otherwise subtract 1. 
  contribs=-1+2*as.numeric(order(times)>Ntips) # each contribution to ltt
  ltt=1+cumsum(contribs)
  return(data.frame(times=sorttimes,topdepths=topdepths[nodeids], lttbelow=ltt,nodeids=nodeids))
}

#' How is the sequential branching distributed across the tips?
#' @param tree Object of class phylo
#' @return Data frame linking tip names to sequential branching along paths to root
#' @examples
#' tipprofiles(rtree(50))
tipprofiles <- function(tree) {
  if (class(tree)!="phylo") tree=as(tree,"phylo")
  N=length(tree$tip.label)
  sstuff=getstattest(tree)$details
  
  sinorder=sstuff$s[order(sstuff$nodeids)] # has the s listing in numerical order: 73, 74, 75; not in temporal order. 
  # first compute the sum of distance * s for all the nodes, along paths from root to that node:
  lttstuff=lttbelow(tree)
  rightorder=order(lttstuff$nodeids)
  
  timesinorder=lttstuff$times[rightorder[-(1:N)]]
  depthsinorder=lttstuff$topdepths[rightorder[-(1:N)]]
  
  timeprod=timesinorder*sinorder 
  topdepprod=depthsinorder*sinorder
  # now for each tip, I want the sum of timeprod and topdepprod and s, along the path from the root to the tip.
  tipsumofs=vector(mode="numeric",length=N)
  tipsumtimeprod=tipsumofs; tipsumtdprod=tipsumofs
  
  sumofs=vector(mode="numeric",length=N-1)
  stimeprod=sumofs; stopdepprod=sumofs # intermediate sums along paths
  for (k in 1:(N-1)) {
    if (tree$edge[k,2]>N) {
      sumofs[tree$edge[k,2]-N]=sumofs[tree$edge[k,1]-N]+sinorder[tree$edge[k,2]-N]
      stimeprod[tree$edge[k,2]-N]=stimeprod[tree$edge[k,1]-N]+timeprod[tree$edge[k,2]-N]
      stopdepprod[tree$edge[k,2]-N]=stopdepprod[tree$edge[k,1]-N]+topdepprod[tree$edge[k,2]-N]
    } else {
      tipsumofs[tree$edge[k,2]]=sumofs[tree$edge[k,1]-N]
      tipsumtimeprod[tree$edge[k,2]]=stimeprod[tree$edge[k,1]-N]
      tipsumtdprod[tree$edge[k,2]]=stopdepprod[tree$edge[k,1]-N]                
    }
  }
  FIXIND <- order(tree$tip.label)
  return(data.frame(SumS=tipsumofs[FIXIND], sWithLengths=tipsumtimeprod[FIXIND],sWithTopDepths=tipsumtdprod[FIXIND],TipNames=tree$tip.label[FIXIND]))
}


#  I DO NOT REMEMBER WHAT THIS WAS FOR. And I don't want to recode with phytools::getDescendants

#getspairs<- function(tree,portion=1/3)  {
#	if (class(tree)!="phylo4") tree=as(tree,"phylo4")
#	IND=bestsplit(tree,portion)
#	if (IND==0) {
#		return(list(spair=c(-1,-1),stimepair=c(-1,-1),sdeppair=c(-1,-1)))
#		} else {
#			sdetails=getstattest(tree)$details
##			ND1=descendants(tree,Ds[1],type="all");ND1=ND1[ND1>nTips(tree)]
#			ND2=descendants(tree,Ds[2],type="all"); ND2=ND2[ND2>nTips(tree)]
##			IND2=match(ND2,sdetails$nodeids)
#spair=c(sum(sdetails$s[IND1]),sum(sdetails$s[IND2]))
#stimepair=c(sum(sdetails$times[IND1]*sdetails$s[IND1]),sum(sdetails$times[IND2]*sdetails$s[IND2]))
#sdeppair=c(sum(sdetails$topdepths[IND1]*sdetails$s[IND1]),sum(sdetails$topdepths[IND2]*sdetails$s[IND2]))
#return(list(spair=spair, stimepair=stimepair,sdeppair=sdeppair))
#		}
# }



i_bl <- function(tree) {
  if (class(tree)!="phylo") tree=as(tree,"phylo")
  
  
  n<- length(tree$tip.label)
  ibl <-  tree$edge.length[tree$edge[,2]>n]
  #mean_ibl <-mean(ibl)
  #median_ibl <- median(ibl)
  #var_ibl <-var(ibl)
  #stdev_ibl <- sqrt(var_ibl)
  
  return(ibl)
}

### This function computes the diameter of a phylo tree without multifurcations
computeDiameter = function(tree) {
  Tab = computeLRDepths(tree)
  diam = max(rowSums(Tab))
  diam
  
}

### This function computes the Wiener index of a phylo tree without multifurcations
computeWienerIndex = function(tree, double = DOUBLE) {
  q = rowSums(computeLRSizes(tree)) + 1
  n = length(q)
  N = 2 * n + 1
  stopifnot(q[1] == N)
  W = (1 + double) * (sum(q * (N - q)) + (n + 1) * (N - 1))
  W
}

### This function computes the betweenness centrality of a phylo tree without multifurcations
computeBetweenness = function(tree) {
  Tab = computeLRSizes(tree)
  n = nrow(Tab)
  rSums = rowSums(Tab)
  Centralities = c(rep(0, n + 1), Tab[,1] * Tab[,2] + rSums * (2 * n - rSums))
  Centralities
}

### This function computes the closeness centrality of a phylo tree without multifurcations
computeCloseness = function(tree) {
  return(1/computeFarness(tree))
}

### This function computes the farness of each node of a phylo tree without multifurcations
computeFarness = function(tree) {
  sizes = rowSums(computeLRSizes(tree))
  n = Ntip(tree)
  N = 2 * n - 1
  Farness = rep(NA, N)
  Farness[n + 1] = sum(sizes)
  edges = tree$edge
  for (ind in 1:(N - 1)) {
    curRow = edges[ind,]
    kid = curRow[2]
    Farness[kid] = N + Farness[curRow[1]] - 2 * (1 + ifelse(kid <= n, 0, sizes[kid - n]))
  }
  Farness
}
#==================================
epitopesites=16+c(2,3,5,25,33,50,53,54,57,62,63,67,75,78,81,82,83,92,94,106,121,122,124,126,131,133,135,
                  137,142,143,144,145,146,155,156,157,158,159,160,163,164,172,173,174,186,188,189,190,192,193,196,197,
                  201,202,207,213,217,222,225,226,227,242,244,248,260,262,271,275,276,278,299,307)
# the HA1 subunit starts at aa 17; this is stated in the NCBI annotation

getEpitopeDist <- function(mytips, MappingData,hdata, Pdata, pastperiod=5, D0=14) {
  if (is.character(mytips)) { tipinds=match(mytips, hdata$tiplab)} else {tipinds=mytips }
  minTime=min(hdata$height[tipinds])
  tiplab_AA=MappingData[mytips,2]
  #find the index of the tips in MappingData
  #ind_tips=match(hdata$tiplab[mytips],MappingData[,1])
  #AAtip_labels=MappingData[ind_tips,2]
  
  reltipindex=which(hdata$height < minTime & hdata$height > minTime-pastperiod)
  
  if (length(reltipindex)==0) reltipindex=tipinds 
  # find relevant tips
  #relTips=hdata$tiplab[reltipindex]
  relTips_AA=MappingData[reltipindex,2]
  #find the AA labels and index of relevant tips
  #AAind=match(relTips,MappingData[,1])
  # AAlabels=MappingData[AAind,3]
  # figure out what to do if empty
  lmt=length(mytips); lrt=length(relTips_AA)
  # if present, remove everything but the epitope sites
  
  if (length(Pdata[[1]]) > 72) { Cur_Pdata= getEpitopeSites(Pdata[c(tiplab_AA, relTips_AA)],epitopesites)}
  
  # get distances between mytips and relTips. NOTE replace with dist.aa when ready
  mydd=dist.aa(Cur_Pdata)
  # access correct part of this info and compute my function:
  mydd=as.matrix(mydd)
  mydd=exp(-mydd[1:lmt,(lmt+1):(lmt+lrt)]/D0) # each row is a tip in my clade, and each column is a relevant other tip
  
  # compute relevant distance: for each tip in my clade I want the sum of all its dists to rel tips
  # which is the sum over the row of exp(-Dij/Do)
  cladeinfo=rowSums(mydd)/lrt # should be as many of these as there are tips in my clade; div by lrt to use mean
  # so that I can compare clades with different sizes of past tips
  names(cladeinfo)=mytips
  return(cladeinfo)
}

getEpitopeSites <- function(sdata,epitopesites){
  newdata=t(sapply(sdata, function(x) x[epitopesites]))
  newdata=as.AAbin(newdata)
  return(newdata)
}

# getEpitopeSites <- function(sdata,epitopesites){
#   newdata=t(sapply(sdata, function(x) as.character(x)[epitopesites]))
#   newdata=as.AAbin(newdata)
#   return(newdata)
# }
#================
getClades2 <- function(rt,MinTotalSize=8, MinTrimSize=8, OverlapCutoff=0.8, TimeFrame=1.4) { 
  # set up
  nTips=length(rt$tip.label)
  myroot=nTips+1
  nnodes=nTips-1
  nodeids=(nTips+1):(nTips+nnodes)
  dfsall=dfs(graph(rt$edge),root=nTips+1)$order # igraph
  dfsnodes=as.vector(dfsall[dfsall>nTips]) # they are already in order but I don't know that this will be true for any tree
  # consider removing that
  RP=NA+(1:nnodes); names(RP)=nodeids
  ChildMatrix=t(sapply(nodeids, function(x) rt$edge[rt$edge[,1]==x, 2]))
  rownames(ChildMatrix)=nodeids # each row lists the two children of a node, names are node ids
  NodeDescOf <- function(node) { mydescs=ChildMatrix[as.character(node),]
  return(mydescs[mydescs > nTips])}
  
  RP[as.character(NodeDescOf(myroot))]=myroot # alternatively could subtract nTips from the node ids toget the row numbers
  # which might be much faster for the whole tree ; doing as.char now for clarity... hmm.
  
  # compute heights 
  allHeights=node.depth.edgelength(rt); # heights from root to node, in units of branch length in the tree. assumes timed tree
  
  # need to define something like allDates, at least for the tips, giving the dates from the metadata, since these can't be 
  # computed from heights of nodes without a timed tree
  
  # compute descendants and clade sizes 
  allD=allDescendants(rt) 
  allCladeSizes=sapply(nodeids,function(x) sum(allD[[x]] <= nTips)) 
  
  # tips within the TimeFrame for each node
  allTrimmedClades = sapply(nodeids, function(x) {  myTipDes=allD[[x]][allD[[x]]<=nTips]
  myTipTimes=allHeights[myTipDes] # here, would need something like allDates[myTipDesc]
  return(myTipDes[myTipTimes <= allHeights[x]+TimeFrame]) }) # here, replace with [myTipTimes == 2006] 
  
  # sizes of trimmed clades 
  allTrimmedSizes = sapply(nodeids, function(x) length(allTrimmedClades[[x-nTips]]))
  
  
  rejectFlag=0*(1:nnodes)
  names(rejectFlag)=nodeids
  
  
  # main loop
  for (k in 2:length(dfsnodes)) {
    ii = dfsnodes[k] 
    
    if (rejectFlag[as.character(ii)] != 1) {
      # does ii have an RP? 
      rpii=RP[as.character(ii)]
      
      # if size is small, reject the node and all its descendants
      if (allCladeSizes[ii-nTips] < MinTotalSize) { 
        rejectFlag[as.character(ii)]= 1 ; iiDescs=allD[[ii]][allD[[ii]]> nTips]; 
        rejectFlag[as.character(iiDescs)]=1 
      }
      # -- if size of Ci(T) too small but Ci is big enough, set flag for non-use of i,
      #  can still use i's descendants
      if (allTrimmedSizes[ii-nTips] < MinTrimSize & allCladeSizes[ii-nTips] >= MinTotalSize) { 
        rejectFlag[as.character(ii)]=1 # as.char or ii-nTips; same effect
        RP[as.character(NodeDescOf(ii))]= ii
      }
      
      # -- if size is big enough, check intersection of clade with relevant parent's clade 
      if (allTrimmedSizes[ii-nTips] >= MinTrimSize & allCladeSizes[ii-nTips] >= MinTotalSize) { 
        # check intersection 
        myintersect = intersect(allTrimmedClades[[rpii-nTips]], allTrimmedClades[[ii-nTips]]) 
        
        #     -- if overlap is "big", set flag for non-use and set relevant parent (RP) of i's children to RP of i. 
        # overlap is the portion of ii's trimmed clade that is contained in the parent's trimmed clade
        if (length(myintersect) > OverlapCutoff*allTrimmedSizes[[ii-nTips]] ) {
          rejectFlag[as.character(ii)]=1
          RP[as.character(NodeDescOf(ii))]=rpii
        } else {   # If overlap small - keep i, set RP of i's children to i; do not reject i
          RP[as.character(NodeDescOf(ii))]=ii
        } # end if - else on the intersection
      } # end if size is big enough
    } # end if not reject flag
  } # end main loop
  rejectFlag[as.character(myroot)]=1
  return(list(nodes=nodeids,RP=RP, sizes= allCladeSizes,trimsize= allTrimmedSizes,rejected= rejectFlag,trimclades=allTrimmedClades))
}
