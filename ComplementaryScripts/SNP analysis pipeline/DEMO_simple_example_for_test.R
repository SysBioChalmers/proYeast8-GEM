
#input the seq and mutation
seq0 <- c(1,2,3,4)
mutation0 <- c(0,0,10,10)


#input the distance matrix
c1 <- c(0,20,5,7)
c2 <- c(20,0,3,6)
c3 <- c(5,3,0,1)
c4 <- c(7,6,1,0)
distance0 <- rbind(c1,c2,c3,c4)


#calculate the WAP
#1,3,4
m1=0
m3=10
m4=10

n1=m1^3/(2^3+m1^3)
n3=m3^3/(2^3+m3^3)
n4=m4^3/(2^3+m4^3)


d13=distance0[1,3]
d14=distance0[1,4]
d34=distance0[3,4]


#constant
t=6

#calculate wap for each pair
wap13=n1*n3*exp(-d13^2/2/t^2)
wap14=n1*n4*exp(-d14^2/2/t^2)
wap34=n3*n4*exp(-d34^2/2/t^2)

wap = wap13+wap14+wap34


getWAP <- function(d0){
  #d0, a vector contains the distance
  t=6
  #calculate wap for each pair
  wap13=n1*n3*exp(-d0[1]^2/2/t^2)
  wap14=n1*n4*exp(-d0[2]^2/2/t^2)
  wap34=n3*n4*exp(-d0[3]^2/2/t^2)
  wap = wap13+wap14+wap34
  return(wap)
}




#change the postion of mutation while keep the mutation number in each postion
n <- 24
wap_sample <- list()
for (i in 1:n){
  sample0 <- sample(seq0, 3, replace = FALSE, prob = NULL)
  print(sample0)
  combine0 <- combn(sample0, 2) #get the random combine of the postions with mutation
  ss <- vector()
  for (j in 1:ncol(combine0)){
    ss[j] = distance0[combine0[,j][1],combine0[,j][2]] # get the ditance based on the index
  }
  
  wap_sample[[i]] <-getWAP(ss)
  
  }


#calculate the p value
wap_dis <- unlist(wap_sample)
min(wap_dis)
max(wap_dis)


plot(density(wap_dis), frame = FALSE, col = "steelblue",
     main = "Density of coverage",
     xlim=c(0, 1.5),
     ylim=c(0,3),
     xlab="coverage",
     ylab="Density")


plot(ecdf(wap_dis),
     main = "Cumulative density of coverage",
     xlim=c(0, 1.5),
     ylim=c(0,1),
     xlab="Coverage",
     ylab="Cumulative")

n <- 1000
seq0 <- c(1,2,3,4,5,6)
wap_sample <- list()
for (i in 1:n){
  sample0 <- sample(seq0, 1, replace = FALSE, prob = NULL)
  wap_sample[[i]] <-  sample0
}
wap_dis <- unlist(wap_sample)

plot(density(wap_dis), frame = FALSE, col = "steelblue",
     main = "Density of coverage",
     xlim=c(0, 7),
     ylim=c(0,0.5),
     xlab="coverage",
     ylab="Density")


plot(ecdf(wap_dis),
     main = "Cumulative density of coverage",
     xlim=c(0, 7),
     ylim=c(0,1),
     xlab="Coverage",
     ylab="Cumulative")


#  A C  5
#  A C  5
#  A B  7
#  C B  8
#  E D  0
#  A A  4
S1 <- c('A','A','A','C','E','A')
S2 <- c('C','C','B','B','D','A')
S3 <- c(5,5,7,8,0,4)

links <- data_frame(from = S1, to = S2, weight = S3)
g <- graph_from_data_frame(d = links, directed = FALSE)
plot(g)
  # library(networkD3)
  # simpleNetwork(links, fontSize = 18, opacity = 1, zoom = TRUE,fontFamily = "sans-serif")
  # split the graph into subgraph and get the unique cluster
  # calculate the closeness centrality for each clust
dg <- decompose.graph(g)
detail_mutant_position0 <- list()
position_combine <- vector()
for (i in seq_along(dg)) {
    clust2 <- dg[[i]][1]
    detail_mutant_position0[[i]] <- names(clust2)
    position_combine[i] <- paste0(names(clust2), collapse = ";")
  }
  
  
closeness0 <- list()
cluster_closeness <- vector()
max_closeness <- vector()
dg <- decompose.graph(g)




#how to calculate the closeness residual
closeness.residual = function (graph, vids = V(graph), mode = c("all", "out", "in"), 
          weights = NULL) 
{
  .cs.checkPreconditions(graph)
  vids <- .cs.as.igraph.vs(graph, vids)
  res <- double()
  sp <- shortest.paths(graph, v = V(graph), mode = mode[1], 
                       weights = weights)
  for (v in V(graph)[vids]) {
    res <- append(res, sum(1/(2^sp[v, ])))
  }
  if (getIgraphOpt("add.vertex.names") && is.named(graph)) {
    names(res) <- V(graph)$name[vids]
  }
  res
}


g <- graph(c(1,2,2,3,3,4,4,2))
plot(g)
closeness.residual(g,mode = 'all')

g <- graph(c(1,2,2,3,3,4,4,2,2,5,2,6,7,8))
plot(g)
closeness.residual(g)
sp <- shortest.paths(g, v = V(g))
sp[1,]

g <- graph(c('A','B','A','E'))
weight0 <- c(10,2)

plot(g)
closeness.residual(g, weights = weight0)
sp <- shortest.paths(g, v = V(g), weights = weight0)
sp[1,]
1+1/(2^10)+1/(2^2)
