## 0 Getting Started
## if you haven't installed the relevant packages, do so with the following code:
## install.packages("statnet")

### Section 1. A Brief R Tutorial

## 1.1 Introduction to basic R syntax
1 + 3 # evaluation
a <- 3 # assignment
a # evaluation
a<-3 # spacing does not matter
a <- 3 # spacing does not matter
a

sqrt(a) # use the square root function
b <- sqrt(a) # use function and save result
b

d # evaluate something that is not there
a == b # is a equivalent to b?
a != b # is a not equal to b?
ls() # list objects in the global environment
rm(a) # remove a single object
ls()
rm(list=ls()) # remove everything from the environment
ls()

# help is here!!
help.start() # get help with R generally
?sqrt # get specific help for a function
??sqrt # looking for help pertaining to sqrt
apropos("sq") # it's on the tip of my tongue...


## 1.2 Vectors and matrices in R
# Creating vectors using the "combine" operator
c
?c
a <- c(1,3,5) # create a vector by combining values
a
a[2] # select the second element
b <- c("one","three","five") # also works with strings
b
b[2]
a <- c(a,a) # can apply recursively
a
a <- c(a,b) # mixing types---what happens?
a # all converted to the same type

# Sequences and replication
a <- seq(from=1,to=5,by=1) # from 1 to 5 the slow way
b <- 1:5 # a shortcut!
a==b # all TRUE
rep(1,times=5) # a lot of ones
rep(1:5,times=2) # repeat an entire sequence
rep(1:5,each=2) # same, but element-wise
rep(1:5,times=5:1) # can vary the count of each element

# any, all, and which (with vectors)
a <- -1:4 # create a vector
a
a>2 # some TRUE, some FALSE
any(a>2) # are any elements TRUE?
all(a>2) # are all elements TRUE?
which(a>2) # which indices are TRUE?

# From vectors to matrices
a <- matrix(data=1:25, nrow=5, ncol=5) # create a matrix the "formal" way
a
a[1,2] # select a matrix element (two dimensions)
a[1,] # just the first row
all(a[1,]==a[1,1:5]) # show the equivalence
a[,2] # can also perform for columns
a[2:3,3:5] # select submatrices
a[-1,] # nice trick: negative numbers omit cells!
a[-2,-2] # get rid of row two, column two

b <- cbind(1:5,1:5) # another way to create matrices
b
d <- rbind(1:5,1:5) # can perform with rows, too
d
cbind(b,d) # no go: must have compatible dimensions!
dim(b) # what were those dimensions, anyway?
dim(d)
NROW(b)
NCOL(b)
cbind(b,b) # combining two matrices

t(b) # can transpose b
cbind(t(b),d) # now it works
rbind(t(b),d) # now it works


## 1.3 Element-wise operations
a <- 1:5
a + 1 # addition
a * 2 # multiplication
a / 3 # division
a - 4 # subtraction
a ^ 5 # you get the idea...

a + a # also works on pairs of vectors
a * a
a + 1:6 # problem: need same length

a <- rbind(1:5,2:6) # same principles apply to matrices
b <- rbind(3:7,4:8)
a + b
a / b
a %*% t(b) # matrix multiplication

# logical operators (generally) work like arithmetic ones
a > 0 # each value greater than zero?
a == b # corresponding values equivalent?
a != b # corresponding values not equivalent?
!(a == b) # same as above
(a>2) | (b>4) # the OR operator
(a>2) & (b>4) # the AND operator
(a>2) || (b>4) # beware the "double-pipe"!
(a>2) && (b>4) # (and the "double-ampersand"!)

# ditto for many other basic transformations
log(a)
exp(b)
sqrt(a+b) # note that we can nest statements!
log((sqrt(a+b)+a)*b) # as recursive as we wanna be


## 1.4 Data Frames
d <- data.frame(income=1:5,sane=c(T,T,T,T,F),name=LETTERS[1:5])
d
d[1,2] # acts a lot like a matrix!
d[,1]*5
d[-1,]
d$sane # can use dollar sign notation
d$sane[3]<-FALSE # making changes
d
d[2,3] # shows factors for string values
d[2,3]<-"F"

# if you want to do without factors
d$name <- LETTERS[1:5] # eliminate evil factors by overwriting
d[2,3]

# or, create it without factors
d <- data.frame(income=1:5,sane=c(T,T,T,T,F),name=LETTERS[1:5],stringsAsFactors=FALSE)
d[2,3]<-"F"
d[2,3]
d

d <- as.data.frame(cbind(1:5,2:6)) # can create from matrices
d
is.data.frame(d) # how can we tell it's not a matrix?
is.matrix(d) # the truth comes out


## 1.5 Finding built-in data sets
# Many packages have built-in data for testing and educational purposes
data() # lists them all
?USArrests # get help on a data set
data(USArrests) # load the data set
USArrests # view the object


## 1.6 Elementary visualization
# R's workhorse is the "plot" command
plot(USArrests$Murder,USArrests$UrbanPop) # using dollar sign notation
plot(USArrests$Murder,USArrests$UrbanPop,log="xy") # log-log scale

# Adding plot title and axis labels
plot(USArrests$Murder,USArrests$Assault,xlab="Murder",ylab="Assault",main="USArrests")

# Can also add text
plot(USArrests$Murder,USArrests$Assault,xlab="Murder",ylab="Assault", main="USArrests",type="n")
text(USArrests$Murder,USArrests$Assault,rownames(USArrests))

# Histograms and boxplots are often helpful
hist(USArrests$Murder)
boxplot(USArrests)



###########################################################
###########################################################
##2 network objects: importing, exploring, and manipulating network data
###########################################################
##2.1 Starting from scratch
##First, in Statnet:
library(statnet)

#Adding Edges
g <- network.initialize(5) # Create an empty graph
g[1,2] <- 1 # Add an edge from 1 to 2
g[2,] <- 1 # Add edges from 2 to everyone else
g # Examine the result

m <- matrix(0, nrow=5, ncol=5) # Create an adjacency matrix
m[3,4:5] <- 1 # Add entries from 3 to 4 and 5
g[m>0] <- 1 # Add more entries
g

#Delete edges
g[3,5] <- 0 # Remove the edge from 3 to 5
g # Its gone!
g[,] <- 0 # Remove all edges
g # Now, an empty graph

##Now, in igraph
library(igraph)

# The function make_graph() can work two ways:
# 1. as a vector of an even number of node IDs where each pair of adjacent nodes is considered an edge
# illustrated here with spacing
make_graph( c(1,2, 2,3, 3,4), directed=FALSE)
#or 2. using a symbolic formula
# `-` undirected tie
# `--+` directed tie (`+` is arrow's head)
# `:` refer to node sets (e.g. `A -- B:C` creates ties `A -- B` and `A -- C`)
# A network is either directed or undirected, it is not possible mix directed and undirected ties
g1 <- make_graph(~ A - B, B - C:D:E)
g2 <- make_graph(~ A --+ B, B +-- C, A --+ D:E, B --+ A)

g2
#To interpret the printout:
#First line includes
# - Class of the object (`IGRAPH`)
# - An ID of the object, not of particular interest (see also `?igraph::graph_id`)
#  A set of four slots for letter codes indicating, in order:
#  `U` or `D` if the network is `U`ndirected or `D`irected
#  `N` if the nodes have names
#  `W` if the network is weighted
#  `B` if the network is bipartite
#   Number of nodes
#   Number of edges
#   Starting with `+ attr:` list of present attributes, each of the form `nameoftheattribute (x/y)` where 
#   `x` informs about the type of an attribute: `v`ertex, `e`dge or `g`raph attribute
#   `y` informs about the mode of an attribute: `n`umeric, `c`haracter, `l`ogical, or e`x`tended (e.g. lists)
#   Starting with `+ edges:` a list of (some of) the edges of the network

###########################################
##2.2 Importing Relational Data
#Be sure to be in the directory where you stored the data for the workshop. If you've not set the working directory, you must do so now. See Section 1.7 for how to do this.
#You can use setwd() to change it, or platform specific shortcuts
getwd() # Check what directory you're in
list.files() # Check what's in the working directory

#Read an adjacency matrix (R stores it as a data frame by default)
relations <- read.csv("relationalData.csv",header=FALSE,stringsAsFactors=FALSE)
relations
#Here's a case where matrix format is preferred
relations <- as.matrix(relations) # convert to matrix format

#Read in some vertex attribute data (okay to leave it as a data frame)
nodeInfo <- read.csv("vertexAttributes.csv",header=TRUE,stringsAsFactors=FALSE)
nodeInfo

#Since our relational data has no row/column names, let's set them now
rownames(relations) <- nodeInfo$name
colnames(relations) <- nodeInfo$name
relations

#Why did we set stringsAsFactors=FALSE?
relations_wFactors<-read.csv("vertexAttributes.csv",header=T,stringsAsFactors=TRUE)
relations_wFactors[1,2]<-"Derek"

numFactor<-as.factor(c(1,3,5))
numFactor+1
as.numeric(numFactor)

#For more information. . .
?list.files
?read.csv
?as.matrix
?rownames

###########################################################
##2.3 Creating network objects

nrelations<-network(relations,directed=FALSE) # Create a network object based on relations
nrelations # Get a quick description of the data
nempty <- network.initialize(5) # Create an empty graph with 5 vertices
nempty # Compare with nrelations

#now in igraph:
library(igraph)
irelations<-graph_from_adjacency_matrix(as.matrix(relations),mode="directed")

#For more information. . .
?network 
?as.network.matrix

###########################################################
##2.4 Description and Visualization

summary(irelations) #summary
ecount(irelations) #How many edges
vcount(irelations) #How many vertices

plot(irelations,
     layout=layout_with_fr,
     vertex.color="red",
     vertex.size=15,
     edge.arrow.size=0.5,
     vertex.label.color="black")

detach("package:igraph")

summary(nrelations) # Get an overall summary
network.dyadcount(nrelations) # How many dyads in nflo?
network.edgecount(nrelations) # How many edges are present?
network.size(nrelations) # How large is the network?
as.sociomatrix(nrelations) # Show it as a sociomatrix
nrelations[,] # Another way to do it

plot(nrelations,displaylabels=T) # Plot with names
plot(nrelations,displaylabels=T,mode="circle") # A less useful layout...

library(sna) # Load the sna library
gplot(nrelations) # Requires sna
gplot(relations) # gplot Will work with a matrix object too

#For more information. . .
?summary.network
?network.dyadcount
?network.edgecount
?as.sociomatrix
?as.matrix.network
?is.directed



###########################################################
##2.5 Network and Vertex Attributes 
#Add some attributes
nrelations %v% "id" <- nodeInfo$id # Add in our vertex attributes
nrelations %v% "age" <- nodeInfo$age
nrelations %v% "sex" <- nodeInfo$sex
nrelations %v% "handed" <- nodeInfo$handed
nrelations %v% "lastDocVisit" <- nodeInfo$lastDocVisit

#Listing attributes
list.vertex.attributes(nrelations) # List all vertex attributes
list.network.attributes(nrelations) # List all network attributes

#Retrieving attributes
nrelations %v% "age" # Retrieve vertex ages
nrelations %v% "id" # Retrieve vertex ids

#For more information. . .
?attribute.methods

####################################
#2.6 Edgelists

edgelist<-read.csv("edgeList.csv")
edgeNet<-network(edgelist,matrix.type="edgelist")
edgeNet

plot(edgeNet,displaylabels=T) #what's missing?

plot(edgeNet,displaylabels=T,edge.lwd=edgeNet%e%"Weight"*10)

##how can you extract the weights?
edgeNet[,]
edgeNet%e%"Weight"
as.sociomatrix.sna(edgeNet,"Weight")
as.edgelist.sna(edgeNet)

##now in igraph
library(igraph)

#this will throw an error - why?
iedge<-graph_from_edgelist(edgelist,directed=T)

iedge<-graph_from_data_frame(edgelist,directed=T)
plot(iedge) #where are the weights?

E(iedge)$Weight
plot(iedge,edge.width=E(iedge)$Weight*5)
detach("package:igraph")

###########################################################
###########################################################
##3 Classical Network Analysis with the SNA package
##3.1 Getting started
library(sna) # Load the sna library
library(help="sna") # See a list of help pages for sna package
load("SNAworkshopData.Rdata") # Load supplemental workshop data


###########################################################
##3.2 Network visualization with gplot()
#Plotting the data we imported earlier
gplot(nrelations)
gplot(nrelations,displaylabels=TRUE) # This time with labels
#Let's color the nodes in sex-stereotypic colors
nodeColors<-ifelse(nodeInfo$sex=="F","hotpink","dodgerblue")
gplot(relations,gmode="graph",displaylabels=TRUE,vertex.col=nodeColors)

#Using data we just loaded in, plot the contiguity among nations in 1993 (from the Correlates of War (CoW)1 project)
gplot(contig_1993) # The default visualization
gplot(contig_1993, usearrows=FALSE) # Turn off arrows manually
gplot(contig_1993, gmode="graph") # Can also tell gplot the data is undirected

#We can add labels to the vertices
gplot(contig_1993, gmode="graph",displaylabels=TRUE,label.cex=0.5,label.col="blue")

#Other ways to specify the labeling
gplot(contig_1993,gmode="graph",label=network.vertex.names(contig_1993),label.cex=0.5,label.col="blue")

#Here's an example of directed data|militarized interstate disputes (MIDs) for 1993
gplot(mids_1993,label.cex=0.5,label.col="blue",displaylabels=TRUE)

#All those isolates can get in the way|we can suppress them using displayisolates
gplot(mids_1993,label.cex=0.5,label.col="blue",displaylabels=TRUE,displayisolates=FALSE)

#The default layout algorithm is that of Frutchterman-Reingold (1991), can use others
gplot(mids_1993,label.cex=0.5,label.col="blue",displaylabels=TRUE,displayisolates=FALSE,mode="circle") # The infamous circle

#or perhaps
gplot(mids_1993,label.cex=0.5,label.col="blue",displaylabels=TRUE,displayisolates=FALSE,mode="mds") # MDS of position similarity

#or perhaps
gplot(mids_1993,label.cex=0.5,label.col="blue",displaylabels=TRUE,displayisolates=FALSE,mode="target") # a target


#When a layout is generated, the results can be saved for later reuse:
coords <- gplot(contig_1993,gmode="graph",label=colnames(contig_1993[,]),label.cex=0.5,label.col="blue") # Capture the magic of the moment
coords # Show the vertex coordinates

#Saved (or a priori) layouts can be used via the coord argument
gplot(mids_1993,gmode="graph",label=colnames(contig_1993[,]),label.cex=0.5,label.col="blue",coord=coords)

#When the default settings are insuficient, interactive mode allows for tweaking
coords <- gplot(contig_1993, interactive=TRUE) # Modify and save
gplot(contig_1993,coord=coords,displaylabels=TRUE,gmode="graph",label.cex=0.5,label.col="blue") # Should reproduce the modified layout

#For more information. . .
?gplot ?gplot.layout

###########################################################
##3.3 Basic centrality indices: degree, betweenness, and closeness
#We begin with the simplest case: degree centrality

#Component information can be obtained in various ways
components(mids_1993) # Strong component count
components(mids_1993, connected="weak") # Weak component count
cd <- component.dist(mids_1993, connected="weak") # Get weak components
cd

#Component sizes
plot(1:length(cd$cdist),cd$cdist,xlab="Size",ylab="Frequency")

#Who's in the largest component?
cl <- component.largest(mids_1993, connected="weak")
cl

#Plot the largest weak component
gplot(mids_1993[cl,cl], 
      boxed.lab=FALSE, 
      label.cex=0.5,
      label.col=4,
      label=network.vertex.names(mids_1993)[cl])

#Likewise, many routines exist for handling isolates
is.isolate(mids_1993, 3) # Is the third vertex (BHM) an isolate?
is.isolate(mids_1993,which(mids_1993%v%"vertex.names"=="BHM")) # The peaceful islands?
is.isolate(mids_1993,which(mids_1993%v%"vertex.names"=="USA")) # Perhaps less so....
isol <- isolates(mids_1993) # Get the entire list of isolates
isol
network.vertex.names(mids_1993)[isol] # Which countries are isolates?

#Another way to remove isolates from sociograms
gplot(mids_1993[-isol,-isol], 
      label.cex=0.5,
      label.col=4,
      label=network.vertex.names(mids_1993)[-isol])

#Geodesic distances
contig.dist<-geodist(contig_1993)

#look at the resulting object
attributes(contig.dist)
contig.dist$gdist
contig.dist$counts


###########################################################
##3.2 Basic centrality indices: degree, betweenness, and closeness
#We begin with the simplest case: degree centrality
degree(mids_1993) # Default: total degree
ideg <- degree(mids_1993, cmode="indegree") # Indegree for MIDs
odeg <- degree(mids_1993, cmode="outdegree") # Outdegree for MIDs
all(degree(mids_1993) == ideg+odeg) # In + out = total?

#Once centrality scores are computed, we can handle them using standard R methods:
plot(ideg, 
     odeg, 
     type="n", 
     xlab="Incoming MIDs", 
     ylab="Outgoing MIDs") # Plot ideg by odeg

abline(0, 1, lty=3)

text(jitter(ideg), 
     jitter(odeg), 
     network.vertex.names(contig_1993), 
     cex=0.75, 
     col=2)

#Plot simple histograms of the degree distribution:
par(mfrow=c(2,2)) # Set up a 2x2 display

hist(ideg, 
     xlab="Indegree",
     main="Indegree Distribution", 
     prob=TRUE)

hist(odeg, 
     xlab="Outdegree", 
     main="Outdegree Distribution", 
     prob=TRUE)

hist(ideg+odeg, 
     xlab="Total Degree", 
     main="Total Degree Distribution", 
     prob=TRUE)

par(mfrow=c(1,1)) # Restore display

#Centrality scores can also be used with other sna routines, e.g., gplot()
gplot(mids_1993, 
      vertex.cex=(ideg+odeg)^0.5/2, 
      vertex.sides=50,
      label.cex=0.4,
      vertex.col=rgb(odeg/max(odeg),0,ideg/max(ideg)),
      displaylabels=TRUE)

#Betweenness and closeness are also popular measures
bet <- betweenness(contig_1993, 
                   gmode="graph") # Geographic betweenness

bet

gplot(contig_1993, 
      vertex.cex=sqrt(bet)/25, 
      gmode="graph") # Use w/gplot

clo <- closeness(contig_1993) # Geographic closeness
clo # A large world after all?

#Can use sna routines to explore alternatives to the common measures. . .
contig_1993_gdist<-geodist(contig_1993)$gdist # matrix of geodesic distances

contig_1993_gdist
1/sum(contig_1993_gdist[1,2:186]) # calculate closeness for country 1

which(contig_1993_gdist[,1]=="Inf") # disconnected matrix

sum(1/contig_1993_gdist[1,2:186]) # alternate closeness for country 1

closeness2 <- function(x){ # Create an alternate closeness function!
  geo <- 1/geodist(x)$gdist # Get the matrix of 1/geodesic distance
  diag(geo) <- 0 # Define self-ties as 0
  apply(geo, 1, sum) # Return sum(1/geodist) for each vertex
}

clo2 <- closeness2(contig_1993) # Use our new function on contiguity data
clo2

hist(clo2, xlab="Alt. Closeness", prob=TRUE) # Much better behaved!
cor(clo2, bet) # Correlate with betweenness
plot(clo2, bet) # Plot the bivariate relationship
all(clo2/185==closeness(contig_1993,cmode="suminvundir")) # Actually, we support this in sna!

#For more information. . .
?betweenness
?bonpow
?closeness
?degree
?evcent
?graphcent
?infocent
?prestige
?stresscent

###########################################################
##3.3 From centrality to centralization
centralization(mids_1993, degree, cmode="indegree") # Do MIDs concentrate?
centralization(contig_1993, evcent) # Eigenvector centralization

#For more information. . .
?centralization
###########################################################
##3.4 Elementary graph-level indices
gden(mids_1993) # Density
grecip(mids_1993) # Dyadic reciprocity
grecip(mids_1993, measure="edgewise") # Edgewise reciprocity
gtrans(mids_1993) # Transitivity

#For more information. . .
?gden 
?grecip 
?gtrans

###########################################################
##3.5 Subgraph census routines, isolates, clusters, and geodesic distance

dyad.census(mids_1993) # M,A,N counts
dyad.census(contig_1993) # No As in undirected graphs
triad.census(mids_1993) # Directed triad census
triad.census(contig_1993, mode="graph") # Undirected triad census

###########################################################
##3.6 Brokerage

#look at the trade data
class(trade)
dim(trade)

gplot(trade[1,,])

tradeat
table(tradeat[,1])

gplot(trade[1,,],vertex.col=rainbow(5)[tradeat[,1]+1]) #color by pop growth
tradeBro<-brokerage(trade[1,,],tradeat[,1])
summary(tradeBro)


###########################################################
##3.7 Basic connectivity/distance measurement, and cohesion
#Several ways to get relatively cohesive groups

kpath.census(mids_1993, maxlen=6, tabulate.by.vertex=FALSE) # Count paths of length <=6
kcycle.census(mids_1993, maxlen=6, tabulate.by.vertex=FALSE) # Count cycles of length <=6
clique.census(mids_1993, tabulate.by.vertex=FALSE, enumerate=FALSE) # Find maximal cliques

#Can also look at more detailed tabulation/co-membership for paths/cycles/cliques
kpath.census(mids_1993, maxlen=4) # Show tabulation by vertex

indirect <- kpath.census(mids_1993, 
                         maxlen=6, 
                         dyadic.tabulation="sum")$paths.bydyad

gplot(indirect,
      label.cex=0.4,
      vertex.cex=0.75,
      displaylabels=TRUE,
      edge.col=rgb(0,0,0,0.25)) # Plot indirect MIDs

#Showing cohesion information can aid visualization. Here, show critical positions
gplot(contig_1993,
      vertex.col=2+cutpoints(contig_1993,mode="graph",return.indicator=T))

#Show the nesting of cores
kc<-kcores(contig_1993,cmode="indegree")
gplot(contig_1993,vertex.col=rainbow(max(kc)+1)[kc+1])

#Showing members of the 5-core only
gplot(contig_1993[kc>4,kc>4],vertex.col=rainbow(max(kc)+1)[kc[kc>4]+1])

#For more information. . .
?bicomponent.dist
?cutpoints
?kcores
?is.connected
?reachability
?symmetrize

####################################################
#3.8 Blockmodeling
# Load the Faux Mesa High data set (aka "faux.mesa.high")
library(statnet)
data(faux.mesa.high)

# Display the data, along with covariates
plot(faux.mesa.high,vertex.col="Sex")        # Plot by gender (covariate is misnamed!)
vals<-sort(unique(faux.mesa.high%v%"Sex"))   # Let's add a legend....
legend("topleft",fill=1:length(vals),legend=vals,bty="n")

plot(faux.mesa.high,vertex.col="Grade")      # Now, plot by grade
legend("topleft",fill=7:12,legend=7:12,bty="n")

plot(faux.mesa.high,vertex.col="Race")       # Finally, plot by race
vals<-sort(unique(faux.mesa.high%v%"Race"))
legend("topleft",fill=1:length(vals),legend=vals,bty="n")

# For more information....
?faux.mesa.high
?legend
#
#Obtaining mixing matrices
# The mixingmatrix command is poorly documented, but very useful.  Let's
# obtain mixing matrices for each of the above attributes: 
mms<-mixingmatrix(faux.mesa.high,"Sex")               # Compute w/a network object and
mmg<-mixingmatrix(faux.mesa.high,"Grade")             # attribute name  
mmr<-mixingmatrix(faux.mesa.high,"Race")
mms                                             # See what we have...
mmg
mmr
names(mms)                                      # This is a complex object
class(mms)
mms$matrix                                      # Extract the mixing matrix

# Can crudely visualize using plot.sociomatrix:
plot.sociomatrix(mms$matrix)                    # Must use matrix, not object
plot.sociomatrix(mmg$matrix)
plot.sociomatrix(mmr$matrix)

########################
##3.9 Structural Equivalence
##use the Mt. St. Helens EMON (treating edge 
# reports as directed and unvalued, for demonstration purposes)
data(emon)
mtsth<-emon[[5]]                           # For convenience, pull the data
plot(mtsth)                                # Remind ourselves what it looks like
plot.sociomatrix(mtsth,drawlab=FALSE)

# We can visualize similarity/distances using multidimensional scaling
d<-sedist(mtsth,method="euclidean")        # Get Euclidean distances
coord<-cmdscale(d)                         # 2-D metric MDS
plot(coord)                          # Plot results
plot(coord,type="n")                 # Cool trick: plot w/labels
text(coord,label=network.vertex.names(mtsth),cex=0.5)


# The new version of sna also supports role similarity using CATREGE
d<-redist(mtsth)                           # Won't work with the old version...
d
coord<-cmdscale(d)
plot(coord)                          # Small number of distinct classes

# For more information....
?sedist
?redist
?cmdscale


###########################################################
#3.10 From equivalence to blockmodels
#
# To perform cluster analysis of network positions, we use the equiv.clust
# command (which is a front-end to both routines like sedist and to hclust)
ec<-equiv.clust(mtsth)          # Complete link SE clustering using Hamming dist
ec                              # Gives not so useful summary info
plot(ec)                        # Shows dendrogram
rect.hclust(ec$cluster,h=15)    # If we cut it at dist=15, what would happen?
rect.hclust(ec$cluster,k=3,border="blue")     # Split into 3 positions; what's the best cut?

# plotting is prettier in igraph:
library(intergraph)
library(igraph)
imtsth<-intergraph::asIgraph(mtsth)
groups<-cutree(ec$cluster,k=3)
listOfGroups<-split(1:27,groups)
plot(imtsth,mark.groups=listOfGroups)

detach("package:igraph")

# Complete link clustering is only one option; Ward's method often does well
# at carving up the MDS space....
ec<-equiv.clust(mtsth,cluster.method="ward")
plot(ec)

# Now that we have our SE clustering, let's try a blockmodel
bm<-blockmodel(mtsth,ec,h=15)               # Cut the tree at 15
bm                                          # Examine the blockmodel
# One way to plot the block image:
gplot(bm$block.model>gden(mtsth),diag=T,edge.lwd=bm$block.model*3) 
# Also useful:
o<-order(bm$block.membership)
plot.sociomatrix(mtsth[o,o],label=list(bm$block.membership[o],
                                       bm$block.membership[o]))

# Fit here seems pretty poor; how many of these blocks even make sense?
bm<-blockmodel(mtsth,ec,h=15,block.content="types")   # Classify the blocks
bm                                                    # Not really SE...



# The new sna version will also support redist (just like sedist, but for
# role similarity)
ec<-equiv.clust(mtsth,equiv.fun="redist",method="catrege") 
plot(ec)       

# For more information....
?equiv.clust
?hclust
?rect.hclust
?blockmodel

##There are a ton of community detection methods in igraph, but none that have been written into Statnet. These are some of the idiosyncratic differences.

library(igraph)
imtsth_undir<-as.undirected(imtsth,mode="collapse")
ceb<-cluster_edge_betweenness(imtsth_undir)
modularity(ceb)

plot(ceb,imtsth_undir)

cfg<-cluster_fast_greedy(imtsth_undir)
plot(cfg,imtsth_undir)
modularity(cfg)

detach("package:igraph")

###########################################################
##4.1 Elementary random graph generation

rgraph(10) # A uniform random digraph of order 10
rgraph(10, tprob=3/9) # Homogeneous Bernoulli w/mean degree 3
rgraph(10, tprob=3/9, mode="graph") # As above, but undirected
rgnm(1, 10, 20) # Uniform conditional on order, edges
rguman(1, 10, mut=0.5, asym=0.25, null=0.25) # Homogeneous multinomial on dyad census
rguman(1, 10, mut=0, asym=1, null=0) # An extreme case: random tournament
gplot(rgws(1,50,1,2,0)) # A Watts-Strogatz process - baseline
gplot(rgws(1,50,1,2,0.05)) # ...with rewiring probability 0.05
gplot(rgws(1,50,1,2,0.2)) # ...with rewiring probability 0.2

##now let's re-create the plot from Watts and Strogatz
##a list is a new type of data structure...
wgNets<-list()
temp<-seq(-3.2,0,.2)
temp
length(temp)

p_vals<-10^(temp)
p_vals

for(i in 1:17){
  wgNets[[i]]<-rgws(n=1,nv=50,d=1,z=2,p=p_vals[i])
}

##lists use 'double-bracket' notation
gplot(wgNets[[1]])
gplot(wgNets[[17]])

l_vals<-vector(length=17)
geodist(wgNets[[1]])$gdist
geodist(wgNets[[17]])$gdist
for(i in 1:17){
  temp<-geodist(wgNets[[i]])$gdist
  temp[which(temp=="Inf")]<-0
  l_vals[i]<-max(temp)
}
l_vals
l_vals<-l_vals/l_vals[1]

##now for the clustering coefficient
triad.census(wgNets[[1]],mode="graph")
c_vals<-vector(length=17)
for(i in 1:17){
  temp<-triad.census(wgNets[[i]],mode="graph")
  c_vals[i]<-temp[4]/(3*temp[4]+temp[3])
}
c_vals
c_vals<-c_vals/c_vals[1]

plot(p_vals,
     l_vals,log="x",
     ylim=c(0,1),
     xlab="Re-Wiring Level",
     ylab="Fraction of Initial Value")

lines(p_vals,c_vals,type="p",pch=2)

##always add legends!
legend("bottomleft",legend=c("Diameter","Clustering Coefficient"),pch=c(1,2))

#For more information. . .
?rgbn
?rgmn
?rgraph
?rguman
?rgws
?list
?lines
?legend

