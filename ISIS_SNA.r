library(igraph)
library(lubridate)
library(dplyr)
library('syuzhet')
require("ggplot2")
library(quanteda)
library(quanteda.textstats)
library(quanteda.textplots)
library(circlize)
  
df <- read.csv("tweets.csv") # import the data
head(df)  

# for the initial graph we are interested in connection between people, so we will create a directed graph
# where an edge exit from a node in case of a mention or retweet to another member
# so we need to create edge that goes from a user to another if the second user is cited with @*username* by the first 
usernames <- unique(df$username)

results <- list()
row_index <- 1
# get a list of pair: user that mention, user mentioned
for (i in 1:nrow(df)) { 
  from_user <- df$username[i]
  tweet_text <- df$tweets[i]
  
  mentions <- regmatches(tweet_text, gregexpr("@\\w+", tweet_text))[[1]] # get the string in the text that are composed by @*username* 
  
  if (length(mentions) > 0) {
    for (mention in mentions) {
      to_user <- sub("@", "", mention)  # remove the @ from the username
      if(to_user %in% usernames)  #  save only if the mentioned user is in our data set
      results[[row_index]] <- data.frame(from = from_user, to = to_user, stringsAsFactors = FALSE)
      row_index <- row_index + 1
    }
  }
}

mentions <- do.call(rbind, results) # transform the list of dataframe in a dataframe
head(mentions)
users <- unique(c(mentions$from, mentions$to)) # list of users mentioned in the edges (we need it later)
# now we have our list of mentions, let's transform it into a weighted edge list:

mention_counts <- as.data.frame(table(mentions$from, mentions$to), stringsAsFactors = FALSE)
colnames(mention_counts) <- c("from", "to", "weight")

mention_counts <- mention_counts[mention_counts$weight > 0, ]
head(mention_counts)

# now we have all we need to create our first graph

net <- graph_from_data_frame(d=mention_counts, vertices = usernames, directed=T)
net <- simplify(net, remove.loops=T)
V(net)$size <- 5
l <- layout.fruchterman.reingold(net)
E(net)$width <- ((E(net)$weight)^(1/4))/3

plot(net, edge.arrow.size=.1, edge.curved=.1, layout=l, vertex.label=NA)

# We can now get some information about the data:

# Node degree:
deg <- igraph::degree(net, mode="in") # Node degree -> most mentioned nodes on the social
plot(net, edge.arrow.size=.1, edge.curved=.1, layout=l, vertex.label=NA, vertex.size=deg/3+3)
hist(deg, breaks=1:vcount(net)-1, main="Histogram of Node Degree") 
deg.dist <- degree_distribution(net, cumulative=T, mode="in")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", xlab="Degree", ylab="Cumulative Frequency")
# we can see that most of the nodes are not mentioned, but can still mention other nodes, so let's also see the out-degree

deg <- igraph::degree(net, mode="out") # Node degree -> most active nodes on the social
plot(net, edge.arrow.size=.1, edge.curved=.1, layout=l, vertex.label=NA, vertex.size=deg/3+3)
hist(deg, breaks=1:vcount(net)-1, main="Histogram of Node Degree")

# It can be interesting to measure also reciprocity, so the proportion of reciprocated ties
reciprocity(net)
dyad_census(net) # we can say that out network is definitely sparse


# Let's get some more insight of the network, but first it's better to remove the nodes that do not appear in the edge table since they are not relevant for this part of the analysis
# we transform the graph net to undirected, in this operation it's important to sum the weights of edges (a,c) (c,a) to have a single edge per connection
und_net <- as.undirected(net, mode = "collapse", edge.attr.comb = list(weight = "sum"))
und_net <- simplify(und_net, remove.loops=T)
E(und_net)$width <- E(und_net)$weight/50
V(und_net)$size <- 5
net2 <- und_net - V(und_net)[igraph::degree(und_net, mode="all")==0]
layout <- layout.fruchterman.reingold(net2)
plot(net2, edge.arrow.size=.1, edge.curved=.1, layout=layout, vertex.label=NA) # in this way we also get a cleaner graph
plot(net2, edge.arrow.size=.1, edge.curved=.1, layout=layout_in_circle, vertex.label=NA)
plot(net2, edge.arrow.size=.1, edge.curved=.1, layout=layout_randomly, vertex.label=NA)


# since our network is basically a citation network, we can measure the centrality of the graph and see the most important vertex of the graph
# we will at first measure eigen centrality, that is a refinement that assigns higher weight to vertex for being connected to vertices which are themselves important
eigenCent <- eigen_centrality(net2)$vector
# let's plot the score on the graph:
bins <- unique(quantile(eigenCent, seq(0,1,length.out=15)))
vals <- cut(eigenCent, bins, labels=FALSE, include.lowest=TRUE)
my_col = heat.colors(length(bins))
colorVals <- rev(my_col)[vals]
V(net2)$color <- colorVals
plot(net2, edge.arrow.size=.1, edge.curved=.1, layout=layout, vertex.label=NA)
# so who are the more relevant user?
sort(eigenCent,decreasing=TRUE)[1:10]

# we will now measure betweenness centrality, to see which are the vertices that connect important parts of the graph
betweenCent <- betweenness(net2)
cor(betweenCent,eigenCent)
bins <- unique(quantile(betweenCent, seq(0,1,length.out=30)))
vals <- cut(betweenCent, bins, labels=FALSE, include.lowest=TRUE)
colorVals <- rev(heat.colors(length(bins)))[vals]
V(net2)$color <- colorVals
plot(net2, edge.arrow.size=.1, edge.curved=.1, layout=layout, vertex.label=NA)
# Which are the most useful users to the graph connection?
sort(betweenCent,decreasing=TRUE)[1:10]


# we now want to do a temporal analysis on the connections, specifically we want to check for triadic closures:
# if A is friend with B and B is friend with C, then A tend to become friend with C
# At first we want to calculate the number of triangles that will create over time and the probability of triadic closure, how many open triads become close triads.

# first we need to sort the data set based on the timestamp:
class(df$time) # we want a datatime value not character
df$time <- mdy_hm(df$time)
class(df$time)
sorted_df <- df[order(df$time),]

# now we want to create the graph progressively by the time, the idea is, since the graph is quite sparse, to see how it evolve month by month
df$month <- format(df$time, "%Y-%m") # to do this we need a new column -> month
net3 <- make_empty_graph(directed=F) # start from an empty graph
net3 <- add_vertices(net3, length(users), name=users)
layout2 <- layout_on_sphere(net3)
#layout <- layout_in_circle(net)
# and for each month, we add the edges for the monthly tweets, measure the number of triadic closure and update the graph
triangles_over_time <- c()
closure_prob_ot <- c()
opt <- c()
graph_list <- list()
months <- sort(unique(df$month))
months <- months[-c(1,2,3,4,5,6,7,8)]
months
par(mfrow=c(3,3))

for (m in months){
  df_month <- df[df$month == m, ]
  
  # this part is the same as before to create the edges:
  edges <- c()
  for (i in 1:nrow(df_month)) { 
    from_user <- df_month$username[i]
    tweet_text <- df_month$tweets[i]
    
    mentions <- regmatches(tweet_text, gregexpr("@\\w+", tweet_text))[[1]] # get the string in the text that are composed by @*username* 
    
    if (length(mentions) > 0) {
      for (mention in mentions) {
        to_user <- sub("@", "", mention)  # remove the @ from the username
        if(to_user %in% usernames)  #  save only if the mentioned user is in our data set
          edges <- c(edges, from_user, to_user)
      }
    }
  }
  # create the new edges
  if(length(edges) > 0){
    new_edges <- matrix(edges, ncol=2, byrow=T)
    for (i in 1:nrow(new_edges)){
      # we add the edge only if it is not a loop(mention to itself), we do this because this type of connection are not useful for our research. 
      # we add an edge only if it wasn't already in the graph, we are interested in the single connection not in weights
      if(new_edges[i,1] != new_edges[i,2] && !are.connected(net3, new_edges[i,1], new_edges[i,2]))
        net3 <- add_edges(net3, c(new_edges[i,1], new_edges[i,2]))
        
    }
  }
  # calculate the triangles and probability of closure
  triangles <- sum(count_triangles(net3))/3
  triangles_over_time <- c(triangles_over_time, triangles)
  # calculate the number of open triangles -> we have 2 of the 3 edges already and we want to close it so add the remaining edge
  deg <- igraph::degree(net3)
  open_triangles <- sum(deg * (deg-1) / 2) # sum of possible combination for each vertex
  if(open_triangles > 0)
    closure_prob <- sum(count_triangles(net3)) / open_triangles # the number of closure triadic compared to the number of possible triadic
  else 
    closure_prob <- 0
  opt <- c(opt, open_triangles)
  closure_prob_ot <- c(closure_prob_ot, closure_prob)
  # save the snapshot of the graphs
  graph_list[[m]] <- as.undirected(net3)
  # plot the graph of the month
  plot(
    net3, 
    main=paste("Rete al mese:", m),
    layout=layout2,
    vertex.label=NA,
    edge.arrow.size=.2,
    vertex.size=deg/10+2
  )
}

par(mfrow = c(1, 1))
plot(
  net3, 
  main=paste("Rete al mese:", m),
  layout=layout2,
  vertex.label=NA,
  edge.arrow.size=.2,
  vertex.size=deg/10+2
)
triangles_over_time
opt
closure_prob_ot
nmonths <- factor(months, levels = c("2015-09", "2015-10", "2015-11", "2015-12", "2016-01", "2016-02", "2016-03", "2016-04", "2016-05"), ordered=T)
par(mfrow=c(1,2))
plot(nmonths, triangles_over_time, type = "a", col = "blue",
     xlab = "Month", ylab = "Triadic Closure", main = "Evolution of Triadic Closure")
plot(nmonths, triangles_over_time, type = "a", col = "blue",
     xlab = "Month", ylab = "Closure Probability", main = "Evolution of Closure Probability")



# How much more likely is a link to form between two people in a social network if they already have a connection in common?
# we want now to check if a connection is created more probably if there is already a common connection -> so if the network evolve more thanks to triadic closure or new casual connection.
# to do this we will track link formation, from a snapshot to another we will:
# - identify the pairs of nodes that have k connection in common in the first snapshot but are not directly connected by an edge
# - find T(k): the fraction of these pairs that have formed an edge by the time of the second snapshot
# we will then plot T as a function to illustrate the effect of common friends on the formation of links

T_k <- matrix(0, nrow = length(months)-1, ncol = 2) # when k is 5 and above in just one column
T_k_frac <- matrix(0, nrow = length(months)-1, ncol = 2)
pairs <- matrix(0, nrow = length(months)-1, ncol=2)
for (m in 1:(length(months)-1)){
  snap1 <- graph_list[[months[m]]]
  snap2 <- graph_list[[months[m+1]]]
  adj1 <- as_adj(snap1, sparse = FALSE)
  adj2 <- as_adj(snap2, sparse=FALSE)
  # first we get the common connections between each pair of nodes 
  com_connection <- adj1 %*% adj1 # this give us a matrix containing in each cell the number of path of length 2 between the nodes of row i and column j, so basically the number of common neighbors between node i and node j
  # we now have to obtain only the pairs that haven't already formed an edge, so we are interested only in pair that aren't connected in adj1 and check if they have formed a connection in adj2, we use com_connection to know how many common connection they have so to save them into the right column
  n <- nrow(adj1)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      if(adj1[i,j]==0){  # if they don't have a connection in the first snaphshot
        k <- com_connection[i,j]
        if(k>1) # we group all the one with above 5 connection together
          k <- 1
        
        pairs[m,k+1] <- pairs[m,k+1] + 1
        if(adj2[i,j]==1)
          T_k[m, k+1] <- T_k[m, k+1] + 1 # increase the function
      }
    }
  }
  T_k_frac[m,] <- T_k[m,] / pairs[m,]
}
# at this point we have a matrix containing the number of connection created every month where there is a common neighbor or not, the number of pair that are not connected by an edge that have a common neighbor or not and the fraction T_k for each month
barplot(t(T_k), names.arg=months[-1], col=c("red", "blue"), legend.text = c("0 common neighbor", "≥1 common neighbor"))
# since we have that in the first few months less connection were made, but the number of possible connection was quite the same, it makes more sense not to do the mean of the months but to sum up the connection made and the not connected pair for each month and then at that point look at the total fraction.
# I decided to do this because low activity months can distort the average if we take the mean, but if we sum numerator and denominators across all months we get a more robust estimator that reflect the actual volume of link formation.
T_k
pairs
T_k_tot <- colSums(T_k)
# since we are calculating the total opportunity over time, we will sum up the denominator, this is because in every month the pairs had the possibility to connect
# in this way we respond to the question: out of all the times a connection could have formed, how often did it actually form across the months?
pairs_tot <- colSums(pairs)
T_k_tot # total number of connection created where there is at least a common neighbor or not
pairs_tot # total number of pairs that can create a connection during the months, so each pair that did not create a connection is retake for the sum each month
T_k_tot <- T_k_tot / pairs_tot
T_k_tot

barplot(T_k_tot, beside=T, names.arg=c("0 common neighbors", "≥1 common neighbors"), main = "T(k) per month",
        ylab = "Probability of nodes connection", col=c("red", "blue"))

# we can see that indeed the triadic closure means a lot in the formation of a connection: when 2 nodes had 0 connection in common at the beginning of a month, the probabilty for them to create a connection in that month was on average 0.007 meanwhile when they had 1 or more connection in common, the probability is doubled



# we are now interested in the most used words in the tweets of the ISIS fan, so from the starting dataset, we are interested in the fields tweet

tweets <- df[, c("username", "tweets")]
# we are interested only on the message in the tweet, not to link or mentions
tweets$tweets <- gsub("ENGLISH TRANSLATION:|rt", "", tweets$tweets, ignore.case = TRUE)
tweets$tweets <- gsub("http[^[:space:]]+", "", tweets$tweets)
tweets$tweets <- gsub("@\\w+", "", tweets$tweets)
tweets$tweets <- gsub("#\\w+", "", tweets$tweets)

head(tweets)

corpus = corpus(tweets, text_field = "tweets")
summary(corpus)
doc.tokens = tokens(corpus) # tokenize the text (split each document into individual tokens)
doc.tokens = tokens(doc.tokens, remove_punct = TRUE, remove_numbers = TRUE) # we remove punctuation and numbers from the token (they are noise)
doc.tokens = tokens_select(doc.tokens, stopwords(language = "en", source = "snowball", simplify = TRUE), selection ='remove') # we remove common english words like "is", "the", "and"
doc.tokens = tokens_tolower(doc.tokens) # convert all words to lower case, making analysis consistent
doc.tokens <- tokens_keep(doc.tokens, pattern = "^[a-z]+$", valuetype = "regex") # we keep only tokens that are entirely lowercase alphabet characters

toks_ngram = tokens_ngrams(doc.tokens, n = 2) # non only single words but also pairs of consecutive words, we now have a richer set of features than just individual words
toks_ngram

dfmat = dfm(toks_ngram) %>% dfm_trim(min_termfreq = 20) # this create a DFM (rows are the documents, columns are the features -> cells contains frequency of feature in a document) and we filter out rare terms (we only want terms with frequency > 10)
dfmat


features_dfm = textstat_frequency(dfmat, n = 100)
features_dfm$feature = with(features_dfm, reorder(feature, -frequency))
features_dfm

textplot_wordcloud(dfmat)

# we now want to check the most used hashtags
hashtags <- df[, c("username", "tweets")]
# this time we are interested only in the hashtags
hashtags$hashtag_list <- regmatches(hashtags$tweets, gregexpr("#\\w+", hashtags$tweets))
hashtags$hashtag <- sapply(hashtags$hashtag_list, paste, collapse=" ")
head(hashtags$hashtag_list)
corpus = corpus(hashtags, text_field = "hashtag")
summary(corpus)
toks_hashtag = tokens(corpus) # tokenize the text (split each document into individual tokens)
toks_hashtag = tokens(toks_hashtag, remove_punct = TRUE) # we remove punctuation and numbers from the token (they are noise)
dfmat2 = dfm(toks_hashtag) %>% dfm_trim(min_termfreq = 20)
dfmat2
features_dfm = textstat_frequency(dfmat2, n = 100)
features_dfm$feature = with(features_dfm, reorder(feature, -frequency))
features_dfm
textplot_wordcloud(dfmat2)


# Sentiment Analysis

# use NRC Emotion Lexicon (list of words and their associations)
sentiment = get_nrc_sentiment(tweets$tweets)
td = data.frame(t(sentiment))
td = data.frame(rowSums(td[-1]))
names(td)[1] <- "count"
tdw <- cbind("sentiment" = rownames(td), td)
rownames(tdw) <- NULL
tdw


# Plot Emotions
ggplot(tdw[1:8, ], aes(x = sentiment, y = count, fill = sentiment)) +
  geom_bar(stat = "identity") +
  labs(x = "emotion") +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.title = element_blank())

# Plot Polarity
ggplot(tdw[9:10, ], aes(x = sentiment, y = count, fill = sentiment)) +
  geom_bar(stat = "identity") +
  labs(x = "polarity") +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.title = element_blank())






# Now that we have complete the analysis of the complete graph let's try to go a little more in the specific and do an analysis of the communities
# let's now try to detect the communities -> disconnected parts of the graph

# thanks to community, it's possible to do an analysis of the content for each community or check the community over time to see if they change or not


communities <- multilevel.community(net2, weights = E(net2)$weight) # create community based on the strength of ties, in this way we group strong mutual connections more tightly
communities
V(net2)$community <- membership(communities)
V(net2)$community
inv_weights <- 1 / E(net2)$weight

centrality_df <- data.frame(
  username = V(net2)$name,
  community = V(net2)$community,
  degree = degree(net2),
  strength = strength(net2, weights = E(net2)$weight),
  centrality = eigen_centrality(net2, weights = E(net2)$weight)$vector,  # how central a community is
  betweenness = betweenness(net2, weights = inv_weights), # how much a node control the flow of information
  closeness = closeness(net2, weights = inv_weights)    # how fast a node can reach other nodes, the inverse of distances
)
centrality_df

# let's see the communities in a graph with a color for each of them:
colorVals <- c("firebrick3", "slateblue1", "yellow1", "olivedrab1", "pink", "seagreen1", "orange", "turquoise1")
V(net2)$color <- colorVals[V(net2)$community]
plot.igraph(net2, vertex.label = NA,layout=layout, vertex.size=5)

# we can now check how the communities behave in the network: the community that are more influent
aggregate(. ~ community, data = centrality_df[-1], FUN = mean)
# we can see that the most important communities under every aspect are community 1, 2 and 3

# let's see more in the detail the structure of the communities, we would like to know the density of the community, the triadic closure and the assortivity (if nodes tends to connect to similiar nodes)
community_id <- unique(V(net2)$community)
edge_density(net2, loops=F)
density <- c()
nnodes <- c()
diameter <- c()
transitivity <- c()
for (i in community_id){
  sub_net <- induced.subgraph(net2, V(net2)[community==i])
  density <- c(density, edge_density(sub_net, loops=F))
  nnodes <- c(nnodes, vcount(sub_net))
  diameter <- c(diameter, length(get_diameter(sub_net)))
  transitivity <- c(transitivity, transitivity(sub_net))
}
density # we can see that the density inside each community is significantly higher than in the general network, we can see that the communities 1-4 have a lower density but it's normal since the have more nodes
nnodes # indeed having 20 nodes and a density of 0.15 means that the network is sparse but not a lot, indeed the diameter of the communities is around 5-6, we can say that the communities are sparse but definitely less than the network
diameter
transitivity
# we are also curious about the assortativity per community : how much does nodes tend to connect to nodes that are of the same community
assortativity_nominal(net2, types = as.factor(V(net2)$community), directed=F)
# the value is positive, so we can say that there is the tendency of nodes to connect with other nodes of the same community
# we can say that nodes interact more with nodes of the same community than with other nodes -> mentions are more common inside a community, probably information circles more inside each community, but to be sure of this let's see if the text of each community are distinct or similiar:

head(tweets)
# we already have a data frame containing the messagges and the users, we need to keep only the users that are in a community and specify the community of each user
comtweets <- tweets %>% inner_join(centrality_df %>% select(username, community))
head(comtweets)
# now we are ready to look at the top 3 communities contents and see if there are differences, we just need to do the classic text and sentiment analysis:
community_ngrams <- list()
dfmat_coms <- list()
community_id <- unique(comtweets$community)
for (i in community_id){
  corpus = corpus(comtweets[comtweets$community==i,], text_field = "tweets")
  summary(corpus)
  doc.tokens = tokens(corpus) # tokenize the text (split each document into individual tokens)
  doc.tokens = tokens(doc.tokens, remove_punct = TRUE, remove_numbers = TRUE) # we remove punctuation and numbers from the token (they are noise)
  doc.tokens = tokens_select(doc.tokens, stopwords(language = "en", source = "snowball", simplify = TRUE), selection ='remove') # we remove common english words like "is", "the", "and"
  doc.tokens = tokens_tolower(doc.tokens) # convert all words to lower case, making analysis consistent
  doc.tokens <- tokens_keep(doc.tokens, pattern = "^[a-z]+$", valuetype = "regex") # we keep only tokens that are entirely lowercase alphabet characters
  
  community_ngrams[[i]] <- tokens_ngrams(doc.tokens, n = 2) # non only single words but also pairs of consecutive words, we now have a richer set of features than just individual words
  dfmat_coms[[i]] <- dfm(community_ngrams[[i]]) %>% dfm_trim(min_termfreq = 20)
}
textplot_wordcloud(dfmat_coms[[1]])

# use NRC Emotion Lexicon (list of words and their associations)
tdw_coms <- list()
for (i in community_id){
  sentiment = get_nrc_sentiment(comtweets[comtweets$community==i,]$tweets)
  td = data.frame(t(sentiment))
  td = data.frame(rowSums(td[-1]))
  names(td)[1] <- "count"
  tdw <- cbind("sentiment" = rownames(td), td)
  rownames(tdw) <- NULL
  tdw_coms[[i]] <-tdw
}


# Plot Emotions
ggplot(tdw[1:8, ], aes(x = sentiment, y = count, fill = sentiment)) +
  geom_bar(stat = "identity") +
  labs(x = "emotion") +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.title = element_blank())


# let's now see who are the users that are more important in each of these community:

top_users <- centrality_df %>%
  filter(community %in% c(1, 2, 3)) %>%
  group_by(community) %>%
  top_n(2, centrality) %>% 
  arrange(community, desc(centrality))
top_users
users_id <- top_users$username

# let's see if what they say and the emotions they transmit are the same of their community or if they differ from it, this is also a way to see if there is some type of homophily and social influence in the communities:

users_ngrams <- list()
dfmat_users <- list()
for (i in users_id){
  corpus = corpus(comtweets[comtweets$username==i,], text_field = "tweets")
  summary(corpus)
  doc.tokens = tokens(corpus) # tokenize the text (split each document into individual tokens)
  doc.tokens = tokens(doc.tokens, remove_punct = TRUE, remove_numbers = TRUE) # we remove punctuation and numbers from the token (they are noise)
  doc.tokens = tokens_select(doc.tokens, stopwords(language = "en", source = "snowball", simplify = TRUE), selection ='remove') # we remove common english words like "is", "the", "and"
  doc.tokens = tokens_tolower(doc.tokens) # convert all words to lower case, making analysis consistent
  doc.tokens <- tokens_keep(doc.tokens, pattern = "^[a-z]+$", valuetype = "regex") # we keep only tokens that are entirely lowercase alphabet characters
  
  users_ngrams[[i]] <- tokens_ngrams(doc.tokens, n = 2) # non only single words but also pairs of consecutive words, we now have a richer set of features than just individual words
  dfmat_users[[i]] <- dfm(users_ngram[[i]]) %>% dfm_trim(min_termfreq = 20)
}

textplot_wordcloud(dfmat_users[1])

# use NRC Emotion Lexicon (list of words and their associations)
tdw_users <- list()
for (i in users_id){
  sentiment = get_nrc_sentiment(comtweets[comtweets$username==i,]$tweets)
  td = data.frame(t(sentiment))
  td = data.frame(rowSums(td[-1]))
  names(td)[1] <- "count"
  tdw <- cbind("sentiment" = rownames(td), td)
  rownames(tdw) <- NULL
  tdw_users[[i]] <- tdw
}

# Plot Emotions
ggplot(tdw[1:8, ], aes(x = sentiment, y = count, fill = sentiment)) +
  geom_bar(stat = "identity") +
  labs(x = "emotion") +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.title = element_blank())

