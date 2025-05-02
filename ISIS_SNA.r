library(igraph)
library(lubridate)
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

net <- graph_from_data_frame(d=mention_counts, vertices = usernames, directed=F)
net <- simplify(net, remove.loops=T)
V(net)$size <- 5
l <- layout.fruchterman.reingold(net)
E(net)$width <- E(net)$weight/50
plot(net, edge.arrow.size=.1, edge.curved=.1, layout=l, vertex.label=NA)

# We can now get some information about the data:

# Node degree:
deg <- degree(net, mode="in") # Node degree -> most mentioned nodes on the social
plot(net, edge.arrow.size=.1, edge.curved=.1, layout=l, vertex.label=NA, vertex.size=deg/3+3)
hist(deg, breaks=1:vcount(net)-1, main="Histogram of Node Degree") 
deg.dist <- degree_distribution(net, cumulative=T, mode="in")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", xlab="Degree", ylab="Cumulative Frequency")
# we can see that most of the nodes are not mentioned, but can still mention other nodes, so let's also see the out-degree

deg <- degree(net, mode="out") # Node degree -> most active nodes on the social
plot(net, edge.arrow.size=.1, edge.curved=.1, layout=l, vertex.label=NA, vertex.size=deg/3+3)
hist(deg, breaks=1:vcount(net)-1, main="Histogram of Node Degree")

# It can be interesting to measure also reciprocity, so the proportion of reciprocated ties
reciprocity(net)
dyad_census(net) # we can say that out network is definitely sparse


# Let's get some more insight of the network, but first it's better to remove the nodes that do not appear in the edge table since they are not relevant for this part of the analysis
net2 <- net - V(net)[degree(net, mode="all")==0]
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

# let's now try to detect the communities -> disconnected parts of the graph
w <- cluster_edge_betweenness(net2)
sort(table(w$membership)) # vector that associate each node with a community that represent the nodes that interact with each other mostly.
# we have 5 significant communities with 9 or more nodes, let's remake the graph with different color for each community -> graphical visualization:
V(net2)$color <- rep("white", length(w$membership))
keepTheseCommunities <- names(sizes(w))[sizes(w) > 4]
matchIndex <- match(w$membership, keepTheseCommunities) # like %in%
colorVals <- rainbow(10)[matchIndex[!is.na(matchIndex)]]
V(net2)$color[!is.na(matchIndex)] <- colorVals
plot.igraph(net2, vertex.label = NA,layout=layout, vertex.size=5)
# thanks to community, it's possible to do an analysis of the content for each community or check the community over time to see if they change or not


# we now want to do a temporal analysis on the connections, specifically we want to check for triadic closures:
# How much more likely is a link to form between two people in a social network if they already have a connection in common?
# if A is friend with B and B is friend with C, then A tend to become friend with C
# we want to calculate the number of triangles that will create over time, how many open triads become close triads and if a connection is created more probably if there is already a common connection -> so if the network evolve more thanks to triadic closure or new casual connection

# first we need to sort the data set based on the timestamp:
class(df$time) # we want a datatime value not character
df$time <- mdy_hm(df$time)
class(df$time)
sorted_df <- df[order(df$time),]

# now we want to create the graph progressively by the time, the idea is, since the graph is quite sparse, to see how it evolve month by month
df$month <- format(df$time, "%Y-%m") # to do this we need a new column -> month
net3 <- make_empty_graph(directed=F) # start from an empty graph
net3 <- add_vertices(net3, length(users), name=users)
#layout <- layout_in_circle(net)
# and for each month, we add the edges for the monthly tweets, measure the number of triadic closure and update the graph
triangles_over_time <- c()
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
  if(length(edges) > 0){
    new_edges <- matrix(edges, ncol=2, byrow=T)
    #new_nodes <- setdiff(unique(c(new_edges)), V(net)$name)
    #if (length(new_nodes) > 0) {
    #  net <- add_vertices(net, length(new_nodes), name=new_nodes)
    #}
    net3 <- add_edges(net3, t(new_edges))
  }
  triangles <- sum(count_triangles(net3))/3
  triangles_over_time <- c(triangles_over_time, triangles)
  net3 <- simplify(net3, remove.multiple = TRUE, remove.loops=T)
  plot(
    net3, 
    main=paste("Rete al mese:", m),
    layout=layout,
    vertex.size=3, 
    vertex.label=NA,
    edge.arrow.size=.2
  )
}

par(mfrow = c(1, 1))

plot(months, triangles_over_time, type = "b", col = "blue",
     xlab = "Mese", ylab = "Triangoli (Triadic Closure)", main = "Evoluzione del Triadic Closure")


# we are now interested in the most used words in the tweets of the ISIS fan, so from the starting dataset, we are interested in the fields tweet
library(quanteda)

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

library(quanteda.textstats)
features_dfm = textstat_frequency(dfmat, n = 100)
features_dfm$feature = with(features_dfm, reorder(feature, -frequency))
features_dfm

library(quanteda.textplots)
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

library('syuzhet')  # use NRC Emotion Lexicon (list of words and their associations)
sentiment = get_nrc_sentiment(tweets$tweets)
td = data.frame(t(sentiment))
td = data.frame(rowSums(td[-1]))
names(td)[1] <- "count"
tdw <- cbind("sentiment" = rownames(td), td)
rownames(tdw) <- NULL
tdw

require("ggplot2")
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



