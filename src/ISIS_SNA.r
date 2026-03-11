# ISIS SOCIAL NETWORK ANALYSIS

# 1. importo le librerie necessarie per il progetto
library(igraph)
library(lubridate)
library(dplyr)
library('syuzhet')
require("ggplot2")
library(quanteda)
library(tidyr)
library(quanteda.textstats)
library(quanteda.textplots)
library(gridExtra)

# 2. importo il dataset da un file CSV
df <- read.csv("../data/tweets.csv")
head(df)  

# ******************************************************************************

## CREAZIONE DELLA RETE

# 3. i vertici del grafo sono gli utenti
usernames <- unique(df$username)

results <- list()
row_index <- 1
# ciclo for utilizzato per ottenere gli archi della rete che corrispondono alle menzioni
for (i in 1:nrow(df)) { 
  # l'utente che menziona
  from_user <- df$username[i]
  tweet_text <- df$tweets[i]
  # vettore di utenti menzioni presenti in un messaggio, li ottengo selezionando le stringhe composte da @*username*
  mentions <- regmatches(tweet_text, gregexpr("@\\w+", tweet_text))[[1]] # get the string in the text that are composed by @*username* 
  
  if (length(mentions) > 0) {
    for (mention in mentions) {
      to_user <- sub("@", "", mention)  # rimuovo la @ dall'username
      if(to_user %in% usernames)  #  salvo solo gli utenti menzionati che sono inclusi nel dataset
      results[[row_index]] <- data.frame(from = from_user, to = to_user, stringsAsFactors = FALSE)
      row_index <- row_index + 1
    }
  }
}
# 4. trasformo la lista di dataframe in un dataframe singolo
mentions <- do.call(rbind, results) 
head(mentions)
users <- unique(c(mentions$from, mentions$to)) # lista di utenti menzionati negli archi (servirà per l'analisi temporale)

# conto le coppie: utente che menziona, utente menzionato per rendere gli archi pesati
mention_counts <- as.data.frame(table(mentions$from, mentions$to), stringsAsFactors = FALSE)
colnames(mention_counts) <- c("from", "to", "weight")

mention_counts <- mention_counts[mention_counts$weight > 0, ] # rimuovo le coppie con peso 0
head(mention_counts)

# 5. ora abbiamo tutto il necessario per la creazione del grafo diretto (il primo grafo)

net <- graph_from_data_frame(d=mention_counts, vertices = usernames, directed=T)
net <- simplify(net, remove.loops=T)
V(net)$size <- 5
l <- layout.fruchterman.reingold(net)
E(net)$width <- ((E(net)$weight)^(2/3))/3

plot(net, edge.arrow.size=.1, edge.curved=.1, layout=l, vertex.label=NA)
title("Grafo diretto completo")

# ******************************************************************************

# ANALISI STRUTTURALE DELLA RETE
# 6. 
vcount(net) # numero di vertici
ecount(net) # numero di archi

# densità dei nodi
edge_density(net)
# numero di nodi connessi mutuali e singole
dyad_census(net)
# la proporzione di connessioni mutuali in confronto alle singole
reciprocity(net)

# Grado dei nodi:

# 7. Incoming -> i nodi più menzionati nella rete sociale
par(mfrow=c(1,2))
indeg <- igraph::degree(net, mode="in") 
plot(net, edge.arrow.size=.1, edge.curved=.1, layout=l, vertex.label=NA, vertex.size=indeg/3+3)
hist(indeg, breaks=1:vcount(net)-1, main="Distribuzione Grado in Entrata",xlab = "Grado in entrata", ylab = "Frequenza") 
# I nodi con grado entrante maggiore
sorted_indeg <- sort(indeg, decreasing = TRUE)
head(sorted_indeg, 10)

# 8. Outgoing -> i nodi più attivi nella rete sociale
outdeg <- igraph::degree(net, mode="out") 
plot(net, edge.arrow.size=.1, edge.curved=.1, layout=l, vertex.label=NA, vertex.size=outdeg/3+3)
hist(outdeg, breaks=1:vcount(net)-1, main="Distribuzione Grado in Uscita",xlab = "Grado in uscita", ylab = "Frequenza") 
# I nodi con grado uscente maggiore
sorted_outdeg <- sort(outdeg, decreasing = TRUE)
head(sorted_outdeg, 10)

# ******************************************************************************

# CENTRALITÀ DEI NODI

# 9. creazione grafo indiretto senza utenti isolati (secondo grafo)
und_net <- as.undirected(net, mode = "collapse", edge.attr.comb = list(weight = "sum")) # l'ultimo parametro serve per sommare i pesi di (a,c) e (c,a) dato che il grafo è indiretto  
und_net <- simplify(und_net, remove.loops=T)
E(und_net)$width <- E(und_net)$weight/50
V(und_net)$size <- 5
net2 <- und_net - V(und_net)[igraph::degree(und_net, mode="all")==0] # rimuovo i nodi isolati
layout <- layout.fruchterman.reingold(net2)
par(mfrow=c(1,1))
plot(net2, edge.arrow.size=.1, edge.curved=.1, layout=layout, vertex.label=NA)

# 10. Eigenvector Centrality
eigenCent <- eigen_centrality(net2)$vector
# visualiziamo i punteggi dei nodi in un grafo dove un colore più scuro implica un punteggio più elevato
bins <- unique(quantile(eigenCent, seq(0,1,length.out=15)))
vals <- cut(eigenCent, bins, labels=FALSE, include.lowest=TRUE)
my_col = heat.colors(length(bins))
colorVals <- rev(my_col)[vals]
V(net2)$color <- colorVals
plot(net2, edge.arrow.size=.1, edge.curved=.1, layout=layout, vertex.label=NA)
# utenti con punteggio di eigenvector centrality più alto
sort(eigenCent,decreasing=TRUE)[1:10]

# 11. Betweenness Centrality
betweenCent <- betweenness(net2)
# correlazione tra i 2 tipi di centralità
cor(betweenCent,eigenCent)
# visualiziamo i punteggi dei nodi in un grafo dove un colore più scuro implica un punteggio più elevato
bins <- unique(quantile(betweenCent, seq(0,1,length.out=30)))
vals <- cut(betweenCent, bins, labels=FALSE, include.lowest=TRUE)
colorVals <- rev(heat.colors(length(bins)))[vals]
V(net2)$color <- colorVals
plot(net2, edge.arrow.size=.1, edge.curved=.1, layout=layout, vertex.label=NA)
# utenti con punteggio di betweenness centrality più alto
sort(betweenCent,decreasing=TRUE)[1:10]

# ******************************************************************************

# ANALISI TEMPORALE

# 12. trasformo il tipo della variabile da carattere a data, in questo modo posso ordinare il dataframe per dataclass(df$time)
df$time <- mdy_hm(df$time)
class(df$time)
sorted_df <- df[order(df$time),] # ordino il dataset per datetime di pubblicazione del tweet

# 13. essendo che misuriamo come cambia la rete mese dopo mese, abbiamo bisogno di una variabile month nel dataframe
df$month <- format(df$time, "%Y-%m")
# net3 è la variabile che rappresenta l'evoluzione del grafo mese per mese
net3 <- make_empty_graph(directed=F) # inizia da un grafo vuoto composto solo dai vertici
net3 <- add_vertices(net3, length(users), name=users)
layout2 <- layout_on_sphere(net3) 
# variabili utilizzate per misurare il numero di chiusure triadiche, probabilità di chiusura e il numero di possibili triangoli (2 lati su 3 presenti)
triangles_over_time <- c()
closure_prob_ot <- c()
opt <- c()
# lista dei grafi mese per mese
graph_list <- list()
months <- sort(unique(df$month))
# considero solo gli ultimi mesi in quanto il dataset diventa significativo dopo l'attentato di parigi
months <- months[-c(1,2,3,4,5,6,7,8)]
months
par(mfrow=c(3,3))

# 14. ciclo di computazione usato per ottenere : probabilità di chiusura triadica e visualizzazione della rete mese per mese
for (m in months){
  # prendo in considerazione soltanto i tweet del mese m
  df_month <- df[df$month == m, ]
  
  # stessa funzione utilizzata per la creazione delle menzioni ad inizio file
  edges <- c()
  for (i in 1:nrow(df_month)) { 
    from_user <- df_month$username[i]
    tweet_text <- df_month$tweets[i]
    
    mentions <- regmatches(tweet_text, gregexpr("@\\w+", tweet_text))[[1]]
    
    if (length(mentions) > 0) {
      for (mention in mentions) {
        to_user <- sub("@", "", mention) 
        if(to_user %in% usernames)
          edges <- c(edges, from_user, to_user)
      }
    }
  }
  # creo i nuovi archi
  if(length(edges) > 0){
    new_edges <- matrix(edges, ncol=2, byrow=T) # matrice che corrisponde a un vettore di coppie
    for (i in 1:nrow(new_edges)){
      # aggiungo l'arco al grafo solo se:
      # - non è una menzione a se stesso (un loop) dato che non è utile per la ricerca
      # - l'arco non è già presente nel grafo, non siamo interessati al peso in questo caso
      if(new_edges[i,1] != new_edges[i,2] && !are.connected(net3, new_edges[i,1], new_edges[i,2]))
        net3 <- add_edges(net3, c(new_edges[i,1], new_edges[i,2]))
        
    }
  }
  # calcolo il numero di triangoli e la probabilità di chiusura:
  triangles <- sum(count_triangles(net3))/3 # numero di triangoli (diviso 3 dato che vengono contati per ogni nodo del triangolo)
  triangles_over_time <- c(triangles_over_time, triangles)
  # calcoliamo il numero di triangoli candidati
  deg <- igraph::degree(net3)
  open_triangles <- sum(deg * (deg-1) / 2)
  if(open_triangles > 0)
    # la frazione di chiusure triadiche rispetto alle candidate
    closure_prob <- sum(count_triangles(net3)) / open_triangles
  else 
    closure_prob <- 0
  # mi salvo i valori per ogni mese
  opt <- c(opt, open_triangles)
  closure_prob_ot <- c(closure_prob_ot, closure_prob)
  graph_list[[m]] <- as.undirected(net3)
  # visualizziamo l'evoluzione della rete mese per mese
  plot(
    net3, 
    main=paste("Rete al mese:", m),
    layout=layout2,
    vertex.label=NA,
    edge.arrow.size=.2,
    vertex.size=deg/10+2
  )
}


# 15. risultati ottenuti dal processo al punto 14
triangles_over_time
opt
closure_prob_ot

# grafico che raffigura l'evoluzione della probabilità di chiusura mese per mese
closure_df <- data.frame(months, closure_prob_ot)
ggplot(closure_df, aes(x = months, y = closure_prob_ot, group = 1)) +
  geom_line() +
  geom_point() +
  labs(title = "Evoluzione Probabilità di Chiusura per Mese",
       x = "Mese",
       y = "Probabilità Chiusura") +
  theme_minimal()

# 16. analisi della formazione dei legami
# matrici con 2 colonne: 0 legami in comune tra i due utenti e 1 o più legami in comune tra i due utenti
T_k <- matrix(0, nrow = length(months)-1, ncol = 2) # numero di legami creati nel mese nei 2 casi
pairs <- matrix(0, nrow = length(months)-1, ncol=2) # numero di coppie di utenti che possono creare una connessione durante il mese nei 2 casi
T_k_frac <- matrix(0, nrow = length(months)-1, ncol = 2) # rapporto tra le 2 misure

for (m in 1:(length(months)-1)){
  # snapshot della rete sociale a inizio e fine mese con rispettive matrici di adiacenze (usate per calcolare le connessioni in comune tra utenti)
  snap1 <- graph_list[[months[m]]]
  snap2 <- graph_list[[months[m+1]]]
  adj1 <- as_adj(snap1, sparse = FALSE)
  adj2 <- as_adj(snap2, sparse=FALSE)
  # connessioni in comune tra coppie di utenti 
  com_connection <- adj1 %*% adj1 # matrice contenente in ogni cella il numero di path di lunghezza 2 tra utente i e j e quindi il numero di vicini in comune tra i 2 utenti, dove il vicino rappresenta il nodo intermedio
  # vogliamo ora ottenere le coppie che non hanno un arco in adj1 e che abbiano formato un arco durante il mese e quindi in adj2 è presente
  n <- nrow(adj1)
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      if(adj1[i,j]==0){  # se non sono connesse nel primo snapshot
        k <- com_connection[i,j] # il numero di vicini in comune
        if(k>1) # se maggiore o uguale a 1 operiamo nella seconda colonna delle matrici
          k <- 1
        
        pairs[m,k+1] <- pairs[m,k+1] + 1 # aumentiamo il numero di coppie di utenti che possono creare una connessione durante il mese nel caso selezionato
        if(adj2[i,j]==1)
          T_k[m, k+1] <- T_k[m, k+1] + 1 # aumentiamo il numero di legami creati durante il mese nel caso selezionato
      }
    }
  }
  T_k_frac[m,] <- T_k[m,] / pairs[m,] # rapporto tra le 2 misure nel mese
}

# visualizziamo i risultati ottenuti:

# 17. visualizziamo il numero assoluto di connessioni data la presenza di vicini in comune
Tk_df <- as.data.frame(T_k, months[-1]) # dataframe per il ggplot
Tk_df$month <- factor(rownames(Tk_df), levels = rownames(Tk_df))
Tk_long <- pivot_longer(Tk_df, cols = c(V1, V2), names_to = "Categoria", values_to = "Valore") # il valore deve essere di tipo long per il ggplot
Tk_long$category <- factor(Tk_long$Categoria, 
                            levels = c("V1", "V2"),
                            labels = c("0 common neighbor", "≥1 common neighbor"))
ggplot(Tk_long, aes(x = month, y = Valore, fill = category)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("azure3", "aquamarine3")) +
  labs(x = "Mese", y = "Connessioni Formate", fill = "Categoria") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("N° Connessioni data presenza Vicini in Comune")

# 18. calcolo il risultato complessivo, sommo il numero di connessioni formate e connessioni possibili di ogni mese e faccio il rapporto
T_k_tot <- colSums(T_k)
pairs_tot <- colSums(pairs)
T_k_tot <- T_k_tot / pairs_tot
# visualizziamo il risultato complessivo, cioè la probabilità di formazione connessione data la presenza di vicini in comune
T_k_tot_df <- data.frame(
  Categoria = c("0 common neighbors", "≥1 common neighbors"),
  Probabilità = T_k_tot
)

ggplot(T_k_tot_df, aes(x = Categoria, y = Probabilità, fill = Categoria)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = c("azure3", "aquamarine3")) +
  labs(title = "Probabilità Formazione Connessione data presenza Vicini in Comune",
       y = "Probabilità Formazione Connessione",
       x = "") +
  theme_minimal() +
  theme(legend.position = "none")

# ******************************************************************************

# PAROLE CHIAVE ED EMOZIONI

# 19. analisi bigrammi:
# creo un dataframe composto da utente e tweet per poi togliere dal messaggio tutto il rumore (menzioni, hashtag, rt, link)
tweets <- df[, c("username", "tweets")]
tweets$tweets <- gsub("ENGLISH TRANSLATION:|rt", "", tweets$tweets, ignore.case = TRUE)
tweets$tweets <- gsub("http[^[:space:]]+", "", tweets$tweets)
tweets$tweets <- gsub("@\\w+", "", tweets$tweets)
tweets$tweets <- gsub("#\\w+", "", tweets$tweets)

head(tweets)

# creo il corpus 
corpus = corpus(tweets, text_field = "tweets")
doc.tokens = tokens(corpus) # tokenizzo i messaggi splitto ogni documento in token singoli
# rimuovo tutte le forme di rumore: punteggiatura, numeri, stopwords
doc.tokens = tokens(doc.tokens, remove_punct = TRUE, remove_numbers = TRUE) 
doc.tokens = tokens_select(doc.tokens, stopwords(language = "en", source = "snowball", simplify = TRUE), selection ='remove') # we remove common english words like "is", "the", "and"
doc.tokens = tokens_tolower(doc.tokens) # tutte le parole convertite in lower case
doc.tokens <- tokens_keep(doc.tokens, pattern = "^[a-z]+$", valuetype = "regex") # teniamo solo token composti da parole con alfabeto classico

toks_ngram = tokens_ngrams(doc.tokens, n = 2) # bigrammi

dfmat = dfm(toks_ngram) %>% dfm_trim(min_termfreq = 20) # creo la DFM (righe sono i document, le colonne this create a DFM (rows are the documents, columns sono parole del dizionario)

features_dfm = textstat_frequency(dfmat, n = 100)
features_dfm$feature = with(features_dfm, reorder(feature, -frequency))
# rappresentazione delle frequenze sotto forma di tabella
features_dfm 
# rappresentazione delle frequenze sotto forma di textplot
textplot_wordcloud(dfmat)

#  20. analisi degli hashtag
hashtags <- df[, c("username", "tweets")]
# interessati solo negli hashtag
hashtags$hashtag_list <- regmatches(hashtags$tweets, gregexpr("#\\w+", hashtags$tweets))
hashtags$hashtag <- sapply(hashtags$hashtag_list, paste, collapse=" ")
# stessa procedura dei bigrammi:
corpus = corpus(hashtags, text_field = "hashtag")
summary(corpus)
toks_hashtag = tokens(corpus) 
toks_hashtag = tokens(toks_hashtag, remove_punct = TRUE) 
dfmat2 = dfm(toks_hashtag) %>% dfm_trim(min_termfreq = 20)
features_dfm = textstat_frequency(dfmat2, n = 100)
features_dfm$feature = with(features_dfm, reorder(feature, -frequency))
features_dfm
textplot_wordcloud(dfmat2)


# 21. analisi dei sentimenti

# uso il lexicon NRC Emotion (lista di parole e il loro valore associato per i vari sentimenti)
sentiment = get_nrc_sentiment(tweets$tweets)
td = data.frame(t(sentiment))
td = data.frame(rowSums(td[-1]))
names(td)[1] <- "count"
tdw <- cbind("sentiment" = rownames(td), td)
rownames(tdw) <- NULL
tdw


# 22. visualizzazione delle emozioni
ggplot(tdw[1:8, ], aes(x = sentiment, y = count, fill = sentiment)) +
  geom_bar(stat = "identity") +
  labs(x = "Emozioni", y="Conteggio") +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.title = element_blank())+
  ggtitle("Distribuzione Emozioni")

# 23. visualizzazione delle polarità
ggplot(tdw[9:10, ], aes(x = sentiment, y = count, fill = sentiment)) +
  geom_bar(stat = "identity") +
  labs(x = "Polarità", y="Conteggio") +
  theme(axis.text.x=element_text(angle=45, hjust=1), legend.title = element_blank()) +
  ggtitle("Distribuzione Polarità")



# ******************************************************************************

# ANALISI DELLE COMUNITÀ

# 24. inizio creando le comunità utilizzando il metodo di Louvain
communities <- multilevel.community(net2, weights = E(net2)$weight) # creo le comunità in base alla forza dei loro legami, in modo da considerare quanto sia stretto un legame per la formazione delle comunità
V(net2)$community <- membership(communities) # assegno ad ogni nodo la comunità di appartenenza
V(net2)$community 
inv_weights <- 1 / E(net2)$weight

# analisi strutturale:
# dataframe che descrive la struttura delle diverse comunità nel grafo completo, inizialmente è un dataframe contenente gli utenti ma poi tramite funzione di aggregazione prendo la media dei valori di ogni comunità
centrality_df <- data.frame(
  username = V(net2)$name,
  community = V(net2)$community,
  degree = igraph::degree(net2),
  strength = strength(net2, weights = E(net2)$weight),
  centrality = eigen_centrality(net2, weights = E(net2)$weight)$vector,  # centralità del nodo
  betweenness = betweenness(net2, weights = inv_weights), # quanto un utente controlla il flusso di informazioni
  closeness = closeness(net2, weights = inv_weights)    # quanto velocemente un utente può raggiungere gli altri utenti
)
centrality_df
aggregate(. ~ community, data = centrality_df[-1], FUN = mean) # tramite funzione di aggregazione prendo la media dei valori per ogni comunità

# 25. visualizziamo la rete sociale dove il colore del nodo rappresenta la comunità di cui fa parte
colorVals <- c("firebrick3", "slateblue1", "yellow1", "olivedrab1", "pink", "seagreen1", "orange", "turquoise1") # i colori sono in ordine di comunità quindi la comunità 1 sarà rossa, la 2 sarà violetta, ...
V(net2)$color <- colorVals[V(net2)$community]
plot.igraph(net2, vertex.label = NA,layout=layout, vertex.size=5)



# 26. struttura interna delle comunità:
community_id <- unique(V(net2)$community)
edge_density(net2, loops=F)
density <- c()
nnodes <- c()
diameter <- c()
transitivity <- c()
# calcolo le misure per ogni comunità e le inserisco in un dataframe
for (i in community_id){
  sub_net <- induced.subgraph(net2, V(net2)[community==i])
  density <- c(density, edge_density(sub_net, loops=F))
  nnodes <- c(nnodes, vcount(sub_net))
  diameter <- c(diameter, length(get_diameter(sub_net)))
  transitivity <- c(transitivity, transitivity(sub_net))
}
community_stats <- data.frame(
  community = community_id,
  density = density,
  nodes = nnodes,
  diameter = diameter,
  transitivity = transitivity
)
community_stats

# assortatività: valore che indica se un utente tende a connettersi a utenti della stessa comunità
assortativity_nominal(net2, types = as.factor(V(net2)$community), directed=F)


# 27. analisi dei contenuti:
# nuovo dataframe simile a quello per l'analisi delle di parole chiave ed emozioni solo che questa volta aggiungiamo anche la variabile comunità
comtweets <- tweets %>% inner_join(centrality_df %>% select(username, community)) 
head(comtweets)
# guardiamo i contenuti delle comunità:
community_ngrams <- list()
dfmat_coms <- list()
community_id <- unique(comtweets$community)
# questa è la stessa procedura usata per l'analisi dei bigrammi:
for (i in community_id){
  corpus = corpus(comtweets[comtweets$community==i,], text_field = "tweets")
  summary(corpus)
  doc.tokens = tokens(corpus) 
  doc.tokens = tokens(doc.tokens, remove_punct = TRUE, remove_numbers = TRUE)
  doc.tokens = tokens_select(doc.tokens, stopwords(language = "en", source = "snowball", simplify = TRUE), selection ='remove')
  doc.tokens = tokens_tolower(doc.tokens) 
  doc.tokens <- tokens_keep(doc.tokens, pattern = "^[a-z]+$", valuetype = "regex")
  community_ngrams[[i]] <- tokens_ngrams(doc.tokens, n = 2)
  dfmat_coms[[i]] <- dfm(community_ngrams[[i]]) %>% dfm_trim(min_termfreq = 8)
}

# 28. vediamo come si comportano le comunità più grandi (le prime 4):
freq_text <- list()
par(mfrow=c(2,1))
for (i in (1:4)){
  features_dfm = textstat_frequency(dfmat_coms[[i]], n = 20)
  features_dfm$feature = with(features_dfm, reorder(feature, -frequency))
  freq_text[[i]] <- features_dfm$feature
  #textplot_wordcloud(dfmat_coms[[i]])
}
freq_text

# 29. passiamo ora all'analisi dell'emozioni delle diverse comunità
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


# 30. visualizziamo le emozioni emanate dai messaggi delle 4 comunità più grandi
plots <- list()
for (i in (1:4)){
  p <- ggplot(tdw_coms[[i]][1:8, ], aes(x = sentiment, y = count, fill = sentiment)) +
    geom_bar(stat = "identity") +
    labs(x = "emotion", title = paste("Community", i)) +
    theme(axis.text.x=element_text(angle=45, hjust=1), legend.title = element_blank())
  #print(tdw_coms[[i]])
  plots[[i]] <- p
}
# dato che con ggplot, par(mfrow()) non funziona, utilizziamo un altra libreria gridExtra
grid.arrange(grobs = plots, ncol = 2)


# 31. analsi degli utenti nelle comunità
# otteniamo l'utente più influente di ciascuna delle 4 comunità più grandi (tramite il valore di eigencentrality)
top_users <- centrality_df %>%
  filter(community %in% c(1, 2, 3, 4)) %>%
  group_by(community) %>%
  top_n(1, centrality) %>% 
  arrange(community, desc(centrality))
top_users # l'utente più influente di ciascuna comunità con le corrispettive caratteristiche
users_id <- top_users$username

# compariamo i contenuti degli utenti e delle corrispettive comunità (senza tenere conto degli utenti stessi):
# 32. bigrammi più utilizzati dai 4 utenti:
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
  dfmat_users[[i]] <- dfm(users_ngrams[[i]]) %>% dfm_trim(min_termfreq = 3)
}

# visualizziamo i bigrammi più utilizzati:
freq_text <- list()
par(mfrow=c(2,1))
for (i in users_id){
  features_dfm = textstat_frequency(dfmat_users[[i]], n = 20)
  features_dfm$feature = with(features_dfm, reorder(feature, -frequency))
  freq_text[[i]] <- features_dfm$feature
}
freq_text

# 33. bigrammi più utilizzati dalle 4 comunità senza tenere conto dell' utente più influente di ciascuna comunità
comtweets <- comtweets %>% filter(!username %in% c("mobi_ayubi", "WarReporter1", "Nidalgazaui", "MaghrebiQM	")) # rimuovo gli utenti più influenti dal dataframe
for (i in community_id){
  corpus = corpus(comtweets[comtweets$community==i,], text_field = "tweets")
  summary(corpus)
  doc.tokens = tokens(corpus) # tokenize the text (split each document into individual tokens)
  doc.tokens = tokens(doc.tokens, remove_punct = TRUE, remove_numbers = TRUE)
  doc.tokens = tokens_select(doc.tokens, stopwords(language = "en", source = "snowball", simplify = TRUE), selection ='remove')
  doc.tokens = tokens_tolower(doc.tokens)
  doc.tokens <- tokens_keep(doc.tokens, pattern = "^[a-z]+$", valuetype = "regex")
  
  community_ngrams[[i]] <- tokens_ngrams(doc.tokens, n = 2)
  dfmat_coms[[i]] <- dfm(community_ngrams[[i]]) %>% dfm_trim(min_termfreq = 8)
}

# visualizziamo i bigrammi più utilizzati
freq_text <- list()
par(mfrow=c(2,1))
for (i in (1:4)){
  features_dfm = textstat_frequency(dfmat_coms[[i]], n = 20)
  features_dfm$feature = with(features_dfm, reorder(feature, -frequency))
  freq_text[[i]] <- features_dfm$feature
  #textplot_wordcloud(dfmat_coms[[i]])
}
freq_text

