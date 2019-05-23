# Strains selection based on metadata


# 0. Libraries  
library(cluster) 
library(fpc)
library(ggplot2)
library(diplyr)
library(dendextend)

# File input name
data<-read.csv("liseq2.csv",sep=";")
 
 
 # 1. Distance
 
 gower.dist <- daisy(data, metric = c("gower"))

 # 2. Clustering
 
 clust_gow <- hclust(gower.dist, method = "complete")
 plot(clust_gow,main = "Dendrogram based on Gower")
 
 # 3. Assess clustering
 
  
# Fonction to cite the origin
 cstats.table <- function(dist, tree, k) {
   clust.assess <- c("cluster.number","n","within.cluster.ss","average.within","average.between",
                     "wb.ratio","dunn2","avg.silwidth")
   clust.size <- c("cluster.size")
   stats.names <- c()
   row.clust <- c()
   output.stats <- matrix(ncol = k, nrow = length(clust.assess))
   cluster.sizes <- matrix(ncol = k, nrow = k)
   for(i in c(1:k)){
       row.clust[i] <- paste("Cluster-", i, " size")
   }
   for(i in c(2:k)){
     print(i)
     stats.names[i] <- paste("Test", i-1)
     
     for(j in seq_along(clust.assess)){
       output.stats[j, i] <- unlist(cluster.stats(d = dist, clustering = cutree(tree, k = i))[clust.assess])[j]
       
     }
     
     for(d in 1:k) {
       cluster.sizes[d, i] <- unlist(cluster.stats(d = dist, clustering = cutree(tree, k = i))[clust.size])[d]
       dim(cluster.sizes[d, i]) <- c(length(cluster.sizes[i]), 1)
       cluster.sizes[d, i]
       
     }
   }
   output.stats.df <- data.frame(output.stats)
   cluster.sizes <- data.frame(cluster.sizes)
   cluster.sizes[is.na(cluster.sizes)] <- 0
   rows.all <- c(clust.assess, row.clust)
   # rownames(output.stats.df) <- clust.assess
   output <- rbind(output.stats.df, cluster.sizes)[ ,-1]
   colnames(output) <- stats.names[2:k]
   rownames(output) <- rows.all
   is.num <- sapply(output, is.numeric)
   output[is.num] <- lapply(output[is.num], round, 2)
   output
 }

 data_fig<-data.frame(t(cstats.table(gower.dist, aggl.clust.c, 105)))
 
 
 # Used of "Silhouette" plot to assess clustering
 
 ggplot(data = data_fig, 
        aes(x=cluster.number, y=avg.silwidth)) + 
   geom_point()+
   geom_line()+
   ggtitle("Agglomerative clustering") +
   labs(x = "Number of clusters", y = "Average silhouette width") +
   theme(plot.title = element_text(hjust = 0.5)) 
 
 
 # 4. Visualisation
 #select k
 k<-100
 dendro <- as.dendrogram(clust_gow)
 dendro.col <- dendro %>%
   set("branches_k_color", k = k, value =rainbow(n=k)    ) %>%
   set("branches_lwd", 0.6) %>%
   set("labels_colors", 
       value = c("darkslategray")) %>% 
   set("labels_cex", 0.5)
 ggclust <- as.ggdend(dendro.col)
 ggplot(ggclust, theme = theme_minimal()) +
   labs(x = "Num. observations", y = "Height", title = "Dendrogram, k")
 

 
 
 clust.num <- cutree(aggl.clust.c, k = k) 
 data.cl <- cbind(data, clust.num)
write.csv(file="liseq_100.csv",data.cl,row.names = FALSE)
 clusplot(data.cl, clust.num, 
          color=TRUE, shade=TRUE, labels=0, lines=0, 
          main = "Customer clusters (k=11)", 
          cex = 0.3)
 
 clust.num <- cutree(aggl.clust.c, k = k) 
 data.cl <- cbind(data, clust.num)
 write.csv(file="liseq_83.csv",data.cl,row.names = FALSE)
 clusplot(data.cl, clust.num, 
          color=TRUE, shade=TRUE, labels=0, lines=0, 
          main = "Customer clusters (k=83)", 
          cex = 0.3)
