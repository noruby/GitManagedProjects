data <- read.csv("./university_ranking.csv", head=T, row.name=1)
pca <-  prcomp(data, scale=T)
png("./university_ranking_PCA.png", width = 1500, height = 1500)  
biplot(pca) 
dev.off() 

