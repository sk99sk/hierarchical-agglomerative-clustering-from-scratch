install.packages("factoextra")
library(factoextra)
library(ISLR)
nci.labs=NCI60$labs
nci.data=NCI60$data
scaled.data=scale(nci.data)
data.dist = dist(scaled.data)

set.seed(123)

fviz_nbclust(scaled.data, kmeans, method = "wss")
fviz_nbclust(scaled.data, kmeans, method = "silhouette")

kresult <- kmeans(scaled.data, centers = 5, nstart = 25)

print(kresult)
