source("model.R")
source("optimization.R")

# (D,K,sigma_t,sigma_0,nb_epochs,nb_MB,nb_NS)
model <- init_model(2,10,1,1,1,1,10)

data <- list()
for(t in 1:3){
  file <- paste(paste("input/dataToy",t,sep=""),".txt",sep = "")
  data[[t]] <- as.matrix(read.csv(file))
}

debug(train)
model <- train(data,model)
