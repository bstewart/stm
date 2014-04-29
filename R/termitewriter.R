#termite writer function
#app.path is the folder to write the json files to
#model is either the model output itself 
#meta is the metadata itself 
#DocID is the optional name of docID within meta

termitewriter <- function(app.path, model, meta, DocID=NULL){

  if(!require(jsonlite)) stop("Install the jsonlite package to use this function")
  
  if(!is.null(DocID)){
    docIDs = meta[,c(DocID)]
  }else{
    docIDs = 1:nrow(meta)
  }


  data.DocTopicMatrix = model$theta
  data.TermTopicMatrix = exp( t( model$beta$logbeta[[1]] ) )

    
# Define output filenames
  app.path.TermTopicMatrix = paste(app.path, "/", "term-topic-matrix.txt", sep="")
  app.path.DocTopicMatrix = paste(app.path, "/", "doc-topic-matrix.txt", sep="")
  app.path.TermIndex = paste(app.path, "/","term-index.json", sep="")
  app.path.DocIndex = paste(app.path, "/", "doc-index.json", sep="")
  app.path.TopicIndex = paste(app.path, "/", "topic-index.json", sep="")

# Document Index
  temp.DocCount <- nrow(model$theta)
  temp.DocIDs <- paste( "Document #", 1:temp.DocCount, sep = "" )
  temp.DocIndex <- 1:temp.DocCount
  temp.DocIndexValues <- cbind( temp.DocIndex, docIDs )
  temp.DocIndexHeader <- c( "index", "docID" )
  colnames( temp.DocIndexValues ) <- temp.DocIndexHeader
  data.DocIndexJSON <- toJSON( as.data.frame( temp.DocIndexValues ), pretty = TRUE, digits = 10 )
  write( data.DocIndexJSON, file = app.path.DocIndex )

# Term Index
  temp.TermCount <- nrow( data.TermTopicMatrix )
  temp.TermFreq <- apply( data.TermTopicMatrix, 1, sum )
  temp.TermText <- model$vocab
  temp.TermIndex <- 1:temp.TermCount
  temp.TermIndexValues = cbind( temp.TermIndex, temp.TermFreq, temp.TermText )
  temp.TermIndexHeader = c( "index", "freq", "text" )
  colnames( temp.TermIndexValues ) <- temp.TermIndexHeader
  data.TermIndexJSON <- toJSON( as.data.frame( temp.TermIndexValues ), pretty = TRUE, digits = 10 )
  write( data.TermIndexJSON, file = app.path.TermIndex )

# Topic Index
  temp.TopicCount <- ncol( data.TermTopicMatrix )
  temp.TopicFreq <- apply( data.TermTopicMatrix, 2, sum )
  temp.TopicIndex <- 1:temp.TopicCount
  temp.TopicIndexValues = cbind( temp.TopicIndex, temp.TopicFreq )
  temp.TopicIndexHeader = c( "index", "freq" )
  colnames( temp.TopicIndexValues ) <- temp.TopicIndexHeader
  data.TopicIndexJSON <- toJSON( as.data.frame( temp.TopicIndexValues ), pretty = TRUE, digits = 10 )
  write( data.TopicIndexJSON, file = app.path.TopicIndex )
  
# Doc-Topic Matrix
# Tab-separated with no headers. Theta (D by K)
  rownames( data.DocTopicMatrix ) <- temp.DocIDs
  colnames( data.DocTopicMatrix ) <- temp.TopicIndex
  data.DocTopicMatrixJSON <- toJSON( data.DocTopicMatrix, pretty = TRUE, digits = 10 )
  write( data.DocTopicMatrixJSON, file = app.path.DocTopicMatrix )


# Term-Topic Matrix
# Tab-separated with no headers. Beta (V by K)
  rownames( data.TermTopicMatrix ) <- temp.TermText
  colnames( data.TermTopicMatrix ) <- temp.TopicIndex
  data.TermTopicMatrixJSON <- toJSON( data.TermTopicMatrix, pretty = TRUE, digits = 10 )
  write( data.TermTopicMatrixJSON, file = app.path.TermTopicMatrix )
  
}
