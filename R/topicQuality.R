topicQuality <- function(model, documents, xlab="Semantic Coherence", ylab="Exclusivity", labels=1:ncol(model$theta), M=10,...){
  if(length(model$beta$logbeta)<2){
    semcoh <- semanticCoherence(model,documents=documents, M=M)
    exclusivity <- exclusivity(model, M=M)
    print(semcoh)
    print(exclusivity)
    plot(0, 0, xlab=xlab, ylab=ylab, col="white", xlim=c(min(semcoh), max(semcoh)),
         ylim=c(min(exclusivity), max(exclusivity)),...)
    for(i in 1:length(labels)){
      text(semcoh[i], exclusivity[i], paste("Topic", labels[i]))
    }
  }
  if(length(model$beta$logbeta)>1){
    semcoh <- semanticCoherence(model,documents=documents, M=M)
    for(i in 1:length(labels)){
      print(paste("Topic", i, "has", semcoh[i], "semantic coherence"))
    }
  }
}
