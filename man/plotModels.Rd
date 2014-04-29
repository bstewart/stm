\name{plotModels}
\alias{plotModels}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots semantic coherence and exclusivity for high likelihood models
outputted from selectModel.}
\description{
  Plots semantic coherence and exclusivity for high likelihood models. 
  In the case of models that include content covariates, prints semantic
  coherence and sparsity.
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotModels(models, xlab="Semantic Coherence",
    ylab="Exclusivity", labels=1:length(models$runout),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{models}{
    output from selectModel.
  }
  \item{labels}{
    labels for each model.
    }
             \item{xlab}{
Character string that is x axis title. This will be semantic coherence.
}
         \item{ylab}{
Character string that is y axis title. This will be exclusivity.
}
  \item{...}{
    Other plotting parameters.
  }
}

\details{
  Each model has semantic coherence and exclusivity values associated
  with each topic.  In the default plot function, the small colored dots
  are associated with a topic's semantic coherence and exclusivity.
  Dots with the same color as topics associated with the same model.
  The average semantic coherence and exclusivity is also plotted in the
  same color, but printed as the model number associated with the output
  from selectModels().

  With content covariates, the model does not output exclusivity because
  exclusivity has been built in with the content covariates.  Instead,
  the user should check to make sure that sparsity is high enough
  (typically greater than .5), and then should select a model based on
  semantic coherence.
%%  ~~ If necessary, more details than the description above ~~
}

\examples{


\dontrun{
#storage is an object created by selectModel
plotModels(storage)
#this selects the first model from selectModel
selected<-storage$runout[[1]]
}
 
}
