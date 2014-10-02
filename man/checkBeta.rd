\name{checkBeta}
\alias{checkBeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Looks for words that load exclusively onto a topic
}
\description{
Checks the log beta matrix for values too close to 0, which reflect words that load onto a single topic. 
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
checkBeta(stmobject, tolerance=0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stmobject}{
     STM Model Output
      }
  \item{tolerance}{
	  User specified input reflecting closeness to 1.  E.g. a tolerance of .01 will flag 
    any values greater than .99.  Tolerance must be above 1e-6. 
      }}

\value{
  \item{problemTopics}{A list of vectors, each vector corresponding to the set of topics
  in the relevant beta matrix that contain words with too high of a loading to that topic
  }
  \item{topicErrorTotal}{A list of integers, each corresponding to the total number of 
  topics with problems in the relevant beta matrix}
  \item{problemWords}{A list of matrices, each corresponding to a relevant beta matrix, 
  which gives the topic and word index of each word with too high of a topic loading}
  \item{wordErrorTotal}{A list of integers, each corresponding to the total words 
  with problems for the relevant beta matrix}
  \item{check}{A boolean representing if the check was passed. If wordErrorTotal is all
  0s (no errors), check is True.}
}
    
\details{
  The function checks the log beta matrix for values that exceed the tolerance threshold, indicating 
  that a word has loaded onto a single topics. The ouput gives the user lists of which topics 
  have problems, which words in which topics have problems, as well as a count of the total problems
  in topics and the total number of problem words.

  Note that if the tolerance value is below 1e-6, this function will throw an error.
}
\author{
Antonio Coppola
}


