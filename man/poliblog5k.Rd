\name{poliblog5k}
\alias{poliblog5k}
\alias{poliblog5k.docs}
\alias{poliblog5k.voc}
\alias{poliblog5k.meta}
\docType{data}
\title{
CMU 2008 Political Blog Corpus
}
\description{
A 5000 document sample from CMU 2008 Political Blog Corpus (Eisenstein and Xing 2010).  Blog posts from 6 blogs during the U.S. 2008 Presidential Election.
}
\usage{poliblog5k.meta}
\format{
  A data frame with 5000 observations on the following 4 variables.
  \describe{
    \item{\code{rating}}{a factor variable giving the partisan affiliation of the blog (based on who they supported for president)}
    \item{\code{day}}{the day of the year (1 to 365).  All entries are from 2008.}
    \item{\code{blog}}{a two digit character code corresponding to the name of the blog. They are: American Thinker (at), Digby (db), Hot Air (ha), Michelle Malkin (mm), Think Progress (tp), Talking Points Memo (tpm)}
    \item{\code{text}}{the first 50 characters (rounded to the nearest full word).}
  }
}
\details{
This is a random sample of the larger CMU 2008 Political Blog Corpus collected by Jacob Eisenstein and Eric Xing.  Quoting from their documentation: "[The blogs] were selected by the following criteria: the Technorati rankings of blog authority, ideological balance, coverage for the full year 2008, and ease of access to blog archives. In the general election for U.S. President in 2008, the following blogs supported Barack Obama:
Digby, ThinkProgress, and Talking Points Memo. John McCain was supported by American
Thinker, Hot Air, and Michelle Malkin. In general, the blogs that supported Obama in the election tend to advocate for similar policies and candidates as the Democratic party; and the blogs that supported McCain tend to advocate Republican policies and candidates. Digby, Hot Air and Michelle Malkin are single-author blogs; the others have multiple authors."
}
\source{
Jacob Eisenstein and Eric Xing (2010) "The CMU 2008 Political Blog Corpus." Technical Report Carnegie Mellon University. http://sailing.cs.cmu.edu/socialmedia/blog2008.html}
\examples{
data(poliblog5k)
head(poliblog5k.meta)
head(poliblog5k.voc)
\dontrun{
stm1 <- stm(poliblog5k.docs, poliblog5k.voc, 3,
prevalence=~rating, data=poliblog5k.meta)
}
}
\keyword{datasets}
