\name{NEWS}
\title{NEWS file for the stm package}
\section{Version 1.1.3}{
\itemize{
\item Fixed the unexported findTopic()
\item Improved some documentation
\item Small finetuning in toLDAViz
}
}
\section{Version 1.1.2}{
\itemize{
\item Fixed a small bug that caused readCorpus to fail on dense document term matrices
\item Fixed a small bug in the random projections algorithm
\item Improved warnings in stm when restarting models (Hat tip to Andrew Goldstone)
\item Added the convertCorpus function for converting stm to other formats. 
\item Formatting changes to the vignette
}
}
\section{Version 1.1.1}{
\itemize{
\item Updated Vignette
}}
\section{Version 1.1.0}{
\itemize{
\item Performance improvements via various optimizations including porting some components to C++
\item Various new experimental features including K=0
\item Improved documentation including a new version of the vignette.
\item Better error messages in several places
\item Experimental options for random projections with spectral initializations
}}
\section{Version 1.0.13}{
\itemize{
\item Fixes a problem in make.heldout where a document could be completely emptied by the procedure. Hat tip to Jesse  Rhodes for the bug report.
}}
\section{Version 1.0.12}{
\itemize{
\item When \code{gamma.prior="L1"} coerce the mu object back to a matrix class object.  Should fix a speed hit introduced in 1.0.10 for this case.
}}
\section{Version 1.0.10}{
\itemize{
\item Prevalence covariates can now use sparse matrices which will result in better performance for large factors.
\item textProcessor() and prepDocuments() now do a better job of preserving labels and keeping track of dropped elements.  Special thanks to Github users gtuckerkellog and elbamos for pull requests.
\item Fixed an edge case in init.type="Spectral" where words appearing only in documents by themselves would throw an error.  The error was correct but hard to address in certain cases, so now it temporarily removes the words and then reintroduces them before starting inference spreading a tiny bit of mass evenly across the topics.  Hat tip to Nathan Sanders for brining this to our attention.  
\item New function findTopic() which helps locate topics containing particular words or phrases.
\item New function topicLasso() helps build predictive models with topics.
}}
\section{Version 1.0.8}{
\itemize{
\item Fixed a minor bug in prepDocuments which arises in cases where there are vocab elements which do not appear in the data.
\item Fixed a minor bug in frex calculation that caused some models not to label.
\item Fixed a minor bug in searchK that caused heldout results to report incorrectly.
\item Rewrite of plot.estimateEffect() which fixed a bug in some interaction models.  Also returns results invisibly for creating custom plots.
\item Increased the stability of the spectral methods for stm initialization.
\item Complete rewrite of plotRemoved() which makes it much faster for larger datasets.
}
}
\section{Version 1.0.1}{
\itemize{
\item A minor patch to deal with textProcessor() in older versions of R.
}
}
\section{Version 1.0.0}{
  \itemize{
    \item Large changes many of which are not backwards compatible.
    \item Numerous speed improvements to the core algorithm.
    \item Introduction of several new options for the core stm function including spectral initalization, memoized inference, and model restarts.
    \item Content covariate models are now estimated using the distributed multinomial formulation which is dramatically faster. Default prior also changed to L1.
    \item Handling of document level convergence was changed to ensure positive definiteness in the document-level covariance matrices
    \item Fixed bug in binary/binary interactions.
    \item Numerous new diagnostic and summary functions
    \item Expanding the console printing of many of the preprocessing functions
    }}
\section{Version 0.6.24}{
  \itemize{
    \item Fix an error with vignettes building on linux machines
    \item sageLabels exported but not documented
    \item factorCheck diagnostic function exported
}}
\section{Version 0.6.19}{
  \itemize{
    \item Bug fix in the semantic Coherence function that affected content covariate models.
    \item Bug fix to plot.STM() where for content covariate models with  only a subset of topics requested the labels would show up as mostly NA. Thanks to Jetson Leder-Luis for pointing this out.
    \item Bug fix for the readCorpus() function with txtorg vocab.  Thanks to Justin Farrell for pointing this out.
    \item Added some diagnostics to notify the user when words have been dropped in preprocessing.
    \item Automatically coerce dates to numeric in spline function.
}}
\section{Version 0.6.16}{
  \itemize{
    \item Very minor change with textProcessor() to accomodate API change in tm version 0.6
    \item New option for plot.STM() which plots the distribution of theta values.  Thanks
          to Antonio Coppola for coauthoring this component.
    \item Deprecated option "custom" in "labeltype" of plot.STM().  Now you can simply specify the
          labels.  Added additional functionality to specify custom topic names rather than the
          default "Topic #:"
    \item Bug fixes to various portions of plot.STM() that would cause labels to not print.
    \item Added numerous error messages.
}}
\section{Version 0.6.12}{
  \itemize{
    \item Added permutationTest() function and associated plot capabilities
    \item Updates to the vignette.
    \item Added functionality to a few plotting functions.
    \item When using summary() and labelTopics() content covariate models now have labels
          thresholded by a small value.  Thus one may see no labels or very few labels
          particularly for topic-covariate interactions which indicates that there are no
          sizable positive deviations from the baseline.
}}
\section{Version 0.6.9}{
  \itemize{
    \item S3 method for findThoughts and ability to threshold by theta.
    \item Allow estimateEffect() to receive a data frame. (Thanks to Baoqiang Cao for pointing this out)
}}
\section{Version 0.6.7}{
  \itemize{
    \item Major updates to the vignette
    \item Minor Updates to several plotting functions
    \item Fixed an error where labelTopics() would mislabel when passed topic numbers out of order
          (Thanks to Jetson Leder-Luis for pointing this out)
    \item Introduction of the termitewriter function.
}}
\section{Version 0.6.1}{
  \itemize{
    \item Version for submission to CRAN (2/28/2014)
    \item Introduced new dataset poliblog5k and shrunk the footprint of the package
    \item Numerous alternate options changed and some slight syntax changes to stm to
          finalize the API.
}}
\section{Version 0.05.04}{
  \itemize{
    \item New build 2/14/2014
    \item Fixing a small bug introduced in the last version which kept defaults
          of manyTopics() from working.
}}
\section{Version 0.05.03}{
  \itemize{
    \item Updated version posted to Github (2/13/2014)

    \item Various improvements to plotting functions.
    
    \item Setting the seed in selectModel() threw an error.  
          This is now corrected. Thanks to Mark Bell for pointing this out.
}}
\section{Initial public version 0.05.01}{
  \itemize{
    \item First public version released on Github (2/5/2014)
    
    \item This is a beta release and we may change some of the API before
    submission to CRAN. 
}}

