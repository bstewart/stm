\name{NEWS}
\title{NEWS file for the stm package}
\section{Version 1.2.1}{
\itemize{
\item Small updates to the documentation and vignette.
}
}
\section{Version 1.2.0}{
\itemize{
\item Thanks to efforts by Ken Benoit stm can now take a quanteda dfm object.
\item Thanks to help from Chris Baker we are now using roxygen for our documentation.
\item Jeffrey Arnold fixed a small bug in toLDAVis.  Thanks!
\item Carsten Schwemmer helped us find a bug where plot.estimateEffect() didn't work when dplyr was loaded.  This is now fixed.
\item In making the change to roxygen we unfortunately break backwards compatability.  The package's generic functions such as plot.estimateEffect() and plot.STM() can now only be called by plot() rather than by their full name.
\item We document and export optimizeDocument which provides access to the document level E-step.
\item We have documented and exported several of the labeling functions including calcfrex, calcscore, calclift and js.estimate.  These are marked with keyword internal because they don't have much error checking and most users will want labelTopics anyhow.  But they can be accessed with ? and are linked from labelTopics
\item After much popular demand we have released a fitNewDocuments() function which will calculate topic proportions for documents not used to fit the models.  There are many different options here.
\item estimateEffect now has a summary function which will make regression tables
\item a number of the internals to plot.estimateEffect have been improved which should eliminate some edge case bugs.
\item Much of the documentation has been updated as has the vignette.
\item Our wrapper s() for the splines package function bs() now has predict functions associated with it so it should work in contexts like lm()
\item We have documented and exported all the metrics for searchK()
\item We made a change in the spectral initialization which ensures that only the top 10000 (a modifiable default) words are used in the initialization.  This allows it to be used effectively with much larger vocab.
\item Added a modifiable max iteration timeout error for the prevalence regression in stm.  This will only matter for people using covariate sets which are very, very large.
}
}
\section{Version 1.1.6}{
\itemize{
\item Added a new recovery algorithm used in spectral initialization from gradient descent to a more accurate and generally faster one based on quadratic programming.  This can be turned on by: control=list(recoverEG=FALSE).  Eventually we may change the default recovery method. Note: while this more
accurately solves the actual problem, we've seen better results with the early stopping produced by
exponential gradient not fully convering.  This has been confirmed by the Arora group as well.
\item Clarified some of the documentation in textProcessor thanks to James Gibbon.
\item Fixed a problem with registration of S3 methods for textProcessor()
\item Added access to the information criterion parameter for L1 mode prevalence covariance in stm.  See the gamma.ic.k option in the control parameters of stm 
}
}
\section{Version 1.1.5}{
\itemize{
\item searchK() can now be used with content covariates thanks to github user rosemm
\item added a querying function based on data.table into findThoughts()
}
}
\section{Version 1.1.4}{
\itemize{
\item Fixed a rare bug in the K=0 feature for spectral initialization where words with the exact same appearance pattern would cause the projection to fail.
}
}
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

