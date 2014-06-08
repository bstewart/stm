\name{NEWS}
\title{NEWS file for the stm package}
\section{Version 0.6.19}{
  \itemize{
    \item 
    \item Bug fix in the semantic Coherence function that affected content covariate models.
    \item Bug fix to plot.STM() where for content covariate models with  only a subset of topics requested the labels would show up as mostly NA. Thanks to Jetson Leder-Luis for pointing this out.
    \item Bug fix for the readCorpus() function with txtorg vocab.  Thanks to Justin Farrell for pointing this out.
    \item Added some diagnostics to notify the user when words have been dropped in preprocessing.
    \item Automatically coerce dates to numeric in spline function.
    \item Test for singular covariate matrices. 
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

