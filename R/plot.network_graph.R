
#' Function to generate an interactive topic network.  Each node is a topic, and
#' each edge shows the correlation between topics.  Figure is generated in .html
#' and javascript, and can be viewed in RStudio's Viewer or in any web browser.
#'
#' @param stmobj Model output from STM.
#'
#' @param file The complete file path for where to save the output to disk.
#'     The filename extension should be `.html`
#'
#' @param method (available options: `community` (default) and  `difference`).
#'     `community`: Node color is defined by community detection
#'     `difference`: Node color is defined by the relative topic prevalence
#'     comparing across two levels of target covariate. This is a color
#'     representation of the \code{\link[stm]{plot.estimateEffect}} forest
#'     plot.
#'
#' @param estimateEffectobj The output of \code{\link[stm]{estimateEffect}}.
#'     (only used when \code{method = 'difference'}).
#'
#' @param main The main title of the figure.
#'
#' @param custom.labels A vector of character labels for each topic.  Must be of
#'     length `K`
#'
#' @param covariate For method "difference", string of the name of the main
#'     covariate of interest. Must be enclosed in quotes. All other covariates
#'     within the formula specified in estimateEffect will be kept at their
#'     median. See \code{\link[stm]{plot.estimateEffect}} for more details.
#'
#' @param cov.value1 For method "difference", the value or set of values of
#'     interest at which to set the covariate. In the case of calculating a
#'     treatment/control contrast, set the treatment to cov.value1.See
#'     \code{\link[stm]{plot.estimateEffect}} for more details.
#'
#' @param cov.value2 For method "difference", the value or set of values of
#'     interest at which to set the covariate. In the case of calculating a
#'     treatment/control contrast, set the treatment to cov.value2.See
#'     \code{\link[stm]{plot.estimateEffect}} for more details.
#'
#' @param cov.value1.col For method "difference", the color to associate with
#'     \code{cov.value1}.
#'
#' @param cov.value2.col For method "difference", the color to associate with
#'     \code{cov.value2}.
#'
#' @param edge.length A numeric length of the edges drawn between correlated
#'     topics. This value is a constant.
#'
#' @param cutoff The minimum Pearson's correlation coefficient threshold to
#'     consider two topics correlated.  This value is passed directly to
#'     \code{\link[stm]{topicCorr}}.  Note that at small numbers of topics, the
#'     correlation between topics is biased to be smaller than 0 on average.
#'     The baseline expected correlation between any two topics is
#'     \code{1 / (K - 1)} where `K` is the number of topics.
#'
#' @param theme.labels A vector of character labels for each topic.  Must be
#'     of the same length as the number of identified themes  This
#'     argument can only be meaningfully set after examining the output of the
#'     figure with this argument not set to identify the number and labels for
#'     each theme.
#'
#' @param node.shape A character string to identify the shape of nodes.  See
#'     \code{\link[visNetwork]{visNetwork}} for more details.
#'
#' @param node.shadow A logical to identify if each node should have a shadow
#'     surrounding the border.See \code{\link[visNetwork]{visNetwork}} for more
#'     details.
#'
#' @param node.color.border The color of the node border.  See
#'     \code{\link[visNetwork]{visNetwork}} for more details.
#'
#' @param node.border.width The width of the node border.  See
#'     \code{\link[visNetwork]{visNetwork}} for more details.
#'
#' @param node.highlight.color.background The color of the node background once
#'     selected. See \code{\link[visNetwork]{visNetwork}} for more details.
#'
#' @param node.highlight.color.border The color of the node border once
#'     selected. See \code{\link[visNetwork]{visNetwork}} for more details.
#'
#' @param topic.label.metric The word weighting to use to select representative
#'     words for each topic (available options include `prob`, `frex`, `lift`,
#'     and `score`. These words appear when the mouse cursor hovers
#'     over a topic.  See \code{\link[stm]{labelTopic}} for more details.
#'
#' @param edge.weight.correlation.scaler A numeric to scale the topic
#'     correlation to edge weight.  By default, edge weights are uniform.
#'     Setting this value will scale the network edge weights.
#'
#' @return a \code{\link[visNetwork]{visNetwork}} visualization.  See
#'     \code{\link[visNetwork]{visNetwork}} for more details.
#'
#' @export
topic_network <- function(
  stmobj,
  file = NULL,
  method = "community",
  estimateEffectobj = NULL,
  main = "Topic Network Graph",
  custom.labels = NULL,
  covariate = NULL,
  cov.value1 = NULL,
  cov.value2 = NULL,
  cov.value1.col = "red",
  cov.value2.col = "blue",
  edge.length = 500,
  cutoff = 0.1,
  theme.labels = NULL,
  node.shape = "dot",
  node.shadow = TRUE,
  node.color.border = "black",
  node.border.width = 2,
  node.highlight.color.background = "orange",
  node.highlight.color.border = "darkred",
  topic.label.metric = "prob",
  edge.weight.correlation.scaler = NULL) {

  # Extract K
  K <- stmobj$settings$dim$K

  # ============================================================================
  # Create data.frame of topics and topic labels
  # ============================================================================
  if (!is.null(custom.labels)) {
    if (length(custom.labels) != K) {
      stop("Number of topic labels does not match the number of topics")
    } else {
      topic.labels <- custom.labels
    }
  } else {
    topic.labels <- c(paste0("Topic_", 1:K))
  }

  # ============================================================================
  #
  # Calculate topic correlation and create igraph nodes and edges data.frames
  #
  # ============================================================================
  topic_corr <- stm::topicCorr(stmobj, method = "simple", cutoff = cutoff)

  edges_mat <- as.matrix(topic_corr$posadj)
  graph <- igraph::graph_from_adjacency_matrix(edges_mat, mode = "undirected")
  graph <- igraph::simplify(graph, remove.multiple = F, remove.loops = T)
  edges <- igraph::as_data_frame(graph, what = "edges")

  # ============================================================================
  # Edge properties:
  #     - Edge thickness proportional to topic correlation
  #     - Edge length set by edge.length parameter
  # ============================================================================
  if (nrow(edges) == 0) {
    warning("No edges detected at the current cutoff, consider decreasing cutoff")
  } else {
    edges$length <- edge.length

    if (is.null(edge.weight.correlation.scaler)) {
      edges$width <- 5
    } else {

      edge_thickness <- NULL
      for (i in 1:nrow(edges)) {
        edge_thickness[i] <- topic_corr$poscor[edges$from[i], edges$to[i]]
      }

      edges$width <- edge_thickness * edge.weight.correlation.scaler
    }
  }

  # ============================================================================
  # Node formatting
  # ============================================================================
  nodes <- igraph::as_data_frame(graph, what = "vertices")

  # Node Identification
  nodes$id <- 1:K
  nodes$label <- topic.labels

  # Node properties
  nodes$shape <- node.shape
  nodes$shadow <- node.shadow
  nodes$borderWidth <- node.border.width

  # Node color
  nodes$color.border <- node.color.border
  nodes$color.highlight.background <- node.highlight.color.background
  nodes$color.highlight.border <- node.highlight.color.border

  # Node size
  nodes$topic.proportion <- apply(stmobj$theta, 2, FUN = mean)
  max_topic_proportion <- max(nodes$topic.proportion)
  nodes$size <- 40 * nodes$topic.proportion / max_topic_proportion

  # Node text on hover: paste representative words together
  representative_words <- stm::sageLabels(stmobj)$marginal
  nodes$title <- apply(
    representative_words[[topic.label.metric]],
    1,
    function(x) paste0(x, collapse = " + ")
  )

  # ============================================================================
  # Community Detection
  # ============================================================================
  set.seed(42)
  cluster_label_prop <- igraph::cluster_label_prop(graph)
  nodes$community <- cluster_label_prop$membership

  if (is.null(theme.labels)) {
    nodes$Themes <- nodes$community
  } else {
    nodes$Themes <- theme.labels[nodes$community]
  }

  # ============================================================================
  #
  # Differential plotting based on method
  #
  # ============================================================================

  # ============================================================================
  # "community" uses a graph community detection algorithm to color nodes
  # ============================================================================
  if (method == "community") {

    # Node colors
    color_spectrum <- RColorBrewer::brewer.pal(12, "Set3")
    nodes$color.background <- color_spectrum[nodes$community]

    # ==========================================================================
    # community labels & colors
    #     - community_color_list is used to create legend in figure
    #     - each object in list is a node in the figure legend
    #     - Community label is different if theme.labels argument was passed
    # ==========================================================================
    if (is.null(theme.labels)) {
      nodes$community_names <- nodes$community

      community_color_list <- NULL
      for (i in 1:max(nodes$community)) {
        color_label <- color_spectrum[i]

        community_color_list[[i]] <- list(
          label = i,
          size = 50,
          color = color_label
        )
      }
    } else {
      nodes$community_names <- theme.labels[nodes$community]

      community_color_list <- NULL
      for (i in 1:max(nodes$community)) {
        color_label <- color_spectrum[i]

        community_color_list[[i]] <- list(
          label = theme.labels[i],
          size = 50,
          color = color_label
        )
      }
    }

    # Create visNetwork graph
    visNetwork_graph <- visNetwork::visNetwork(
      nodes = nodes,
      edges = edges,
      width = "100%",
      height = "600px",
      main = main
    ) %>%
      visNetwork::visInteraction(navigationButtons = FALSE) %>%
      visNetwork::visIgraphLayout(randomSeed = 42) %>%
      visNetwork::visLegend(addNodes = community_color_list, useGroups = F) %>%
      visNetwork::visOptions(
        highlightNearest = list(
          enabled = TRUE,
          algorithm = "hierarchical"
        ),
        selectedBy = "Themes"
      )
  }


  # ============================================================================
  # "difference" colors the nodes using according to the relative topic
  #  prevalence commonly visualized in the plot.estimateEffect()
  # ============================================================================
  else if (method == "difference") {
    estimateEffect_plot_obj <- stm::plot.estimateEffect(
      model = stmobj,
      x = estimateEffectobj,
      method = "difference",
      covariate = covariate,
      cov.value1 = cov.value1,
      cov.value2 = cov.value2,
      omit.plot = TRUE
    )

    means <- data.frame(unlist(estimateEffect_plot_obj$means))
    colnames(means) <- "means"
    color1_ramp <- colorRamp(c("white", cov.value1.col))(abs(means$means) / 0.05)
    means$color1 <- grDevices::rgb(
      color1_ramp[, 1],
      color1_ramp[, 2],
      color1_ramp[, 3],
      maxColorValue = 255
    )

    color2_ramp <- colorRamp(c("white", cov.value2.col))(abs(means$means) / 0.05)
    means$color2 <- grDevices::rgb(
      color2_ramp[, 1],
      color2_ramp[, 2],
      color2_ramp[, 3],
      maxColorValue = 255
    )

    nodes$color.background <- ifelse(means$means > 0, means$color1, means$color2)

    # ==========================================================================
    # List of lists to create legend
    # ==========================================================================
    difference.legend.list <- list(
      list(
        label = cov.value1,
        color = cov.value1.col,
        shape = node.shape,
        color.border = node.color.border,
        shadow = node.shadow,
        borderWidth = node.border.width
      ),
      list(
        label = cov.value2,
        color = cov.value2.col,
        shape = node.shape,
        color.border = node.color.border,
        shadow = node.shadow,
        borderWidth = node.border.width
      )
    )

    # Plot the required graph
    visNetwork_graph <- visNetwork::visNetwork(
      nodes = nodes,
      edges = edges,
      width = "100%",
      height = "600px",
      main = main
    ) %>%
      visNetwork::visInteraction(navigationButtons = FALSE) %>%
      visNetwork::visIgraphLayout(randomSeed = 42) %>%
      visNetwork::visLegend(addNodes = difference.legend.list, useGroups = F) %>%
      visNetwork::visOptions(
        highlightNearest = list(
          enabled = TRUE,
          algorithm = "hierarchical"
        ),
        selectedBy = "Themes"
      )
  } else {
    print("Do not recognize method!  Check the available options")
    return(NULL)
  }

  if (!is.null(file)) {
    visNetwork_graph %>% visNetwork::visSave(file = file)
  }

  return(visNetwork_graph)
}
