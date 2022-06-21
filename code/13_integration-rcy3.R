library(RCy3)
library(igraph)
cytoscapePing()

load("Data/databases.RData") # mirna-tf, tf-gene and gene-gene interactions ####
load("Data/d.fgl.Rdata") # first level mirna-mirna (GLASSO) ####
load("Data/targets.Rdata") # second level mirna-gene (mirnatap) ####

# 1 miRNA - miRNA ####
igraph::as_data_frame(d.fgl$CCR$Blood$d.fgl$Delta.graph.connected) |> 
  dplyr::mutate(source = from, target = to, interaction = "interacts") |> 
  dplyr::select(source, target, interaction) |> as_tibble() -> d.ccr.edge

d.ccr.node <- tibble(id = vertex.attributes(d.fgl$CCR$Blood$d.fgl$Delta.graph.connected)$name, 
                     group = "miRNA")

# 2 miRNA - gene ####
targets.unique |> 
  as_tibble() |> dplyr::filter(mirna %in% d.ccr.edge$source) |> unique() -> d.ccr.targets

d.ccr.targets |> dplyr::mutate(source = mirna, 
                               target = symbol, 
                               interaction = "inhibits") |> 
  dplyr::select(source, target, interaction) |> unique() |> 
  as_tibble() -> d.ccr.targets.edge # edge list for mirna-gene

tibble(id = c(d.ccr.targets$mirna, d.ccr.targets$symbol),
       group = c(rep("miRNA", nrow(d.ccr.targets)),
                 rep("gene", nrow(d.ccr.targets)))) |> 
  unique() -> d.ccr.targets.node # node list for mirna-gene

# 3 miRNA - TF ####
transmir.edge |> filter(target %in% d.ccr.edge$source) |> unique() -> d.ccr.transmir.edge

tibble(id = c(d.ccr.transmir.edge$source, d.ccr.transmir.edge$target), 
       group = c(rep("TF", nrow(d.ccr.transmir.edge)), 
                 rep("miRNA", nrow(d.ccr.transmir.edge)))) |> unique() -> d.ccr.transmir.node

# 4 Gene - gene ####
mppi.edge |> filter(source %in% d.ccr.targets.edge$target) |> unique() -> d.ccr.mppi.edge

tibble(id = c(d.ccr.mppi.edge$source, d.ccr.mppi.edge$target), 
       group = "gene") |> unique() -> d.ccr.mppi.node

# 5 TF - gene ####
trrust.edge |> filter(source %in% d.ccr.transmir.edge$source) |> unique() -> d.ccr.trrust.edge

tibble(id = c(d.ccr.trrust.edge$source, d.ccr.trrust.edge$target), 
       group = c(rep("TF", nrow(d.ccr.trrust.edge)), 
                 rep("gene", nrow(d.ccr.trrust.edge)))) |> unique() -> d.ccr.trrust.node

edgelist <- rbind(d.ccr.edge, d.ccr.targets.edge, d.ccr.transmir.edge, d.ccr.mppi.edge, d.ccr.trrust.edge)
nodelist <- rbind(d.ccr.node, d.ccr.targets.node, d.ccr.transmir.node, d.ccr.mppi.node, d.ccr.trrust.node)


createNetworkFromDataFrames(nodes = nodelist, edges = edgelist, title = "d.CCR", collection = "Blood 1")
createVisualStyle(style.name = "Regulatory", 
                  mappings = list(nodeLabels = mapVisualProperty('node label', 
                                                                 'id', 
                                                                 'p'),
                                  nodeFills = mapVisualProperty('node fill color',
                                                                'group',
                                                                'd',
                                                                c("miRNA", "gene", "TF"),
                                                                c("#a37a3e", "#5cc166","#343a69")), 
                                  arrowShapes = mapVisualProperty('Edge Target Arrow Shape', 
                                                                  'interaction', 
                                                                  'd',
                                                                  c("activates", "inhibits", "interacts"),
                                                                  c("Arrow", "T", "None"))))
setVisualStyle(style.name = "Regulatory")

g <- igraph::graph_from_data_frame(edgelist, directed = FALSE, vertices = nodelist)


table <- read_delim("Data/d.CCR default node.csv", delim = ",")







nodes <- data.frame(
  id = c("node 0", "node 1", "node 2", "node 3"),
  group = c("A", "A", "B", "B"),
  # categorical strings
  score = as.integer(c(20, 10, 15, 5)),
  # integers
  stringsAsFactors = FALSE
)
edges <- data.frame(
  source = c("node 0", "node 0", "node 0", "node 2"),
  target = c("node 1", "node 2", "node 3", "node 3"),
  interaction = c("inhibits", "interacts", "activates", "interacts"),
  # optional
  weight = c(5.1, 3.0, 5.2, 9.9),
  # numeric
  stringsAsFactors = FALSE
)

createNetworkFromDataFrames(nodes, edges, title = "my first network", collection =
                              "DataFrame Example")

style.name = "myStyle"
defaults <- list(
  NODE_SHAPE = "diamond",
  NODE_SIZE = 30,
  EDGE_TRANSPARENCY = 120,
  NODE_LABEL_POSITION = "W,E,c,0.00,0.00"
)
nodeLabels <- mapVisualProperty('node label', 'id', 'p')
nodeFills <-
  mapVisualProperty('node fill color',
                    'group',
                    'd',
                    c("A", "B"),
                    c("#FF9900", "#66AAAA"))
arrowShapes <-
  mapVisualProperty('Edge Target Arrow Shape', 'interaction', 'd',
    c("activates", "inhibits", "interacts"),
    c("Arrow", "T", "None"))
edgeWidth <- mapVisualProperty('edge width', 'weight', 'p')

createVisualStyle(style.name,
                  defaults,
                  list(nodeLabels, nodeFills, arrowShapes, edgeWidth))
setVisualStyle(style.name)