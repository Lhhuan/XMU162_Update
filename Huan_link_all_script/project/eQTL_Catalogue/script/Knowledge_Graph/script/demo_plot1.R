# library(networkdata)
library(igraph)
library(ggraph)
library(graphlayouts)

actors <- data.frame(name=c("hotspot", "gene1","gene2","gene3","gene4", "enhancer1", "H3K27ac",
                            "H3K4me1","tissue1","tissue2"),
                     Node=c("Hotspot","Gene","Gene","Gene","Gene","CRE","Epigenetic marker","Epigenetic marker","Tissue","Tissue"))

relations <- data.frame(from=c("hotspot", "gene1", "gene1", "gene1","hotspot", "hotspot","hotspot","enhancer1","gene1"),
                        to=c("gene1", "gene2", "gene3", "gene4", "enhancer1", "H3K27ac","H3K4me1","tissue2","tissue1"),
                        Edge=c("Hotspot-Target","Co-expression","Functional interaction","Functional interaction",
                                "Hotspot annotation","Hotspot annotation","Hotspot annotation","Tissue","Tissue"))

g <- graph_from_data_frame(relations, directed=T, vertices=actors)

nodes <-c(relations$from,relations$to)
# require(extrafont)

p <- ggraph(g, "stress", bbox = 15) +
  geom_edge_link2(aes(edge_colour = Edge), edge_width = 0.5) +
  geom_node_point(aes(fill = Node), shape = 21, size = degree(g)) +
  geom_node_text(aes(label = name, size = 3), repel = TRUE) +
  scale_edge_colour_brewer(palette = "Set1") +
  # scale_fill_manual(values = c("grey66", "#EEB422", "#424242")) +
  scale_size(range = c(2, 5), guide = "none") +
  theme_graph(base_family = "serif") +
  theme(legend.position = "right")

pdf("demo.pdf")
print(p)
dev.off()