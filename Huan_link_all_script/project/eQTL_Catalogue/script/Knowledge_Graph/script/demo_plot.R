# library(networkdata)
library(igraph)
library(ggraph)
library(graphlayouts)

actors <- data.frame(name=c("hotspot2","hotspot1", "gene1","gene2","gene3","gene4","gene5","gene6","gene7","gene8"),
                     Node=c("Hotspot","Hotspot","Gene","Gene","Gene","Gene","Gene","Gene","Gene","Gene"))

# rep()
relations <- data.frame(from=c("hotspot1","hotspot1", "gene1", "gene1", "gene1","gene1","hotspot2","hotspot2","gene7","gene8","gene6","gene3"),
                        to=c("gene7","gene1", "gene2", "gene3", "gene4","gene5","gene5","gene6","gene8","gene1","gene2","gene4"),
                        Edge=c("Hotspot-Target","Hotspot-Target","Co-expression","Functional interaction","Functional interaction","Co-expression","Hotspot-Target","Hotspot-Target","Co-expression","Functional interaction","Co-expression","Functional interaction"))



# relations <- relations[-c(5:9),]
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



