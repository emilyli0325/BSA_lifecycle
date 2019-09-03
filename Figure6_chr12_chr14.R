#### chr12 and chr14

G.A <- read.csv("df_filt42.csv", header = TRUE)
G.B <- read.csv("df_filt44.csv", header = TRUE)

G.A.chr12 <- G.A[G.A$CHROM == "12", ]
G.B.chr12 <- G.B[G.B$CHROM == "12", ]


G.A.chr14 <- G.A[G.A$CHROM == "14", ]
G.B.chr14 <- G.B[G.B$CHROM == "14", ]

pChr12 <- ggplot(data = G.A.chr12,aes(x = POS, y = G)) +scale_x_continuous(breaks = seq(from = 0,to = max(G.A.chr12$POS), by = 10^(floor(log10(max(G.A.chr12$POS))))), name = "Genomic Position (Mb)") +ylab("G'")+ geom_point(aes_string(x = "POS", y = "G"),color = "grey",size=2)+ geom_point(data = G.B.chr12, aes_string(x = "POS", y = "G"),color = "orange",size=2) +geom_line(aes_string(x = "POS", y = "Gprime"),color = "red",size=0.8) + geom_line(data = G.B.chr12, aes_string(x = "POS", y = "Gprime"),color = "black",size=0.8) +theme_classic()+geom_vline(xintercept = 1102148, linetype="dashed", color = "black", size=0.8)+geom_vline(xintercept = 1327968, linetype="dashed", color = "black", size=0.8)


pChr14 <- ggplot(data = G.A.chr14,aes(x = POS, y = G)) +scale_x_continuous(breaks = seq(from = 0,to = max(G.A.chr12$POS), by = 10^(floor(log10(max(G.A.chr12$POS))))), name = "Genomic Position (Mb)") +ylab("G'") + geom_point(aes_string(x = "POS", y = "G"),color = "grey",size=2)+ geom_point(data = G.B.chr14, aes_string(x = "POS", y = "G"),color = "orange",size=2) +geom_line(aes_string(x = "POS", y = "Gprime"),color = "red",size=0.8)+ geom_line(data = G.B.chr14, aes_string(x = "POS", y = "Gprime"),color = "black",size=0.8)+theme_classic()+geom_vline(xintercept = 2378002, linetype="dashed", color = "black", size=0.8)+geom_vline(xintercept = 2541869, linetype="dashed", color = "black", size=0.8)

library(gridExtra)
pdf('Figure6_chr12_chr14_01152019.pdf', width=6, height=6)
grid.arrange(pChr12,pChr14,
widths = c(2.5, 1.2),heights = c(1,1), layout_matrix = rbind(c(1,NA),c(2,2))
)
dev.off()

