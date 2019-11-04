#Paquetes requeridos
require(rWikiPathways)
require(biomaRt)


print('Importing House Keeping Genes ...')
data("hkGenes")
ensembl = useEnsembl(biomart = "ensembl", dataset="hsapiens_gene_ensembl")
hkGenes<- getBM(attributes=c('ensembl_gene_id','hgnc_symbol'), filters ='ensembl_gene_id', values =as.vector(hkGenes), mart = ensembl)
hkGenes<-as.vector(unique(hkGenes$hgnc_symbol))   ## 565 HKGenes

#Cell cycle genes
print('Importing Cell Cycle Genes ...')
wiki_cc<-getXrefList(pathway = 'WP179', systemCode = 'H')    ## 120 wikicc_genes

#Breast Cancer genes
print('Importing BRCA Genes ...')
wiki_brca<-getXrefList(pathway = 'WP1984', systemCode = 'H')  ##Wikipathways

#Lung Adenocarcinome  WP2512
print('Importing LUAD Genes ...')
wiki_luad<-getXrefList(pathway = 'WP2512', systemCode = 'H')  #Wikipathways

genes = unique(c(hkGenes, wiki_cc, wiki_brca, wiki_luad))
print(paste('Total imported:',length(genes), 'genes', sep = ' '))
