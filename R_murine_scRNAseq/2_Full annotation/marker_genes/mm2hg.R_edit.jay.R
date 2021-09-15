##Convert mouse genes to human orthologs
#SWK 09/28/15

####SMP EDIT FOR CONVERSION WHEN SERVERS ARE DOWN 11/13/15

#currently translates gene symbols to gene symbols
#WILL IMPLEMENT ENSEMBL ID EVENTUALLY

require(biomaRt)
listMarts(host="www.ensembl.org")

mouse=useMart("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl",host="www.ensembl.org")
human=useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",host="www.ensembl.org")

mm2hg = function(gene.list=NULL){
  mm_symbols = as.vector(gene.list)
  
  hg_symbols = getLDS(attributes = 'external_gene_name', filters = 'external_gene_name', values = mm_symbols, mart = mouse,
                      attributesL = 'external_gene_name', martL = human)
  
  return(unique(hg_symbols[,2]))
}

hg2mm = function(gene.list=NULL){
  hg_symbols = as.vector(gene.list)
  
  mm_symbols = getLDS(attributes = 'external_gene_name', filters = 'external_gene_name', values = hg_symbols, mart = human,
                      attributesL = 'external_gene_name', martL = mouse)
  
  return(unique(mm_symbols[,2]))
}

