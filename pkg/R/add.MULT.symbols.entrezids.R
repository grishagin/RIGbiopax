add.MULT.symbols.entrezids <-
function(df_pw_proteins=NULL
             ,keytypes=c("entrezgene"
                         ,"mim"
                         ,"uniprot"
                         ,"unigene"
                         ,"ensemblgene"
                         ,"hprd"
                         ,"hgnc")
             ,filter_keytypes=TRUE){
        
        require(dplyr)
        
        if (is.null(df_pw_proteins)){
            stop("add_symbols_entrezids: df_pw_proteins is missing!")
        }
        
        #replace db ids name with proper annotation ids
        df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type %chin% "LL"]<-"entrezgene"
        df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type %chin% "EntrezGene"]<-"entrezgene"
        df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type %chin% "Entrez Gene"]<-"entrezgene"
        df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type %chin% "UniProt"]<-"uniprot"
        df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type %chin% "UniGene"]<-"unigene"
        df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type %chin% "Ensembl"]<-"ensemblgene"
        df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type %chin% "OMIM"]<-"mim"
        df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type %chin% "HPRD"]<-"hprd"
        df_pw_proteins$biopax.Gene.ID.Type[df_pw_proteins$biopax.Gene.ID.Type %chin% "HGNC"]<-"hgnc"
        
        
        #ensure correct types of Gene Symbol and Entrezid columns
        df_pw_proteins$biopax.Gene.Symbol<-
            as.character(df_pw_proteins$biopax.Gene.Symbol)
        
        df_pw_proteins$ENTREZID<-
            as.integer(df_pw_proteins$ENTREZID)
        
        #if there are multiple comma|semicolon|pipe-separated values,
        #split them and expand the df
        df_pw_proteins<-
            df_pw_proteins %>%
            expand.df.via.split.col(colToSplit="biopax.Gene.ID")
        
        #annotate data table iteratively
        for(KEYTYPE in keytypes){
            df_pw_proteins<-
                df_pw_proteins %>%
                add.symbols.entrezids.mygene(KEYTYPE)
        }
        
        #if multiple comma|semicolon|pipe-separated values were added when annotating,
        #split them and expand the df
        df_pw_proteins<-
            df_pw_proteins %>%
            expand.df.via.split.col(colToSplit=c("ENTREZID","biopax.Gene.Symbol"))
        
        
        if(filter_keytypes){
            df_pw_proteins<-
                df_pw_proteins %>%
                filter(biopax.Gene.ID.Type %chin% keytypes)
        }
        
        #return the annotated df
        return(df_pw_proteins)
        
    }
