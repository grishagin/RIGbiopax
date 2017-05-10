internal_keep_class_list<-
    function(){
        keep_class_list<-
            list(c("Protein"
                   ,"SmallMolecule"
                   ,"PhysicalEntity"
                   ,"RNA"
                   ,"DNA")
                 ,"Complex"
                 # ,c("BiochemicalReaction"
                 #    ,"TemplateReaction"
                 #    ,"ComplexAssembly"
                 #    ,"Transport"
                 #    ,"TransportWithBiochemicalReaction"
                 #    ,"Degradation")
                 # ,c("Catalysis"
                 #    ,"TemplateReactionRegulation" 
                 #    ,"Control"               
                 #    ,"Modulation")
                 )
        
        return(keep_class_list)
    }