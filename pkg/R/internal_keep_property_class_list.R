internal_keep_property_class_list<-
    function(keep_list_name=c("property","class")){
        if(keep_list_name=="property"){
            keep_list<-
                list(c("term"
                       ,"db"
                       ,"id"
                       ,"positionStatus"
                       ,"sequencePosition")
                     ,c("term"
                        ,"xref"
                        ,"evidenceCode"
                        ,"featureLocation"
                        ,"modificationType"
                        ,"sequenceIntervalBegin"
                        ,"sequenceIntervalEnd"))
            
        }else if(keep_list_name=="class"){
            keep_list<-
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
        }

        
        return(keep_list)
    }