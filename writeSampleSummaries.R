#!/nfs/goldstein/software/R-3.0.1/bin/Rscript
library(rtf)

writeDNM <- function(dnm,rtf){
    for(i in 1:dim(dnm)[1]){
        setFontSize(rtf,12)
        gene = gsub("'","",dnm[i,]$Gene.Name)
        addHeader(rtf,gene,subtitle=dnm[i,]$Variant.ID,font.size=14)
        startParagraph(rtf)
        addText(rtf,paste0("This is a ",dnm[i,]$Flag," ",gsub("_"," ",tolower(dnm[i,]$Function))," variant in ",gene,". "))
        if(dnm[i,]$Ctrl.MAF == 0 & (dnm[i,]$Evs.All.Maf == 0 | is.na(dnm[i,]$Evs.All.Maf)) & (dnm[i,]$ExAC.global.maf == 0 | is.na(dnm[i,]$ExAC.global.maf))){
            addText(rtf,"This variant is absent from internal and external control samples. ")}
        else{addText(rtf,paste0("This variant has a control MAF of ",dnm[i,]$Ctrl.MAF*100,"% in IGM controls, ",dnm[i,]$Evs.All.Maf*100,"% in EVS, and ",dnm[i,]$ExAC.global.maf,"% in ExAC. "))}
        if(dnm[i,]$Function == "NON_SYNONYMOUS_CODING"){addText(rtf,paste0("It is a ",gsub("_"," ",dnm[i,]$Polyphen.Humvar.Prediction)," missense variant with a PolyPhen2 score of ",dnm[i,]$Polyphen.Humvar.Score,". "))}
        adj=" ";if(!is.na(dnm[i,]$X0.05._anypopn_RVIS.tile.ExAC.) & dnm[i,]$X0.05._anypopn_RVIS.tile.ExAC. <= 25){adj="n in"}
        addText(rtf,paste0(gene," is a",adj,"tolerant gene with an RVIS score of ",signif(as.double(dnm[i,]$X0.05._anypopn_RVIS.tile.ExAC.),digits=2),". "))
        if(!is.na(dnm[i,]$LoF.pLI.ExAC.) & dnm[i,]$LoF.pLI.ExAC. >= .9){addText(rtf,paste0(gene," is LoF intolerant gene with a PLI score of ",signif(as.double(dnm[i,]$LoF.pLI.ExAC.),digits=2),". "))}
        if(is.na(dnm[i,]$Gerp.RS.Score)){adj=" not"}else{adj=" very strongly";if(dnm[i,]$Gerp.RS.Score<5){adj=" strongly"};if(dnm[i,]$Gerp.RS.Score<4){adj=""};if(dnm[i,]$Gerp.RS.Score<2){adj=" weakly"};if(dnm[i,]$Gerp.RS.Score<0){adj=" not"}}
        addText(rtf,paste0("This site is",adj," conserved with a GERP++ RS score of ",dnm[i,]$Gerp.RS.Score,". "))
        if(!is.na(dnm[i,]$OMIM.Disease.Name)){
           adj="";if(dnm[i,]$MGI.Essential){adj=", and an essential gene"}
           addText(rtf,paste0(gene," is an OMIM disease gene associated with ",gsub(" \\|",", and", dnm[i,]$OMIM.Disease.Name),adj,". "))}
        else if(dnm[i,]$MGI.Essential){addText(rtf,paste0(gene," is an essential gene. "))}
        if(!is.na(dnm[i,]$ClinVar.Clinical.Significance)){
           addText(rtf,paste0("This variant is listed as ",tolower(dnm[i,]$ClinVar.Clinical.Significance)," in ClinVar. "))}
        if(!is.na(dnm[i,]$HGMD.Variant.Class)){
            addText(rtf,paste0("This variant is listed as ",dnm[i,]$HGMD.Variant.Class," in HGMD. "))}
        addText(rtf,"\n")
        endParagraph(rtf)
    }
}  

writeCHET <- function(chet,rtf){
    for(i in 1:dim(chet)[1]){
        setFontSize(rtf,12)
        gene = gsub("'","",chet[i,]$Gene.Name)
        addHeader(rtf,gene,paste0(subtitle=chet[i,]$Variant.ID.1," , ",subtitle=chet[i,]$Variant.ID.2),font.size=14)
        startParagraph(rtf)
        addText(rtf,paste0("These are ",chet[i,]$Flag," variants in ",gene,". The first is a ",gsub("_"," ",tolower(chet[i,]$Function.1))," and the second is a ",gsub("_"," ",tolower(chet[i,]$Function.2)), " variant. "))
        adj=" ";if(!is.na(chet[i,]$X0.05._anypopn_RVIS.tile.ExAC..1) & chet[i,]$X0.05._anypopn_RVIS.tile.ExAC..1 <= 25){adj="n in"}
        addText(rtf,paste0(gene," is a",adj,"tolerant gene with an RVIS score of ",signif(as.double(chet[i,]$X0.05._anypopn_RVIS.tile.ExAC..1),digits=2),". "))
        if(!is.na(chet[i,]$LoF.pLI.ExAC..1) & chet[i,]$LoF.pLI.ExAC..1 >= .9){addText(rtf,paste0(gene," is LoF intolerant gene with a PLI score of ",signif(as.double(chet[i,]$LoF.pLI.ExAC..1),digits=2),". "))}
        if(!is.na(chet[i,]$OMIM.Disease.Name.1)){
           adj="";if(chet[i,]$MGI.Essential){adj=", and an essential gene"}
           addText(rtf,paste0(gene," is an OMIM disease gene associated with ",gsub(" \\|",", and", chet[i,]$OMIM.Disease.Name.1),adj,". "))}
        else if(chet[i,]$MGI.Essential){addText(rtf,paste0(gene," is an essential gene. "))}
        addText(rtf,"\n")
        if(chet[i,]$Ctrl.MAF.1 == 0 & (chet[i,]$Evs.All.Maf.1 == 0 | is.na(chet[i,]$Evs.All.Maf.1)) & (chet[i,]$ExAC.global.maf.1 == 0 | is.na(dnm[i,]$ExAC.global.maf))){
            addText(rtf,"The first variant is absent from internal and external control samples. ")}
        else{addText(rtf,paste0("The first variant has a control MAF of ",chet[i,]$Ctrl.MAF.1*100,"% in IGM controls, ",chet[i,]$Evs.All.Maf.1*100,"% in EVS, and ",chet[i,]$ExAC.global.maf.1,"% in ExAC. "))}
        if(chet[i,]$Function.1 == "NON_SYNONYMOUS_CODING"){addText(rtf,paste0("It is a ",gsub("_"," ",chet[i,]$Polyphen.Humvar.Prediction.1)," missense variant with a PolyPhen2 score of ",chet[i,]$Polyphen.Humvar.Score.1,". "))}
        if(is.na(chet[i,]$Gerp.RS.Score.1)){adj=" not"}else{adj=" very strongly";if(chet[i,]$Gerp.RS.Score.1<5){adj=" strongly"};if(chet[i,]$Gerp.RS.Score.1<4){adj=""};if(chet[i,]$Gerp.RS.Score.1<2){adj=" weakly"};if(chet[i,]$Gerp.RS.Score.1<0){adj=" not"}}
        addText(rtf,paste0("This site is",adj," conserved with a GERP++ RS score of ",chet[i,]$Gerp.RS.Score.1,". "))
        if(!is.na(chet[i,]$ClinVar.Clinical.Significance.1)){
           addText(rtf,paste0("This variant is listed as ",tolower(chet[i,]$ClinVar.Clinical.Significance.1)," in ClinVar. "))}
        if(!is.na(chet[i,]$HGMD.Variant.Class.1)){
            addText(rtf,paste0("This variant is listed as ",chet[i,]$HGMD.Variant.Class.1," in HGMD. "))}
        addText(rtf,"\n")
        if(chet[i,]$Ctrl.MAF.2 == 0 & (chet[i,]$Evs.All.Maf.2 == 0 | is.na(chet[i,]$Evs.All.Maf.2)) & (chet[i,]$ExAC.global.maf.2 == 0 | is.na(dnm[i,]$ExAC.global.maf))){
            addText(rtf,"The second variant is absent from internal and external control samples. ")}
        else{addText(rtf,paste0("The second variant has a control MAF of ",chet[i,]$Ctrl.MAF.2*100,"% in IGM controls, ",chet[i,]$Evs.All.Maf.2*100,"% in EVS, and ",chet[i,]$ExAC.global.maf.2,"% in ExAC. "))}
        if(chet[i,]$Function.2 == "NON_SYNONYMOUS_CODING"){addText(rtf,paste0("It is a ",gsub("_"," ",chet[i,]$Polyphen.Humvar.Prediction.2)," missense variant with a PolyPhen2 score of ",chet[i,]$Polyphen.Humvar.Score.2,". "))}
        if(is.na(chet[i,]$Gerp.RS.Score.2)){adj=" not"}else{adj=" very strongly";if(chet[i,]$Gerp.RS.Score.2<5){adj=" strongly"};if(chet[i,]$Gerp.RS.Score.2<4){adj=""};if(chet[i,]$Gerp.RS.Score.2<2){adj=" weakly"};if(chet[i,]$Gerp.RS.Score.2<0){adj=" not"}}
        addText(rtf,paste0("This site is",adj," conserved with a GERP++ RS score of ",chet[i,]$Gerp.RS.Score.2,". "))
        if(!is.na(chet[i,]$ClinVar.Clinical.Significance.2)){
           addText(rtf,paste0("This variant is listed as ",tolower(chet[i,]$ClinVar.Clinical.Significance.2)," in ClinVar. "))}
        if(!is.na(chet[i,]$HGMD.Variant.Class.2)){
            addText(rtf,paste0("This variant is listed as ",chet[i,]$HGMD.Variant.Class.2," in HGMD. "))}
        addText(rtf,"\n")
        addText(rtf,"\n")
        endParagraph(rtf)
    }
}   

writeSummary <- function(dnm,hom,hem,chet,dir){
    dir.create(file.path(dir,"Sample_Summaries"),showWarnings=F)
    allChilds <- c(as.vector(dnm$Child),as.vector(hom$Child),as.vector(hem$Child),as.vector(chet$Child))
    if(length(allChilds) == 0){return}
    out <- file.path(dir,"Sample_Summaries",paste0(allChilds[1],".doc"))
    rtf <- RTF(out,width=8.5,height=11,omi=c(1,1,1,1),font.size=18)
    addHeader(rtf,allChilds[1],font.size=18)
    if(dim(dnm)[1] > 0){writeDNM(dnm,rtf)}
    if(dim(hem)[1] > 0){writeDNM(hem,rtf)}
    if(dim(hom)[1] > 0){writeDNM(hom,rtf)}
    if(dim(chet)[1] > 0){writeCHET(chet,rtf)}
    done(rtf)
    return
}

