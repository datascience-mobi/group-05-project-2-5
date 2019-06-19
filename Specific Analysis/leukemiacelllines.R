
FCbiomarkerswithtissue <- as.data.frame(FCbiomarkerswithtissue)
Leukemiacelllines <- colnames(FCbiomarkerswithtissue[Metadatavorinostattissue== "Leukemia"])
FCleukemia <- FCbiomarkerswithtissue[,Leukemiacelllines]