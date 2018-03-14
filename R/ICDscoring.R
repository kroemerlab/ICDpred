#' ICD inducing prediction
#'
#' Compute an ICD score
#'
#' @import rcdk
#' @import ChemmineR
#' @import RCurl
#'
#' @param CID Pubchem CID#'
#' @return A list with computed scores and other useful information
#'
#' @author Allan Sauvat, \email{allan.sauvat@gustaveroussy.fr}
#'
#' @examples
#' #Not running
#' test=ICDscoring(CID=4212) #Compute score for MTX
#'
#' @export
#'

ICDscoring = function(CID){
  
  #IMPORT MOLECULAR INFOS
  SDF=getIds(CID)
  SMI=sdf2smiles(SDF)
  PAR=parse.smiles(as.character(SMI))
  
  #CALCULATE DESCRIPTORS
  dn=sapply(get.desc.categories(),function(x)get.desc.names(x))
  descs=do.call('cbind',lapply(dn, function(x)eval.desc(PAR,x)))
  
  #ADD PRECOMPUTED DESCRIPTORS  
  miss=c("Molecular.Weight","Hydrogen.Bond.Donor.Count","Hydrogen.Bond.Acceptor.Count",
         "Rotatable.Bond.Count","Exact.Mass","Monoisotopic.Mass","Topological.Polar.Surface.Area","Heavy.Atom.Count","Complexity")
  pugn=gsub('[.]','',miss);pugn[c(2,3,7)]=c('HBondDonorCount','HBondAcceptorCount','TPSA')
  precp=sapply(pugn,function(x)as.numeric(getURL(paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',as.character(CID),'/property/',x,'/TXT'))))
  descs[miss]=precp #For compatibility with ancient use of rpubchem
  
  #EVALUATE SCORE
  Score=predict(Best.mod$Model,data.frame(pred.transform(as.matrix(descs[setdiff(colnames(Train.Descs),Txt)])))[1:(Best.mod$ncp)])
  
  #RETRUN RESULTS
  if(Score<2){
    ccl='Weak chance of inducing ICD'
  }else if(Score>6){
    ccl='Artefactual risk'
  }else{
    ccl='Potential ICD inducer'
  }
  return(list(Descriptors=descs[setdiff(colnames(Train.Descs),Txt)],Score=Score,Conclusion=ccl,sdf=SDF))
  
}

