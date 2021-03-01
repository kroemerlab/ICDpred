#' ICD inducing prediction
#'
#' Compute an ICD score
#'
#' @import rcdk
#' @import ChemmineR
#' @import RCurl
#'
#' @param CID Pubchem CID#'
#' @param SDF Molecule 2D structure (must be sdf format)
#' @return A list with computed scores and other useful information
#'
#' @author Allan Sauvat, \email{allan.sauvat@gustaveroussy.fr}
#'
#' @examples
#' test=ICDscoring(CID=4212) #Compute score for MTX
#'
#' @export
#'

ICDscoring = function(CID,SDF){
  
  #Error handling
  if(missing(SDF)){SDF=getIds(CID)}
  if(missing(CID)){CID=as.numeric(sdfid(SDF))}
  if(length(CID)!=1|length(SDF)!=1){
    return('Function can only be applied on a single molecule at a time, please iterate!')
  }
  #IMPORT MOLECULAR INFOS
  SMI=sdf2smiles(SDF)
  PAR=parse.smiles(unname(as.character(SMI)))[[1]]
  
  ## Need to configure the molecule
  invisible({do.aromaticity(PAR);do.typing(PAR);do.isotopes(PAR)})
  
  #CALCULATE DESCRIPTORS
  dn=sapply(get.desc.categories(),function(x)get.desc.names(x))
  dn$constitutional=dn$constitutional[-7] #Buggy descriptor to be removed
  descs=do.call('cbind',lapply(dn, function(x)eval.desc(PAR,x)))
  descs[,'constitutional.nAtomLAC']=NA #Add buggy descriptor name
  
  #ADD PRECOMPUTED DESCRIPTORS
  miss=c("Molecular.Weight","Hydrogen.Bond.Donor.Count","Hydrogen.Bond.Acceptor.Count",
         "Rotatable.Bond.Count","Exact.Mass","Monoisotopic.Mass","Topological.Polar.Surface.Area","Heavy.Atom.Count","Complexity")
  
  PUB=any(grepl('PUBCHEM_CACTVS_COMPLEXITY',colnames(datablock2ma(datablocklist=datablock(SDF)))))
  
  if(PUB){
    nm=c("PUBCHEM_MOLECULAR_WEIGHT","PUBCHEM_CACTVS_HBOND_DONOR","PUBCHEM_CACTVS_HBOND_ACCEPTOR","PUBCHEM_CACTVS_ROTATABLE_BOND",
         "PUBCHEM_EXACT_MASS","PUBCHEM_MONOISOTOPIC_WEIGHT","PUBCHEM_CACTVS_TPSA","PUBCHEM_HEAVY_ATOM_COUNT","PUBCHEM_CACTVS_COMPLEXITY")
    
    precp=data.frame(datablock2ma(datablocklist=datablock(SDF)))[nm]
    precp=sapply(colnames(precp),function(x)precp[,x]=as.numeric(as.character(precp[,x])))
  }else{
    print('No CACTVS records, calculating an approximation...')
    precp=cbind(get.natural.mass(PAR),descs[,c('electronic.nHBDon','electronic.nHBAcc','constitutional.nRotB')],
                get.exact.mass(PAR),descs[c('constitutional.MW','topological.TopoPSA')],
                descs$constitutional.nAtom-get.total.hydrogen.count(PAR),
                exp(log(descs[,'topological.fragC']+1)*cdk2cac[2]+cdk2cac[1])-1) #CACTVS complexity approximation
  }

  descs[miss]=precp #For compatibility with ancient use of rpubchem
  
  #CLEAN DATA
  descs=descs[setdiff(colnames(Train.Descs),Txt)] #Keep relevant features
  if(any(is.na(descs))){
    nm = colnames(descs)[is.na(descs)]
    for(x in nm){descs[x]=median(Train.Descs[,nm],na.rm=T)}
  }
  
  #EVALUATE SCORE
  Score=as.numeric(predict(Best.mod$Model,data.frame(pred.transform(as.matrix(descs)))[1:(Best.mod$ncp)]))
  
  #RETURN RESULTS
  ccl='unable to calculate score'
  if(!is.na(Score)){
    if(Score<2){
      ccl='Weak chance of inducing ICD'
    }else if(Score>6){
      ccl='Artefactual risk'
    }else{
      ccl='Potential ICD inducer'
    }}
  return(list(cid=CID,Descriptors=descs[setdiff(colnames(Train.Descs),Txt)],Score=Score,Conclusion=ccl,sdf=SDF))
  
}

