This project also used other general scripts to generate an unfolded SFS using EstSFS output- these are available in GeneralScripts.
This project also used GeneralScripts/AnnotateRef_CodonDegeneracy.py for site annotation.

Details of per-gene estimates of pi0 and pi0/pi4 calculation are available in CapsellaGenomeAnalysis.R and ArabidopsisGenomeAnalysis.R

Primary figures filename records:

Grouping folder:
Files for Arabidopsis thaliana

GroupingCapsella folder:
Files for Capsella grandiflora

Data for the per-gene estimates of genome biology features used in this analysis, plus per gne estimates of pi0 (piN) and pi4 (piS) are:
Pi_Biology_byGene_Athaliana.csv and Pi_Biology_byGene_Capsella.csv.
These files were written using the steps shown in ArabidopsisGenomeAnalysis.R and CapsellaGenomeAnalysis.R

###Note on names
Original names
A thaliana: CoeffVarfpkm, connec, express, NetworkConnectivity, fpkm, GOTermCount (PolyDFE10 only)

C grandiflora: aranetconnec, aranetconneccorrec, express, norm.exp, GOTermCount (PolyDFE10 only)

Renamed for downstream analysis:
A thaliana: CoeffVarFpkm, conneccorrec, expresscorrec, NetworkConnectivity, fpkm, GOTermCount (PolyDFE10 only)

C grandiflora: NetworkConnectivity, conneccorrec, expresscorrec, NormExpression, GOTermCount (PolyDFE10 only)



