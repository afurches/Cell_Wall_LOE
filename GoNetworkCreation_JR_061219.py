#script to create GO-term functional networks
#created by Jonathon Romero
#as used in the manuscript:
#Furches, A., Kainer, D., Weighill, D., Large, A., Jones, P., Walker, A.M., Romero, J., Gazolla, JGFM,
#    Joubert, W., Shah, M., Streich, J., Ranjan, P., Schmutz, J., Sreedasyam, A., Macaya-Sanz, D., Zhao, N.,
#    Martin, M.Z., Rao, X., Dixon, R.A., DiFazio, S., Tschaplinski, T.J., Chen, J-C., Tuskan, G.A., Jacobson, D.
#    Finding New Cell Wall Regulatory Genes in Populus trichocarpa Using Multiple Lines of Evidence.
#    Front. Plant Sci. (in review).



import numpy as np
import pandas as pd

longanno=pd.read_csv("PlantRegMap_Ptr_GO_annotation", delimiter='\t')
GoList=pd.read_csv("cellWallGo.txt", delimiter='\t',header=None)

GoSet=set(longanno.iloc[:,1].values.tolist())
GeneSet=set(longanno.iloc[:,0].values.tolist())
GeneList=list(GeneSet)

adjGenes=pd.DataFrame(np.identity(len(GeneSet)),columns=GeneSet,index=GeneSet)

#make Specific GO adj matrix
count=0
GoCount=0
GoNum=[]
for Go in GoList[0].values.tolist():
    tempGenes=longanno[longanno.iloc[:,1] == Go]
    tempGeneNum=len(tempGenes)
    GoNum.append(tempGeneNum)
    count+=1
    print(count)
    if tempGeneNum > 1:
        GoCount+=1
        print("GoCount:",GoCount, "GoID:", Go, "Number of Genes:", tempGeneNum)
        #if tempGeneNum > 1000:
        #    continue
        geneList=tempGenes.iloc[:,0].values.tolist()
        goValue=(1.0/(tempGeneNum-1))
        for outerGene in range(len(geneList)):
            for innerGene in range(outerGene,len(geneList)):
                if adjGenes.loc[geneList[innerGene],geneList[outerGene]] < goValue:
                    adjGenes.at[geneList[innerGene],geneList[outerGene]]=goValue
                    adjGenes.at[geneList[outerGene],geneList[innerGene]]=goValue
                    
GoEdge=adjGenes.rename_axis("Source").reset_index().melt(id_vars="Source", var_name="Sink",value_name="Weight")
temp=GoEdge[GoEdge['Weight']!=0]
temp[temp['Weight']!=1].to_csv("pTryCellWallGo.Network.csv",index=False)
