import pandas as pd
from tabulate import tabulate
import os

'''
Read in values of number of reads processed in total, number of discared reads after processeing,
and number of successfully processed reads (total minus discarded)
'''
# df=pd.read_csv(snakemake.input[0], index_col=False, names=["processed"]).div(2)
df=pd.read_csv(snakemake.input[1], index_col=False, names=["processed"]).div(2)
df_list = list(dict.fromkeys(list(df.processed)))
total_processed = df_list[0]
discarded = df_list[1]
discarded_perc= round(((df_list[1]/df_list[0])*100), 2)

print("HERE")
print(df_list)
successfull = df_list[2]
successfull_perc= round(((df_list[2]/df_list[0])*100), 2)

libraryName= os.path.basename(snakemake.input[0])

libraryName=libraryName.replace('_read_count.tsv','')

'''
List of possible barcodes
'''
df3=pd.read_csv(snakemake.input[4], index_col=False, names=['barcode_combination'])



'''
Read counts of reads assigned to each barcode following demultiplexing with flexbar
'''
# df2=pd.read_csv(snakemake.input[1], index_col=False, names=["PE_read_count"])
df2=pd.read_csv(snakemake.input[0], index_col=False, names=["PE_read_count"])
df2['%_of_processed']=df2.PE_read_count.div(df_list[2])
df2['%_of_processed']=df2['%_of_processed'] * 100
df2=df2.round(2)
df2['%_of_processed']=df2['%_of_processed'].astype(str) + '%'
df2 = pd.concat([df3, df2], axis=1).sort_values(by ='barcode_combination' ).reset_index(drop=True)
print(df2)


df4=pd.read_csv(snakemake.input[3], delim_whitespace=True, index_col=False, skip_blank_lines=True, names=['barcode_combination', 'individualID'], converters={'individualID' : str})
result = pd.merge(df3,df4, how='outer').fillna('none').drop_duplicates(subset=['barcode_combination'], keep='first').reset_index(drop=True).sort_values(by ='barcode_combination').reset_index(drop=True)
print(result)#

result= pd.concat([result, df2['PE_read_count'], df2['%_of_processed']], axis=1)


renamingLoc = pd.read_csv(snakemake.input[2], index_col=False, names=["location"])
renamingLoc = list(dict.fromkeys(list(renamingLoc.location)))
renamingLocation = renamingLoc[0]
runNum=renamingLoc[1]

renamingR1 = result[['barcode_combination', 'individualID']].copy()
renamingR2 = result[['barcode_combination', 'individualID']].copy()
renamingR1['barcode_combination'] = renamingLocation + libraryName + '_barcode_' + renamingR1['barcode_combination'].astype(str) + '_1.fastq'
renamingR2['barcode_combination'] = renamingLocation + libraryName + '_barcode_' + renamingR2['barcode_combination'].astype(str) + '_2.fastq'
renamingR1['individualID'] = renamingLocation + libraryName + '_' + renamingR1['individualID'].astype(str) + '_' + runNum + '_1.fastq'
renamingR2['individualID'] = renamingLocation + libraryName + '_' + renamingR2['individualID'].astype(str) + '_' + runNum + '_2.fastq'

renamingR1.to_csv(snakemake.output[0], header=False, sep='\t', index=False)
renamingR2.to_csv(snakemake.output[1], header=False, sep='\t', index=False)

total_processed=str(total_processed)
total_processed=total_processed.replace(".0",'')
discarded=str(discarded)
discarded=discarded.replace(".0",'')
successfull=str(successfull)
successfull=successfull.replace(".0",'')

libraries=['#' + libraryName]
totalP= ["Total Reads Processed (Paired):        " + total_processed + "   ( 100 %)", "Discarded reads (Paired):              " + discarded + "    ( "+str(discarded_perc)+"%)", "Successfully Processed reads (Paired): " + successfull + "   ( "+str(successfull_perc)+"%)"]
loadinTop=pd.DataFrame(totalP, columns=['#' + libraryName])

with open(snakemake.output[2], 'w') as outputfile:
#    print('#' + libraryName, file=outputfile)
#    print("Total Reads Processed (Paired):        " + total_processed + "   ( 100 %)", file=outputfile)
#    print("Discarded reads (Paired):              " + discarded + "    ( "+str(discarded_perc)+"%)", file=outputfile)
#    print("Successfully Processed reads (Paired): " + successfull + "   ( "+str(successfull_perc)+"%)", file=outputfile)
    print(tabulate(loadinTop, headers='keys',tablefmt="rst", showindex=False), file=outputfile)
#    print('', file=outputfile)
    print(tabulate(result, headers='keys', tablefmt="psql", showindex=False), file=outputfile)
    print('', file=outputfile)
