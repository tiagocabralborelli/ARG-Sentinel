def pca (df):
    """-> Criar um dataframe para plotar um PCA"""
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA
    import pandas as pd

    del df.index.name
    df.reset_index(inplace = True)
    df.rename(columns = {'index':'Grupos:'}, inplace = True)
    grupos = df['Grupos:']
    df.drop('Grupos:', inplace = True, axis = 1)
    genes = df.columns
    df['Grupos:'] = grupos
    features = list(genes)
    #Separar as caracteristicas
    x = df.loc[:, features].values
    #Separar os alvos
    y = df.loc[:,['Grupos:']].values
    #Normalização das caracteristicas
    x = StandardScaler().fit_transform(x)
    #PCA
    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(x)
    principalDf = pd.DataFrame(data = principalComponents, columns = ['PC1','PC2'])
    finalDf = pd.concat([principalDf, df[['Grupos:']]], axis = 1)
    finalDf.drop(0, axis = 0, inplace = True)                            
    
    
    return finalDf
def extract (df,sp,gene,inp):
    from Bio import SeqIO
    from Bio.Seq import Seq
    df = df.loc[sp,:]
    fasta = open(f'{sp}.fasta','a')
    for items in df.itertuples():
        if gene in items.gene:
            for f in SeqIO.parse(f'{inp}\\{items.file}','fasta'):
                if items.sequence == f.id:
                    fasta.write(f'>{items.gene}_{f.id}\n')
                    if items.strand == '+':
                        fasta.write(f'{f.seq[items.start-1:items.end+1]}\n')
                    else:
                        fasta.write(f'{f.seq[items.start-1:items.end+1].reverse_complement()}\n')
    fasta.close()
def parsnp(df,org, dirr, name):
    "Copia os genomas de uma mesma especies para um diretorio separado e executa o parsnp"
    import os
    import pandas as pd
    os.system(f'mkdir {dirr}/{name}')
    for c in df.itertuples():
        if c.specie == org:
            os.system(f'cp {c.path} {dirr}/{name}')
def convert_to_fasta(path):
    """Convert genbanks files into fasta
    ->Pam path: genbank files path""" 
    from Bio import SeqIO
    import glob
    gbfiles = glob.glob(f"{path}/*.gbff")
    for c in gbfiles:
        f = SeqIO.convert(f"{c}","genbank",f"{c}.fasta","fasta")
def SwitchMonths(x):
    """Troca os meses pelos seus respectvios números.
    ->Pam x: Dados para serem alterados """
    
    #Listas para fazer a troca
    monthsnum = ["01","02","03","04","05","06","07","08","09","10","11","12"]
    monthsstr = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
    xsplit = x.split("-")
    setx = set(xsplit)
    set_months = set(monthsstr)
    if len(setx & set_months) > 0:
        InComum = list(setx & set_months) 
        return x.replace(InComum[0],monthsnum[monthsstr.index(InComum[0])])
    else:
        return x
def cataloger(path):
    """Cria um catálogo com informações adicionais dos genomas no formato gbff para adicionar nos bancos de dados.
    -> Pam path: diretório onde os arvivos estão"""
    
    #Complexidade:  O(n³) pelo menos. É possível reduzir?
        
    #Dependências
    from Bio import SeqIO
    import glob
    import pandas as pd
    
    #Listas para criar o dataframe 
    accession = [] #Número de acesso
    colecdate = [] #Data da coleta
    host = []      #Hospedeiro de onde a bactéria foi isolada 
    source = []    #Fonte
    coord = []     #Coordenadas
    country = []   #Local 
    strain = []  
    plasmids = []
    organism = []  #Espécie
    
    gbfiles = glob.glob(path+"/*.gbff") #Pegar os arquivos no diretório dado pelo parâmetro path
    for files in gbfiles:
        for record in SeqIO.parse(f"{files}","genbank"): #Abre os arquivos gbff para retirar o resto das informações
            for f in record.features:
                if f.type == 'source':
                    
                    #Atribuir número de acesso
                    accession.append(files.split("/")[-1])
    
                     
                    #Atribuir data de coleta
                    try:
                        colecdate.append(f.qualifiers["collection_date"][0])
                    except:
                        colecdate.append("Na")
                    
                    #Atribuir hospedeiro
                    try:
                        host.append(f.qualifiers["host"][0])
                    except:
                        host.append("Na")
                        
                    #Atribuir fonte de coleta
                    try:
                        source.append(f.qualifiers["isolation_source"][0])
                    except:
                        source.append("Na")
                        
                    #Atribuir coordenadas da coleta
                    try:
                        coord.append(f.qualifiers["lat_lon"][0])
                    except:
                        coord.append("Na")
                        
                    #Atribuir local da coleta
                    try:
                        country.append(f.qualifiers["country"][0])
                    except:
                        country.append("Na")
                        
                    #Atribuir espécie
                    try:
                        organism.append(f.qualifiers["organism"][0])
                    except:
                        organism.append("Na")
                        
                    #Atribuir plasmídeo
                    try:
                        plasmids.append(f.qualifiers["plasmid"][0])
                    except:
                        plasmids.append("Na")
                        
                    #Atribuir cepa
                    try:
                        strain.append(f.qualifiers["strain"][0])
                    except:
                        strain.append("Na")
    
    df = pd.DataFrame({"accession":accession,
                       "colection_date":colecdate,
                       "host":host,
                       "source":source,
                       "coord":coord,
                       "country":country,
                       "organism":organism,
                       "strain":strain,
                       "plasmid":plasmids
                      })
    df.drop_duplicates('accession', inplace = True)
    df.to_csv('catalogo_teste.csv')
    return df.head()
def filldf(target,source):
    
    """Preenche os dataframes criados pelo abricate com informações adicionais dos arquivos genbank
    ->Pam target: dataframe a ser preechido.
    ->Pam source: dataframe com as informações adicionais"""
    target.set_index('file', inplace = True)

    target['colection_date'] = None #ok
    target['host'] = None
    target["source"] = None #ok
    target['coordenates'] = None #ok
    target['country'] = None #ok
    target['organism'] = None #ok
    target['strain'] = None #ok
    target['plasmid'] = None #ok

    for i in source.itertuples():
        if i.accession in target.index:
            target.loc[i.accession,['host']] = i.host
            target.loc[i.accession,['source']] = i.source
            target.loc[i.accession,['coordenates']] = i.coord
            target.loc[i.accession,['country']] = i.country
            target.loc[i.accession,['organism']] = i.organism
            target.loc[i.accession,['strain']] = i.strain
            target.loc[i.accession,['plasmid']] = i.plasmid
            target.loc[i.accession,['colection_date']] = i.colection_date
    return target