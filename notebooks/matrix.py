def heatmatrix (df, mode ):
    import pandas as pd
    if mode ==  'class':
        df['class'] = df.gene.str.split(')').apply(lambda x: x[0][1:])
        dfgenomes = df[['file','class','species']]
        
        files = []
        index = 0
        for sp in dfgenomes['species']:
            newfile = f'{sp}_{dfgenomes.file[index]}'
            newfile = newfile.split('.')[0]
            files.append(newfile)
            index = index + 1
        dfgenomes['species'] = files
        
        dfindex = pd.DataFrame(dfgenomes['species'])
        dfgenomes = pd.get_dummies(dfgenomes['class'])
        matrixgenomes = pd.concat([dfindex,dfgenomes], axis = 1, sort = False)
        matrixgenomes = matrixgenomes.groupby(['species']).sum()
        
    elif mode == 'gene':
        dfgenomes = df[['file','gene','species']]
        
        files = []
        index = 0
        for sp in dfgenomes['species']:
            newfile = f'{sp}_{dfgenomes.file[index]}'
            newfile = newfile.split('.')[0]
            files.append(newfile)
            index = index + 1
        dfgenomes['species'] = files
        
        dfindex = pd.DataFrame(dfgenomes['species'])
        dfgenomes = pd.get_dummies(dfgenomes['gene'])
        matrixgenomes = pd.concat([dfindex,dfgenomes], axis = 1, sort = False)
        matrixgenomes = matrixgenomes.groupby(['species']).sum()
    return matrixgenomes

def matrix (df,x_axis, y_axix):
    import pandas as pd
    redes = df[[y_axix,x_axis]]
    dummy = pd.get_dummies(df[y_axix])
    files = pd.DataFrame(df[x_axis])
    dummy = pd.concat([files,dummy], axis=1 )
    matrix = dummy.groupby([x_axis]).sum()
    return matrix

def matrix_box(data,specie, gen = False):
    import pandas as pd
    df = data[:]
    if gen == False:
        df['class'] = df.gene.str.split(')').apply(lambda x: x[0][1:])
        df.set_index('species', inplace = True)
        df = df.loc[f'{specie}',['file','class']]
        df.reset_index(inplace =True)
        df.drop('species', axis = 1, inplace = True)
        df.set_index('file', inplace =  True)
        df = pd.get_dummies(df['class'])
        df = df.groupby('file').sum()
    else:
        df['specie'] = df.specie.str.split(' ').apply(lambda x: x[0])
        df['class'] = df.gene.str.split(')').apply(lambda x: x[0][1:])
        df.set_index('species', inplace = True)
        df = df.loc[f'{specie}',['file','class']]
        df.reset_index(inplace =True)
        df.drop('species', axis = 1, inplace = True)
        df.set_index('file', inplace =  True)
        df = pd.get_dummies(df['class'])
        df = df.groupby('file').sum()
    return df