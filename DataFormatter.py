from nltk.corpus import stopwords
import pandas as pd
import numpy as np
from sklearn.preprocessing import normalize
from sklearn.feature_extraction.text import TfidfTransformer, TfidfVectorizer
import scipy.sparse as sp
import pickle

#text : list de string
#years: list de years
def create_dataset(doc_set,years,en_stop,max_df=0.25, min_df=10):

    vectorizerTF = TfidfVectorizer(lowercase=True, analyzer="word", stop_words=en_stop, max_df=max_df, min_df=min_df, norm=None, use_idf=False)
    tf = vectorizerTF.fit_transform(doc_set)
    ndocs = tf.shape[0]
    vocabulary = vectorizerTF.get_feature_names()
    nvoc = len(vocabulary)
      

    temp = list(set(years))
    temp.sort()

    nyears = len(temp)
    idi = list(range(nyears))
    year2id = dict(zip(temp,idi))

    doc_year = sp.dok_matrix((ndocs,nyears),dtype = np.int)
    for i in range(ndocs):
        j = year2id.get(years[i])
        doc_year[i,j] = 1
        
    doc_year = doc_year.tocsc()
    print("Extracted Corpus with %d documents and %d words, spanning %d years" % (ndocs,nvoc,nyears) )
    
    X = {}
    print("Creating word pairs", end ="")
    for t in idi:
        print(".", end ="")
        temp_tf = sp.diags(doc_year[:,t].toarray().reshape((ndocs, )),0) @ tf
        
        X[t] = (temp_tf.T @ temp_tf).tocoo()

        
    print()
    data = {}
    data['voc'] = vocabulary
    data['data'] = X
    data['years'] = temp

    return data
