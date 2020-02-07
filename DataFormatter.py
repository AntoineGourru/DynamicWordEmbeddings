from nltk.corpus import stopwords
from nltk.tokenize import RegexpTokenizer
import pandas as pd
import numpy as np
from sklearn.preprocessing import normalize
from sklearn.feature_extraction.text import TfidfTransformer, TfidfVectorizer
import scipy.sparse as sp
import pickle

#text : list de string
#years: list de years
def create_dataset(doc_set,years,en_stop,l = 5,max_df=0.25, min_df=10):

    vectorizerTF = TfidfVectorizer(lowercase=True, analyzer="word", stop_words=en_stop, max_df=max_df, min_df=min_df, norm=None, use_idf=False)
    tf = vectorizerTF.fit_transform(doc_set)
    ndocs = tf.shape[0]
    vocabulary = list(vectorizerTF.get_feature_names())
    nvoc = len(vocabulary)
    word2id = dict(zip(vocabulary,list(range(nvoc))))
    
    tokenizer = RegexpTokenizer(r'\w+')

    texts = []
    # loop through document list
    for i in doc_set:
        
        # clean and tokenize document string
        raw = i.lower()
        tokens = tokenizer.tokenize(raw)

        # remove stop words from tokens
        stopped_tokens = [i for i in tokens if i in vocabulary]
        
        # add tokens to list
        texts.append(stopped_tokens)

        
    temp = list(set(years))
    temp.sort()

    nyears = len(temp)
    idi = list(range(nyears))
    year2id = dict(zip(temp,idi))

    doc_year = {}
    
    for y in idi:
        doc_year[y] = []
        
    for i in range(ndocs):
        j = year2id.get(years[i])
        doc_year[y].append(i)
        
    print("Extracted Corpus with %d documents and %d words, spanning %d years" % (ndocs,nvoc,nyears) )
    
    X = {}
    print("Creating word pairs", end ="")
    for t in idi:
        print(".", end ="")
        list_doc = doc_year[t]
        
        X[t] = sp.dok_matrix((nvoc,nvoc),dtype = np.int)
        
        for i in list_doc:
            doc = texts[i]
            l_tot = len(doc)
            for i in range(len(doc)):
                w_i = doc[i]
                for elle in range(i-l, i+l):
                    if elle > -1 and elle < l_tot:
                        w_j = doc[elle]
                        X[t][word2id[w_i],word2id[w_j]] += 1
                            
        
        X[t] = X[t].tocoo()

        
    print()
    data = {}
    data['voc'] = vocabulary
    data['data'] = X
    data['years'] = temp

    return data
