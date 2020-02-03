# -*- coding: utf-8 -*-

from nltk.tokenize import RegexpTokenizer
from nltk.corpus import stopwords
import pandas as pd
import numpy as np
from gensim import corpora
from gensim.matutils import corpus2csc
import scipy.sparse as sp
import json

data = pd.read_csv('export_articles_EGC_2004_2018.csv',sep = "\t",encoding = "utf-8")
data['txt'] = data['title'].astype(str) + ". " +data['abstract'].astype(str)
doc_set = list(data['txt'])
ndoc = len(doc_set)

years = sorted(list(set(data['year'])))

print(years)


#requirements
en_stop = set(stopwords.words('french'))
en_stop.add('les')
en_stop.add('a')
en_stop.add('ce')
en_stop.add('cet')
en_stop.add('cette')
en_stop.add('article')
en_stop.add('approche')
en_stop.add('données')

tokenizer = RegexpTokenizer(r'\w+')



texts = []

# loop through document list
for i in doc_set:
    
    # clean and tokenize document string
    raw = i.lower()
    tokens = tokenizer.tokenize(raw)

    # remove stop words from tokens
    stopped_tokens = [i for i in tokens if not i in en_stop]
    
    # add tokens to list
    texts.append(stopped_tokens)

# turn our tokenized documents into a id <-> term dictionary
dictionary = corpora.Dictionary(texts)
dictionary.filter_extremes(no_below=10)
dictionary[0]
X = {}
yearmap = {}
k = 0
for year in years:
    corpus_k = [dictionary.doc2bow(text) for text in [texts[ind] for ind in [ind[0] for ind in np.argwhere(data['year'] == year)]]]
    term_doc_mat = corpus2csc(corpus_k,num_terms=len(dictionary))
    X[k] = np.dot(term_doc_mat, term_doc_mat.T)
    yearmap[k] = year
    print(X[k].shape)
    k+=1

with open('vocabmap.txt', 'w',encoding = "utf-8") as outfile:
    json.dump(dictionary.id2token, outfile, indent=2)

with open('yearmap.txt', 'w',encoding = "utf-8") as outfile:
    json.dump(yearmap, outfile, indent=2)
    
for i in range(k):
    sp.save_npz('X_'+str(i)+'.npz',X[i])
    
    
    

