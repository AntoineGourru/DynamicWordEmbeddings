import pandas as pd
import numpy as np
import pickle

from DataFormatter import create_dataset
from nltk.corpus import stopwords

# EGC
en_stop = set(stopwords.words('french'))
data = pd.read_csv('Data/egc.csv',sep = "\t")

data['txt'] = data['title'].astype(str) + ". " +data['abstract'].astype(str)
doc_set = list(data['txt'])
years = np.array(data['year'])
years = years.flatten().tolist()

dataset = create_dataset(doc_set,years,en_stop)
pickle.dump(dataset, open( "Data/egc.dwe", "wb" ) )

inpu = pickle.load( open( "Data/egc.dwe", "rb" ) )
voc = inpu['voc']
X = inpu['data']
years = inpu['years']

print(years)
