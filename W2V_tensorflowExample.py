from __future__ import absolute_import, division, print_function, unicode_literals
from tensorflow.keras import layers,Model
from tensorflow.keras.initializers import Constant
from tensorflow.keras.initializers import RandomUniform
from sklearn.preprocessing import normalize
import scipy.sparse as sp
import numpy as np
import random
import tensorflow as tf
import logging
logger = logging.getLogger()

#X  : sp matrice
#k  : number of negatiev examples
def neg_example(X,k):

    f_t = normalize(np.sum(X,axis = 0), norm='l1')
    freqNeg = np.squeeze(normalize(np.power(f_t,3/4), norm='l1'))

    NN = np.sum(X,axis = 0).T

    Neg = sp.dok_matrix((X.shape[0],X.shape[1]),dtype = np.int)
    for i in range(X.shape[0]):
        Ka = int(NN[i]) * k
        Neg[i,:] = np.random.multinomial(Ka, freqNeg, size=1)
    Neg = Neg.tocoo()
    Neg.eliminate_zeros()
    return Neg

def cooc_to_list(cooc,val):    
    cooc = cooc.tocoo()
    cooc.eliminate_zeros()
    ind = list(zip(cooc.row.tolist(),cooc.col.tolist(),cooc.data.tolist()))
    X = []
    for dat in ind:
        X += [(dat[0],dat[1],val) for i in range(int(dat[2]))]
    return X        
    
class DWE(Model):
    def __init__(self,v,d,T):
        super(DWE, self).__init__()
        self.U = layers.Embedding(v,d)
        #self.V = layers.Embedding(v,d)
        
    def call(self, pair):
        i,j,de = tf.split(pair, [1,1,1], 1)
        de = tf.cast(de, tf.float32)
        u = tf.squeeze(self.U(i))
        #v = tf.squeeze(self.V(j))
        v = tf.squeeze(self.U(j))
        x = tf.reduce_sum(tf.multiply(u, v),axis = 1)
        x = tf.expand_dims(x,1)
        x = tf.multiply(de,x)
        x = tf.math.sigmoid(x)
        return -x
    
@tf.function
def train(model, dataset, optimizer,train_loss):
    for triple in dataset:
        with tf.GradientTape() as tape:
        # training=True is only needed if there are layers with different
        # behavior during training versus inference (e.g. Dropout).
            loss = model(triple, training=True)
        gradients = tape.gradient(loss, model.trainable_variables)
        optimizer.apply_gradients(zip(gradients, model.trainable_variables))
        train_loss(loss)
        
def compute_nn(i,U,V,voc,n_nn):
    
    U = normalize(U)
    V = normalize(V)
    print("NN for : ", end="")
    print(voc[i])
    sim = (U[i,:] @ V.transpose()).flatten()
    
    indi = np.argsort(sim * -1)
    
    d = sim[indi][1:(n_nn+1)]
    g = np.asarray(voc)[indi][1:(n_nn+1)]
    out = dict(zip(g,d))
    print(out)
    

            
import pickle
    
inpu = pickle.load( open( "Data/egc.dwe", "rb" ) )
voc = inpu['voc']
voc2id = dict(zip(voc,list(range(len(voc)))))
cooc = inpu['data']
years = inpu['years']

cooc = sum([*cooc.values()])
#cooc = sum([cooc[12],cooc[13],cooc[14]])
#cooc = cooc[14]

print("Vocabulary of size %d with %d observed cooccurences" % (len(voc),np.sum(cooc)), flush=True)

X_pos = cooc_to_list(cooc,1)
cooc_neg = neg_example(cooc,5)
X_neg = cooc_to_list(cooc_neg,-1)


cooc_f = cooc - cooc_neg
cooc_f = normalize(cooc_f)

compute_nn(voc2id["classification"],cooc_f.todense(),cooc_f.todense(),voc,5)

X = X_pos + X_neg
data = tf.data.Dataset.from_tensor_slices(np.asarray(X)).shuffle(35000000).batch(256)

dwe = DWE(len(voc),200,1)
  
train_loss = tf.keras.metrics.Sum(name='train_loss')
optimizer = tf.keras.optimizers.Adagrad(learning_rate=0.01)

nepochs = 10

print('Starting Learning', flush=True)
ll = []
for epoch in range(nepochs):

    # Reset the metrics at the start of the next epoch
    train_loss.reset_states()
    train(dwe,data,optimizer,train_loss)
    template = 'Epoch {}, Loss: {}'
    print(template.format(epoch+1,train_loss.result()), flush=True)
    ll.append(train_loss.result())

U = np.array(dwe.U(tf.constant(list(range(len(voc))))))
#V = np.array(dwe.V(tf.constant(list(range(len(voc))))))
compute_nn(voc2id["classification"],U,U,voc,5)

import matplotlib.pyplot as plt
plt.plot(ll)
plt.ylabel('loss')
plt.show()

