from __future__ import absolute_import, division, print_function, unicode_literals
from tensorflow.keras import layers,Model
from tensorflow.keras.initializers import Constant
from sklearn.preprocessing import normalize
import scipy.sparse as sp
import numpy as np
import random
import tensorflow as tf
tf.enable_eager_execution()

class DWE(Model):
    def __init__(self,v,d,T):
        super(DWE, self).__init__()
        self.embeddings = layers.Embedding(v,d)
        self.d = d
        
    def call(self, pair):
        i,j,de = tf.split(pair, [1,1,1], 1)
        de = tf.cast(de, tf.float32)
        u = tf.squeeze(self.embeddings(i))
        v = tf.squeeze(self.embeddings(j))
        x = tf.reduce_sum(tf.multiply(u, v),axis = 1)
        x = tf.expand_dims(x,1)
        x = tf.multiply(de,x)
        x = tf.math.sigmoid(x)
        return -x
    
    @tf.function
    def train_step(self,triple,optimizer,train_loss):     
        with tf.GradientTape() as tape:
        # training=True is only needed if there are layers with different
        # behavior during training versus inference (e.g. Dropout).
            loss = self(triple, training=True)
        gradients = tape.gradient(loss, self.trainable_variables)
        optimizer.apply_gradients(zip(gradients, self.trainable_variables))

        train_loss(loss)


    def fit(self,train_ds,nepochs):

        train_loss = tf.keras.metrics.Sum(name='train_loss')
        optimizer = tf.keras.optimizers.Adam()
        
        for epoch in range(nepochs):
            # Reset the metrics at the start of the next epoch
            train_loss.reset_states()

            for triple in train_ds:
                self.train_step(triple,optimizer,train_loss)

            template = 'Epoch {}, Loss: {}'
            print(template.format(epoch+1,train_loss.result()))

def compute_nn(i,embeddings,voc,n_nn):
    
    embeddings = normalize(embeddings, axis=1)
    print("NN for : ", end="")
    print(voc[i])
    dist = (embeddings[i,:] @ embeddings.transpose()).flatten()
    
    indi = np.argsort(dist * -1)[1:(n_nn + 1)]
    
    d = dist[indi]
    g = []
    for i in indi:
        g.append(voc[i])
    out = dict(zip(g,d))    
    print(out)

            
import pickle
    
inpu = pickle.load( open( "Data/egc.dwe", "rb" ) )
voc = inpu['voc']
print(voc)
voc2id = dict(zip(voc,list(range(len(voc)))))
X = inpu['data']
years = inpu['years']

dwe = DWE(len(voc),10,1)
data = tf.data.Dataset.from_tensor_slices(np.asarray(X[0])).batch(512)

dwe.fit(data,100)

U = np.array(dwe.embeddings(tf.constant(list(range(len(voc))))))

compute_nn(voc2id["clusters"],U,voc,3)
