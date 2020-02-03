import copy
import numpy as np
import scipy.sparse as sp
import math
from numpy.linalg import inv
import warnings
warnings.filterwarnings("ignore")

sigmoid_uv = lambda x: 1 / (1 + math.exp(-x))
sigmoid = np.vectorize(sigmoid_uv)
log_sigmoid_uv = lambda x: -math.log(1 + math.exp(-x))
log_sigmoid = np.vectorize(log_sigmoid_uv)
n_log_sigmoid_uv = lambda x: -math.log(1 + math.exp(x))
n_log_sigmoid = np.vectorize(n_log_sigmoid_uv)

class Variationnal_Parameters:

    def __init__(self,D, K,tau):

        self.u_mu = np.random.multivariate_normal(np.zeros(K), np.diag(np.full(K,1/tau)), size = D)
        self.v_mu = np.random.multivariate_normal(np.zeros(K), np.diag(np.full(K,1/tau)), size = D)
        self.u_sigma = np.full((D, K), 1/tau)
        self.v_sigma = np.full((D, K), 1/tau)
        self.xi = np.zeros([D, D]) 
    
class BDWE:

    def __init__(self,D, K,tau,sigma_t,T,nb_epochs):

        self.D = D
        self.K = K
        self.sigma_t = sigma_t
        self.tau = tau
        self.T = T
        self.nb_epochs = nb_epochs

        self.vp = self.draw_vp()

        self.ll = [0 for i in range(T)]
        
    def draw_vp(self):
        
        vp = {}
        draw = Variationnal_Parameters(self.D,self.K,self.tau)
        
        for t in range(self.T):
            vp[t] = draw
            
        return vp 
    
    def draw_ne(X):
        N = X.shape[0]
        
        summe = np.sum(X,axis = 1)
        freqNeg = np.power((summe/np.sum(summe)),(3/4))
        freqNeg = freqNeg /np.sum(freqNeg)
        
        freqNeg = np.array(freqNeg).flatten()
        
        Neg = X.copy()
        
        for i in range(N):
            Neg[i] = np.random.multinomial(summe[i],freqNeg, size=1)

        Neg.eliminate_zeros()    
        return Neg
    
    def likelihood(self,data,t):

        vp = self.vp[t]
        x = data.X[t]
        
        x = vp.u_mu @ vp.v_mu.transpose()
        
        sig = log_sigmoid(x)

        gauche = data.X[t] * sig

        sig = n_log_sigmoid(x)

        droite = BDWE.draw_ne(data.X[t]) * sig

        ll = np.sum(gauche + droite)

        #ll += sum(apply(vp.u_mu,1,function(x){log(dmvnorm(x,rep(0,model.K),diag(rep(1/model.tau,model.K))))}))
        #ll += sum(apply(vp.v_mu,1,function(x){log(dmvnorm(x,rep(0,model.K),diag(rep(1/model.tau,model.K))))}))
        return ll
    
    def fit(self,data):
        for t in range(self.T):
            self.optim(data,t)

          
    
    def build_C(X):
        
        C = X.copy()

        freq = np.sum(X,axis = 1)
        freq = freq /np.sum(freq)

        #A vérifier !!!!
        freq =  1 - np.sqrt(0.0001/freq)
        #freq =  1 - np.sqrt(0.1/freq)
        
        freq[np.isnan(freq)] = 0
        freq[freq<0] = 0
        #print(freq)
        N = X.shape[0]
        nz = X.nonzero()
        nz = np.concatenate((nz[0].reshape((nz[0].shape[0],1)),nz[1].reshape((nz[1].shape[0],1))),axis = 1)
        for sete in nz:
            i = sete[0]
            j = sete[1]
            C[i,j] = X[i,j] - np.random.binomial(X[i,j],freq[j])

        C -= BDWE.draw_ne(X)
        C.eliminate_zeros()
        return C
        


    def optim(self,data,t):

        if t == 0:
            vp = self.vp[0]
        else:
            vp = copy.deepcopy(self.vp[(t-1)])
            temp_pri = copy.deepcopy(vp)
  
            #Temporal prior
            for i in range(self.D):
                temp_pri.u_sigma[i,:] = 1/(vp.u_sigma[i,:] + 1 + self.tau)
                temp_pri.v_sigma[i,:] = 1/(vp.v_sigma[i,:] + 1 + self.tau)
                temp_pri.u_mu[i,:] = (np.diag(temp_pri.u_sigma[i,:]) @ np.diag(vp.u_sigma[i,:] + 1)) @ vp.u_mu[i,:]
                temp_pri.v_mu[i,:] = (np.diag(temp_pri.v_sigma[i,:]) @ np.diag(vp.v_sigma[i,:] + 1)) @ vp.v_mu[i,:]

        #Init P
        oldPu = {}
        oldPv = {}
        oldRu = {}
        oldRv = {}
        beta = 0
                
        X = data.X[t]
        
        for epo in range(self.nb_epochs):
            
            print("eopch :%d" % epo)
            print("LL :%d" % self.likelihood(data,t))
            
            C = BDWE.build_C(data.X[t])
            print('C done')

            nnz = np.zeros(self.D)
            
            for i in range(self.D):
                
                
                nnz[i] = C[i].nnz
            
            for i in range(self.D):
                
                if nnz[i]>0:

                    a_i = vp.u_sigma[i,:] + vp.u_mu[i,:] * vp.u_mu[i,:]
                    b_i = vp.v_sigma[i,:] + vp.v_mu[i,:] * vp.v_mu[i,:]


                    P_u = np.zeros((self.K,self.K))
                    R_u = np.zeros((self.K,))

                    P_v = np.zeros((self.K,self.K))
                    R_v = np.zeros((self.K,))
                    
                    for j in range(self.D):
                        
                        a_j = vp.u_sigma[j,] + vp.u_mu[j,] * vp.u_mu[j,]
                        b_j = vp.v_sigma[j,] + vp.v_mu[j,] * vp.v_mu[j,]

                        xij = np.sqrt(np.dot(a_i,b_j))

                        deux_lambda_xi = (1/xij)  * (sigmoid(xij)  - 0.5)


                        E = np.diag(vp.v_sigma[j,:]) + ( vp.v_mu[j,:].reshape((self.K,1)) @ vp.v_mu[j,:].reshape((1,self.K)) ) 

                        P_u += np.absolute(C[i,j]) * ( deux_lambda_xi *  E )

                        
                        R_u += 0.5 * C[i,j] * vp.v_mu[j,:]

                        
                        xji = np.sqrt(np.dot(b_i,a_j))
                        deux_lambda_xi = (1/xji)  * (sigmoid(xji)  - 0.5)

                        E = np.diag(vp.u_sigma[j,:]) + ( vp.u_mu[j,:].reshape((self.K,1)) @ vp.u_mu[j,:].reshape((1,self.K)) ) 

                        P_v += np.absolute(C[j,i]) * ( deux_lambda_xi *  E )

                        
                        R_v += 0.5 * C[i,j] * vp.u_mu[j,:]

                    
                    if t == 0:

                        P_u += self.tau
                        P_v += self.tau
                        
                    else:

                        P_u += self.tau + np.diag(1/(temp_pri.u_sigma[i,:]))
                        P_v += self.tau + np.diag(1/(temp_pri.v_sigma[i,:]))

                        R_u += (np.diag(1/(temp_pri.u_sigma[i,:])) @ temp_pri.u_mu[i,:])
                        R_v += (np.diag(1/(temp_pri.v_sigma[i,:])) @ temp_pri.v_mu[i,:])

                    if epo > 4 :
                        beta = np.power((epo-4),(-0.7))
                        P_u =  beta * P_u + (1-beta) * oldPu[i]
                        P_v =  beta * P_v + (1-beta) * oldPv[i]
                        R_u =  beta * R_u + (1-beta) * oldRu[i]
                        R_v =  beta * R_v + (1-beta) * oldRv[i]
                        
                    # Update U
                    
                    Pm1 = inv(P_u)
                    vp.u_mu[i,:] = Pm1 @ R_u

                    vp.u_sigma[i,:] = np.diag(Pm1)
                
                    # Update V
                    
                    Pm1 = inv(P_v)
                    vp.v_mu[i,:] = Pm1 @ R_v
                    
                    vp.v_sigma[i,:] = np.diag(Pm1)
                
                    #print(vp.v_mu[i,:])
                    
                    oldPu[i] = P_u
                    oldPv[i] = P_v
                    oldRu[i] = R_u
                    oldRv[i] = R_v


                else:
                    
                    oldPu[i] = np.zeros((self.K,self.K))
                    oldPv[i] = np.zeros((self.K,self.K))
                    oldRu[i] = np.zeros((self.K,))
                    oldRv[i] = np.zeros((self.K,))

            logL = self.likelihood(data,t)

            
            # print(beta)
        
        

        self.vp[t] = vp


class Data:
    def __init__(self, vocab,X, T):
        self.vocab = vocab
        self.X = X
        self.T = T
        self.D = X[0].shape[0]

'''
D = 20
K = 4
T = 2
X = {}
for i in range(T):
    X[i] = sp.csr_matrix(np.random.randint(D * D, size=(D, D)))
    
vocab = ["un"]*10
data = Data(vocab,X, T)
#model = BDWE(D, K,tau,sigma_t,T,nb_epochs)
model = BDWE(D,K,1,1,T,20)
ll = model.likelihood(data,0)
ll = model.likelihood(data,1)

model.fit(data)
'''

import json

with open('vocabmap.txt') as f:
    dicti = json.load(f)
  
X = {}
for i in range(15):
    X[i] = sp.load_npz('X_'+str(i)+'.npz')
T = 15    
data = Data([*dicti.values()],X, T)
model = BDWE(X[0].shape[0],10,1,1,T,8)
model.fit(data)

for i in range(15):
    np.savetxt('U_'+str(i)+'.emb',model.vp[i].u_mu)
