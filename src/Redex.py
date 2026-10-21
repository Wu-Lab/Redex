import torch
from preprocess import preprocess_adj, preprocess, construct_interaction, construct_interaction_KNN, get_feature, fix_seed
import dgl
import numpy as np
import numpy as np
from tqdm import tqdm
from torch import nn
import torch.nn.functional as F
import pandas as pd
from CCAmodel import CCA_SSG, LogReg
from aug import random_aug
    
class Redex():
    def __init__(self, 
        adata,
        dataset_name,
        gene_path,
        device= torch.device('cpu')
        ):
        '''\

     
      
        '''

        self.epochs=600
        self.dim_input=3000
        self.gene_path=gene_path
        self.alpha=9
        self.der=0.00001 
        self.beta=1
        self.lambd=0.01
        self.dfr=0.01
        self.num_heads=6
        self.learning_rate=0.001
        self.weight_decay=0.0001
        self.dataset_name=dataset_name
        self.dim_hidden=128
        self.n_layers=2
        self.dim_output123=64

         

        print( 'self.der',self.der)
        print( 'self.alpha',self.alpha)
        print( 'self.beta',self.beta)
        
        print('   self.num_heads',   self.num_heads)

       
        self.random_seed = 2024
        self.adata = adata.copy()
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.gene_path=gene_path
        self.dim_input=3000
      
        fix_seed(self.random_seed)
        
        if 'highly_variable_gene_yuan' not in adata.var.keys():
           preprocess(self.adata,dataset_name=self.dataset_name,gene_path=self.gene_path)
        
        if 'adj' not in adata.obsm.keys():
            #construct_interaction_KNN(self.adata)
            construct_interaction(self.adata)
           
        if 'feat' not in adata.obsm.keys():
           get_feature(self.adata)
       

        def adata_to_dgl(adata):
            adj = adata.obsm['adj']
            src, dst = np.nonzero(adj)
            graph = dgl.graph((src, dst))
            if 'feat' in adata.obsm:
                graph.ndata['feat'] = torch.tensor(adata.obsm['feat'], dtype=torch.float32)
            if 'label' in adata.obs:
                graph.ndata['label'] = torch.tensor(adata.obs['label'].values, dtype=torch.int64)
            return graph

        # 假设 adata 是你的 AnnData 对象
        self.graph = adata_to_dgl(self.adata )
        self.features = torch.FloatTensor(self.adata.obsm['feat'].copy()).to(self.device)
        self.adj = self.adata.obsm['adj']
        self.n_spot = adata.n_obs
        self.num_edges = np.sum(np.triu(self.adj) > 0)
        self.dim_input = self.features.shape[1]

        self.adj = preprocess_adj(self.adj)
        self.adj = torch.FloatTensor(self.adj).to(self.device)
        
        
            
    def train(self): 
        self.model = CCA_SSG(self.dim_input,self.dim_hidden,self.dim_output123,self.n_layers,self.num_heads).to(self.device)
        self.optimizer = torch.optim.Adam(self.model.parameters(), self.learning_rate, 
                                          weight_decay=self.weight_decay)
        
        print('Begin to train ST data...')
        self.model.train()
        
        for epoch in tqdm(range(self.epochs)): 
            self.model.train()
            self.optimizer.zero_grad()
    
            graph1, feat1 = random_aug(self.graph, self.features,self.dfr, self.der)
            graph2, feat2 = random_aug(self.graph, self.features,self.dfr, self.der )

            graph1 = graph1.to(self.device)
            graph2 = graph2.to(self.device)

            feat1 = feat1.to(self.device)
            feat2 = feat2.to(self.device)

            z1, z2, decode1,decode2= self.model(graph1, feat1, graph2, feat2)

            c = torch.mm(z1.T, z2)
            c1 = torch.mm(z1.T, z1)
            c2 = torch.mm(z2.T, z2)

            c = c / self.n_spot
            c1 = c1 / self.n_spot
            c2 = c2 / self.n_spot

            self.loss_inv = -torch.diagonal(c).sum()
            iden = torch.tensor(np.eye(c.shape[0])).to(self.device)
            self.loss_dec1 = (iden - c1).pow(2).sum()
            self.loss_dec2 = (iden - c2).pow(2).sum()


            self.loss_reconstruction1= F.mse_loss(feat1, decode1)
            self.loss_reconstruction2= F.mse_loss(feat2, decode2)

            self.loss =   self.beta*self.loss_inv + self.lambd * (self.loss_dec1 + self.loss_dec2)+ self.alpha*(self.loss_reconstruction1+self.loss_reconstruction2)
            
            self.loss.backward() 
            self.optimizer.step()
        print("Optimization finished for ST data!")
        
        with torch.no_grad():
            self.model.eval()
            graph = self.graph.to(self.device)
            feat = self.features.to(self.device)
            self.hiden_emb= self.model.get_embedding(graph, feat).detach().cpu().numpy()
            self.reconsfeat= self.model.get_reconstructed_features(graph, feat).detach().cpu().numpy()
        self.adata.obsm['hideemb'] = self.hiden_emb
        self.adata.obsm['reconsfeat'] = self.reconsfeat
                    
        return self.adata
     
