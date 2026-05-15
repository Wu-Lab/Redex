import torch.nn as nn
import torch.nn.functional as F
from dgl.nn import GraphConv 
from dgl.nn import GATConv 
import torch
from dgl.nn import SAGEConv

class LogReg(nn.Module):
    def __init__(self, hid_dim, out_dim):
        super(LogReg, self).__init__()
        self.fc = nn.Linear(hid_dim, out_dim)

    def forward(self, x):
        ret = self.fc(x)
        return ret


class MLP(nn.Module):
    def __init__(self, nfeat, nhid, nclass, use_bn=True):
        super(MLP, self).__init__()

        self.layer1 = nn.Linear(nfeat, nhid, bias=True)
        self.layer2 = nn.Linear(nhid, nclass, bias=True)

        self.bn = nn.BatchNorm1d(nhid)
        self.use_bn = use_bn
        self.act_fn = nn.ReLU()

    def forward(self, _, x):
        x = self.layer1(x)
        if self.use_bn:
            x = self.bn(x)

        x = self.act_fn(x)
        x = self.layer2(x)

        return x


class GraphSAGE(nn.Module):
    def __init__(self, in_dim, hid_dim, out_dim, n_layers, aggregator_type='mean'):
        super(GraphSAGE, self).__init__()

        self.n_layers = n_layers
        self.convs = nn.ModuleList()

        # 第一层 GraphSAGE
        self.convs.append(SAGEConv(in_dim, hid_dim, aggregator_type=aggregator_type))

        # 中间层 GraphSAGE
        if n_layers > 1:
            for i in range(n_layers - 2):
                self.convs.append(SAGEConv(hid_dim, hid_dim, aggregator_type=aggregator_type))
            
            # 最后一层 GraphSAGE
            self.convs.append(SAGEConv(hid_dim, out_dim, aggregator_type=aggregator_type))

    def forward(self, graph, x):
        # 第一层 GraphSAGE
        x = F.relu(self.convs[0](graph, x))

        # 中间层 GraphSAGE
        for i in range(1, self.n_layers - 1):
            x = F.relu(self.convs[i](graph, x))

        # 最后一层 GraphSAGE
        x = self.convs[-1](graph, x)

        return x

    
class GAT(nn.Module):
    def __init__(self, in_dim, hid_dim, out_dim, n_layers, num_heads):
        super().__init__()

        self.n_layers = n_layers
        self.convs = nn.ModuleList()

        # 添加第一层GATConv
        self.convs.append(GATConv(in_dim, hid_dim, num_heads=num_heads, attn_drop=0.01, activation=F.elu))

        # 如果有多于一层的情况
        if n_layers > 1:
            for i in range(n_layers - 2):
                self.convs.append(GATConv(hid_dim * num_heads, hid_dim, num_heads=num_heads, attn_drop=0.01, activation=F.elu))
            self.convs.append(GATConv(hid_dim * num_heads, out_dim, num_heads=num_heads, attn_drop=0.01, activation=None))

    def forward(self, graph, x):
        for i in range(self.n_layers - 1):
            x = self.convs[i](graph, x).flatten(1)
        x = self.convs[-1](graph, x).mean(1)  # 输出层通常使用mean操作以结合多个头

        return x


    
class GCN(nn.Module):
    def __init__(self, in_dim, hid_dim, out_dim, n_layers):
        super().__init__()

        self.n_layers = n_layers
        self.convs = nn.ModuleList()

        self.convs.append(GraphConv(in_dim, hid_dim, norm='both'))

        if n_layers > 1:
            for i in range(n_layers - 2):
                self.convs.append(GraphConv(hid_dim, hid_dim, norm='both'))
            self.convs.append(GraphConv(hid_dim, out_dim, norm='both'))

    def forward(self, graph, x):

        for i in range(self.n_layers - 1):
            x = F.relu(self.convs[i](graph, x))
        x = self.convs[-1](graph, x)

        return x
    
class Decoder(nn.Module):
    def __init__(self, in_dim, out_dim):
        super().__init__()
        self.fc = nn.Linear(in_dim, out_dim)

    def forward(self, z):
        return self.fc(z)
    
# class GCNDecoder(nn.Module):
#     def __init__(self, in_dim, out_dim):
#         super(GCNDecoder, self).__init__()
#         # 使用 DGL 的 GraphConv 层
#         self.gcn = GraphConv(in_dim, out_dim)

#     def forward(self, g, x):
#         # g: DGL 图对象
#         # x: 节点特征 (num_nodes, in_dim)

#         # 应用图卷积层
#         x = self.gcn(g, x)
        
#         return x  
    
class CCA_SSG(nn.Module):
   # def __init__(self, in_dim, hid_dim, out_dim, n_layers,aggregator_type='mean'):
    def __init__(self, in_dim, hid_dim, out_dim, n_layers, num_heads):
        super().__init__()
        
        #self.backbone = GraphSAGE(in_dim, hid_dim, out_dim, n_layers,aggregator_type='mean')
        #self.backbone = GCN(in_dim, hid_dim, out_dim, n_layers)
        self.backbone = GAT(in_dim, hid_dim, out_dim, n_layers, num_heads)
        

        self.decoder = Decoder(out_dim, in_dim)

    def get_embedding(self, graph, feat):
        out = self.backbone(graph, feat)
        return out.detach()
    
    def get_reconstructed_features(self, graph, feat):
        out1 = self.backbone(graph, feat)
        out=self.decoder(out1)
        return out.detach()
    
    def forward(self, graph1, feat1, graph2, feat2):
        h1 = self.backbone(graph1, feat1)
        h2 = self.backbone(graph2, feat2)

        z1 = (h1 - h1.mean(0)) / h1.std(0)
        z2 = (h2 - h2.mean(0)) / h2.std(0)

        decoded1 = self.decoder(z1)
        decoded2 = self.decoder(z2)

        return z1, z2, decoded1, decoded2
        #return z1, z2

