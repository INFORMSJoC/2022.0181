## GAT
import os
import numpy as np
import pandas as pd
import torch
from torch.nn import Linear, MSELoss
import torch.nn.functional as F
from torch_geometric.data import Data
from torch_geometric.nn import global_mean_pool, GATConv
from torch_geometric.nn import GCNConv
from sklearn.model_selection import ParameterGrid
from scipy.spatial.distance import pdist, squareform
from torch.nn import MSELoss
from sklearn.preprocessing import StandardScaler
import time

 
class GATNet(torch.nn.Module):
    def __init__(self, hidden_dim):
        super(GATNet, self).__init__()
        self.conv1 = GATConv(5, hidden_dim)
        self.conv2 = GATConv(hidden_dim, 1) 

    def forward(self, data):
        x, edge_index = data.x, data.edge_index
        x = self.conv1(x, edge_index)
        x = F.relu(x)
        x = self.conv2(x, edge_index) 
        return x
    
def construct_graph(X, threshold):
    distances = pdist(np.transpose(X), 'euclidean')
    adjacency = squareform(distances)
    adjacency[adjacency > threshold] = 0
    adjacency = torch.tensor(adjacency, dtype=torch.long)
    edges = adjacency.nonzero(as_tuple=False).t().contiguous()
    return edges

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
# Parameters for the GNN model
#%%
path_result = "Results-Sim23_7_25_Taobao_TS"
# Grid of hidden dimensions to search
grid = {'hidden_dim': [2, 4, 8, 16, 32,64]}
#grid = {'hidden_dim': [32]}
M = 90 # specify the number of sets here 
results = [] 
epochs = 100
loss_func = MSELoss()
scaler = StandardScaler()
threshold = 50

start_time = time.time() 
for m in range(1, M+1):
    print(m)  
    # Read data from CSV files
    X = pd.read_csv(f"{path_result}/Xused{m}.csv").values
    Y = pd.read_csv(f'{path_result}/Yused{m}.csv').values 
    
    # Only use the first 120 timestamps to construct the graph
    X_trained = X[:120]  
    Y_trained = Y[:120] 
    # edge = torch.tensor([[i, j] for i in range(150) for j in range(150) if i != j], dtype=torch.long).t().contiguous() # Create edge index for a complete graph
    edges = construct_graph(Y_trained,threshold= threshold)

    # Normalize the training data
    X = scaler.fit_transform(X)
    Y = scaler.fit_transform(Y)

    # Reshape the X data to have 150 nodes, each with 5 features
    reshaped = X.reshape(-1, 150, 5)
    X = torch.tensor(reshaped, dtype=torch.float)
    Y = torch.tensor(Y, dtype=torch.float)

    for params in ParameterGrid(grid):
        # Instantiate the network and optimizer 
        model = GATNet(params['hidden_dim']).to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=0.01, weight_decay=5e-4)

        # Train the network with the data from the first 120 timestamps
        model.train()
        for epoch in range(epochs):
            for t in range(120):
                data = Data(x=X[t].to(device), edge_index=edges.to(device))
                optimizer.zero_grad()
                out = model(data)
                loss = loss_func(Y[t].unsqueeze(1).to(device),out) 
                loss.backward()
                optimizer.step()
     
        # Use the trained GNN to predict the outcomes at t = 121
        model.eval()
        with torch.no_grad():
            data = Data(x=X[120].to(device), edge_index=edges.to(device))
            pred = model(data)
            pred = np.transpose(pred.cpu().detach().numpy())

            # inverse transform the prediction and ground truth to the original scale
            pred_original_scale = scaler.inverse_transform(pred)
            Y_original_scale = scaler.inverse_transform(Y[120].unsqueeze(0).cpu().detach().numpy())

            f_norm = np.linalg.norm(pred_original_scale - Y_original_scale)
            relative_f_norm = f_norm / np.linalg.norm(Y_original_scale)
            # Record the results
            results.append({'dataset': m, 'optimal_hidden_dim': params['hidden_dim'], 'relative_f_norm': relative_f_norm})
        


end_time = time.time()    
elapsed_time = end_time - start_time
total_time_minutes = elapsed_time /60
print(f"Total running time in minutes: {total_time_minutes}") 

results_df = pd.DataFrame(results)
print(results_df)
