import os
import scipy.io
import mat73
import numpy as np
import matplotlib.pyplot as plt
import torch
import torch.nn as nn
import torch.nn.functional as F
import random
from torch.optim.lr_scheduler import StepLR

# =============================================================================
# 1. Initialization and Reproducibility
# =============================================================================
def set_seed(seed=42):
    """
    Sets the seed for all pseudo random number generators to ensure 
    reproducible model training and initialization.
    """
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

# Apply the seed function immediately after imports
set_seed(42)

# Set computation device
device = torch.device("cpu")
# torch.cuda.set_device(3)

epochs = 2000

# =============================================================================
# 2. Load Training Data
# =============================================================================
# Load the dictionary and labels prepared by the MATLAB simulation
mat = scipy.io.loadmat("training_data_MLP_with_QSM_modeling.mat")
data = torch.FloatTensor(mat.get("data"))
labels = torch.FloatTensor(mat.get("labels"))
N_TE = mat.get("N_TE")

# =============================================================================
# 3. Model Definition
# =============================================================================
class MLP(nn.Module):
    """
    Multi Layer Perceptron for Multi Mechanistic Source Separation (MMSS).
    Takes MRI magnitude and QSM inputs and outputs microstructural parameters.
    """
    def __init__(self, in_features=N_TE[0][0]+1, h1=64, h2=64, h3=64, out_features=5):
        super().__init__() 
        self.fc1 = nn.Linear(in_features, h1)
        self.fc2 = nn.Linear(h1, h2)
        self.fc3 = nn.Linear(h2, h3)
        self.out = nn.Linear(h3, out_features)

    def forward(self, x):
        x = (F.leaky_relu(self.fc1(x), negative_slope=0.1))
        x = (F.leaky_relu(self.fc2(x), negative_slope=0.1))
        x = (F.leaky_relu(self.fc3(x), negative_slope=0.1))
        x = F.leaky_relu(self.out(x), negative_slope=0.1)
        return x
  
model = MLP()

# =============================================================================
# 4. Training Setup
# =============================================================================
def loss_fn(pred, label):
    """Calculates the L1 Loss between predictions and ground truth labels."""
    return F.l1_loss(pred, label) 

optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)
scheduler = StepLR(optimizer, step_size=500, gamma=0.1) 

# Move data and model to the selected device
data = data.to(device)
labels = labels.to(device)
model = model.to(device)

# =============================================================================
# 5. Training Loop
# =============================================================================
for i in range(epochs):
    # Forward pass
    y_pred = model.forward(data)
    
    # Calculate loss
    loss = loss_fn(y_pred, labels)

    # Backpropagation
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()

    # Clean up memory
    del y_pred
    del loss
    
    # Update learning rate
    scheduler.step()  

# Free up memory before inference
del data
del labels

# =============================================================================
# 6. Inference on In Vivo Data
# =============================================================================
# Load the processed in vivo MRI data
data = mat73.loadmat("data.mat")
input = torch.FloatTensor(data.get("input"))

# Set model to evaluation mode
model.eval()
input = input.to(device)

# Perform inference without tracking gradients
with torch.no_grad():  
    result = model.forward(input)

# =============================================================================
# 7. Save Results
# =============================================================================
# Structure the results and save back to a .mat file for MATLAB processing
struct = {"MMSS_results": result.data.cpu().numpy()}
scipy.io.savemat("MMSS_results.mat", struct)

# Final memory cleanup
del input
del result