import torch
import pandas as pd

"""
This file contains the program for the neural network model training to compute the drug effect of Truvada. 
Aim is to replace the interpolation function which runs super slow...
Program in this file will not be used within this package, and it will be run separately only if new model is needed.
"""

# c: input, concentration matrix, 10000 * 2
# eta: output, drug effect 10000
conc_effect_matrix = pd.read_csv('modMMOA_FTC_TDF_zero_extended.csv')
c_ftc = torch.tensor(conc_effect_matrix['FTC'], dtype=torch.float32)
c_tfv = torch.tensor(conc_effect_matrix['TFV'], dtype=torch.float32)
c = torch.log10(torch.stack((c_ftc, c_tfv), dim=1) + 1e-20)
eta = torch.tensor(conc_effect_matrix['eps'], dtype=torch.float32).unsqueeze(dim=1)
# concatenate concentration and drug effect together so that the matrix can be shuffled later.
data = torch.cat((c, eta), dim=1)

# current nn model. If new model is trained, the code in pd.py must be changed respectively.
model = torch.nn.Sequential(
    torch.nn.Linear(2, 500),
    torch.nn.Sigmoid(),
    torch.nn.Linear(500, 100), 
    torch.nn.Sigmoid(),
    torch.nn.Linear(100, 50),
    torch.nn.Sigmoid(),
    torch.nn.Linear(50, 1)
)

# use Mean Squared Error (MSE) as our loss function.
loss_fn = torch.nn.MSELoss()
model.train()
# at beginning 1e-4, for further training the rate can be reduced to 1e-6 or 5e-7.
learning_rate = 1e-4
# Use the optim package to define an Optimizer that will update the weights of
# the model. Here weuse RMSprop;
optimizer = torch.optim.RMSprop(model.parameters(), lr=learning_rate)
for t in range(2000):
    index = torch.randperm(10000)
    data_shuffled = data[index]
    n = 1
    for i in range(n):
        c = data_shuffled[i * 10000 // n:(i + 1) * 10000 // n, :2]
        eta = data_shuffled[i * 10000 // n:(i + 1) * 10000 // n, [2]]

        # Forward pass: compute predicted eta by passing c to the model.
        eta_pred = model(c)

        # Compute loss.
        loss = loss_fn(eta, eta_pred)

        # Before the backward pass, use the optimizer object to zero all of the
        # gradients for the variables it will update (which are the learnable
        # weights of the model). This is because by default, gradients are
        # accumulated in buffers( i.e, not overwritten) whenever .backward()
        # is called.
        optimizer.zero_grad()

        # Backward pass: compute gradient of the loss with respect to model parameters
        loss.backward()

        # Calling the step function on an Optimizer makes an update to its parameters
        optimizer.step()
