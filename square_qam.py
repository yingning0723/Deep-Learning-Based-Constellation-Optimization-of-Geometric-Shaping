# bitAE for QAM 只训练编码器，用于形成square qam
import torch
import matplotlib.pyplot as plt
from torch import nn
import numpy as np
import scipy.io as sio
from torch.optim import Adam,lr_scheduler
import torch.utils.data as Data
import os
from scipy.io import loadmat

os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
# 调制阶数
symbol_dim = 2
bits_per_symbol = 10
modulation_order = 2 ** bits_per_symbol
pattern = "square-qam1024"
# 模型参数
device = torch.device("cpu:0")
BATCH_SIZE = 2000
NUM_EPOCHS = 80000
train_num = BATCH_SIZE
test_num = 20000
learning_rate = 1e-3
USE_CUDA = True

class bit_AE(nn.Module):
    def __init__(self, nn_points):
        super(bit_AE, self).__init__()
        self.encoder4 = nn.Sequential(
            nn.Linear(bits_per_symbol, nn_points),
            nn.ReLU(),
            nn.Linear(nn_points, nn_points),
            nn.ReLU(),
            nn.Linear(nn_points, nn_points),
            nn.ReLU(),
            nn.Linear(nn_points, nn_points),
            nn.ReLU(),
            nn.Linear(nn_points, nn_points),
            nn.ReLU(),
            nn.Linear(nn_points, 2),
        )

    def encode_signal(self, x):
        return self.encoder(x)

    def forward(self,b,snr = 7):
        b = b.to(device)
        x = self.encoder4(b)
        batchsize = len(b)
        norm = torch.sum(x ** 2) / batchsize
        x_norm = x / torch.sqrt(norm)

        return x_norm

def train(path , snr ):
    model = bit_AE(4096)
    for layer in model.modules():
        if isinstance(layer, nn.Linear):
            nn.init.xavier_normal_(layer.weight)  # 使用Xavier均匀分布初始化权重
            nn.init.constant_(layer.bias, 0.002)  # 使用常量值初始化偏差
    model.to(device)

    train_labels = torch.randint(low = 0,high=2,size=[train_num,bits_per_symbol],dtype=torch.float)
    train_data = train_labels
    train_labels = torch.zeros([BATCH_SIZE,2],dtype=torch.float)
    mat_data = loadmat(path + 'square_qam' +str(modulation_order) +'.mat')['m1']
    for i in range(BATCH_SIZE):
        index = int(''.join(map(str, train_data[i].int().tolist())), 2)
        train_labels[i, 0] = mat_data[0, index]
        train_labels[i, 1] = mat_data[1, index]

# DataBase in Pytorch
    dataset = Data.TensorDataset( train_data,  train_labels)
    train_loader = Data.DataLoader(dataset = dataset, batch_size = BATCH_SIZE, shuffle = True, num_workers = 0)
#optmizer & Loss
    optimizer = Adam(model.parameters(),lr=learning_rate)
    loss_fn = nn.MSELoss()
#Training
    for epoch in range(NUM_EPOCHS):
       for step, (x, y) in enumerate(train_loader):
            y = y.float()
            x = x.to(device)
            y = y.to(device)
            mapped = model(x,snr)
            loss = loss_fn(mapped, y)
            optimizer.zero_grad()               # clear gradients for this training step
            loss.backward()                     # backpropagation, compute gradients
            optimizer.step()                    # apply gradients

            if epoch%10 ==0:
                print('Epoch: ', epoch, 'loss: ', loss.item())
                torch.save(model, path + pattern + 'SNR{}.pth'.format(snr))
                torch.save(model.state_dict(), path +'square_qam'+str(modulation_order)+'_{}dB_goodbase.pt'.format(snr))
                # 画星座图
                index = torch.linspace(0, modulation_order - 1, modulation_order).long()
                bits = [bin(i)[2:].zfill(bits_per_symbol) for i in index]
                test_data = np.zeros([len(bits), bits_per_symbol])
                for i in range(len(bits)):
                    for j in range(bits_per_symbol):
                        test_data[i, j] = int(bits[i][j])
                test_data = torch.tensor(test_data, dtype=torch.float)
                x = model(test_data.to(device), snr)
                plot_data = x.to("cpu").data.numpy()
                plot_Constellation(path+"plot",plot_data,snr, epoch,modulation_order)

    return ber

def plot_Constellation(path ,plot_data,snr, epoch,modulation_order):
    plt.figure()
    plt.scatter(plot_data[:, 0], plot_data[:, 1])
    plt.axis((-2, 2, -2, 2))
    plt.grid()
    if epoch > -1:
        plt.savefig(path + "qam{}epoch{}_snr{}dB_bit_wise.png".format(modulation_order,epoch,snr))
    else:
        index = torch.linspace(0, modulation_order - 1, modulation_order).long()
        label_text = [bin(i)[2:].zfill(bits_per_symbol) for i in index]
        plt.axis((-1.5, 1.5, -1.5, 1.5))
        for i in range(len(plot_data)):
            plt.text(plot_data[i, 0], plot_data[i, 1], label_text[i])
        plt.show()
    plt.close()
def frange(x, y, jump):
    while x < y:
        yield x
        x += jump


if __name__ == "__main__":
    ber = []
    ser = []
    path = r"D:\\graduation\\qam1024\\"
    train(path, 0)
