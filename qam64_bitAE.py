# bitAE for QAM64
from scipy.io import savemat
import torch
import matplotlib.pyplot as plt
from torch import nn
import numpy as np
import scipy.io as sio
from torch.optim import Adam,lr_scheduler
import torch.utils.data as Data
import os

os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
# 调制阶数
symbol_dim = 2
bits_per_symbol = 6
modulation_order = 2 ** bits_per_symbol
pattern = "qam16"
# 模型参数
device = torch.device("cpu:0")
BATCH_SIZE = 2000
NUM_EPOCHS = 200000
train_num = BATCH_SIZE
test_num = 20000
learning_rate = 1e-4
USE_CUDA = True

class bit_AE(nn.Module):
    def __init__(self, ):
        super(bit_AE, self).__init__()
        self.encoder4 = nn.Sequential(
            nn.Linear(10, 4096),
            nn.ReLU(),
            nn.Linear(4096, 4096),
            nn.ReLU(),
            nn.Linear(4096, 4096),
            nn.ReLU(),
            nn.Linear(4096, 4096),
            nn.ReLU(),
            nn.Linear(4096, 4096),
            nn.ReLU(),
            nn.Linear(4096, 2)
        )
        self.decoder = nn.Sequential(
            nn.Linear(2, 4096),
            nn.ReLU(),
            nn.Linear(4096, 4096),
            nn.ReLU(),
            nn.Linear(4096, 4096),
            nn.ReLU(),
            nn.Linear(4096, 4096),
            nn.ReLU(),
            nn.Linear(4096, 4096),
            nn.ReLU(),
            nn.Linear(4096, 10)
        )

    def decode_signal(self, x):
        return self.decoder(x)

    def encode_signal(self, x):
        return self.encoder(x)

    def forward(self,b,snr = 7):
        x = self.encoder4(b.to(device))
        batchsize = len(b)

        x_norm = x/torch.sqrt(torch.sum(x ** 2) / batchsize)
        snr_w = 10.0 ** (snr / 10.0)
        noise = torch.randn([batchsize,2]) / ((2* bits_per_symbol * snr_w) ** 0.5)
        x_noise = x_norm + noise.to(device)
        x_ = self.decoder(x_noise)

        return x_,x_norm

def train(path , snr , use_pre_model):
    model = bit_AE()
    pre_state_dict = torch.load(path +'square_qam1024_0dB_goodbase.pt')
    model.load_state_dict(pre_state_dict, strict=0)
    # 重新训练 or 在之前基础上训练
    if use_pre_model > 0:
        model.load_state_dict(torch.load(path + "qam1024_{}dB_goodbase.pt".format(use_pre_model)))
    model.to(device)
    ber = 1
    train_labels = torch.randint(low = 0,high=2,size=[train_num,bits_per_symbol],dtype=torch.float)
    train_data = train_labels
# DataBase in Pytorch
    dataset = Data.TensorDataset( train_data,  train_labels)
    train_loader = Data.DataLoader(dataset = dataset, batch_size = BATCH_SIZE, shuffle = True, num_workers = 0)
#optmizer & Loss
    optimizer = Adam(model.parameters(),lr=learning_rate)
    loss_fn = nn.BCEWithLogitsLoss()
#Training
    for epoch in range(NUM_EPOCHS):
       for step, (x, y) in enumerate(train_loader):
            y = y.float()
            x = x.to(device)
            y = y.to(device)
            decoded,_ = model(x,snr)
            loss = loss_fn(decoded, y)
            optimizer.zero_grad()               # clear gradients for this training step
            loss.backward()                     # backpropagation, compute gradients
            optimizer.step()                    # apply gradients
            error_bits = torch.sum(torch.abs((torch.sign(decoded)+1)/2 - y))
            if error_bits/BATCH_SIZE / bits_per_symbol < ber :
                ber = error_bits/BATCH_SIZE / bits_per_symbol
                path = path
                torch.save(model, path + pattern + 'SNR{}.pth'.format(snr))
                torch.save(model.state_dict(), path +'qam1024_{}dB_goodbase.pt'.format(snr))
                print("BER" , ber)
                # 画星座图
                index = torch.linspace(0, modulation_order - 1, modulation_order).long()
                bits = [bin(i)[2:].zfill(bits_per_symbol) for i in index]
                test_data = np.zeros([len(bits), bits_per_symbol])
                for i in range(len(bits)):
                    for j in range(bits_per_symbol):
                        test_data[i, j] = int(bits[i][j])
                test_data = torch.tensor(test_data, dtype=torch.float)
                _, x = model(test_data.to(device), snr)
                plot_data = x.to("cpu").data.numpy()
                plot_Constellation(path+"result//",plot_data,snr, epoch,modulation_order)
       if epoch%100 ==0 :
           print('Epoch: ', epoch, 'loss: ', loss.item() )
    file_name = 'qam1024_AE.mat'
    savemat(file_name, {'complex_values':plot_data })
    return ber

def plot_Constellation(path ,plot_data,snr, epoch,modulation_order):
    plt.figure()
    plt.scatter(plot_data[:, 0], plot_data[:, 1])
    plt.axis((-2, 2, -2, 2))
    plt.grid()
    if epoch > 0:
        plt.savefig(path + "qam{}epoch{}_snr{}dB_bit_wise.png".format(modulation_order,epoch,snr))

    plt.close()
def frange(x, y, jump):
    while x < y:
        yield x
        x += jump
def test(path ,mode, snr):
    # path = 'E://DF//Deep-Learning-for-the-Physical-Layer-master//bit_AE_model//'
    if mode ==0:  # 展示星座图
        path = path
        model = torch.load(path + pattern + 'SNR{}.pth'.format(snr))

        index = torch.linspace(0, modulation_order - 1,modulation_order).long()
        bits = [bin(i)[2:].zfill(bits_per_symbol) for i in index]
        test_data = np.zeros([len(bits) , bits_per_symbol])
        for i in range(len(bits)):
            for j in range(bits_per_symbol):
                test_data[i,j] = int(bits[i][j])
        test_data = torch.tensor(test_data , dtype=torch.float)
        _,x = model(test_data.to(device),snr)
        plot_data = x.to("cpu").data.numpy()
        plot_Constellation(path ,plot_data, snr, 0, modulation_order)
        # plt.show()


if __name__ == "__main__":
    train_SNR = 12
    path = 'D:\\graduation\\qam1024\\'
    r = 2/3
    use_pre_model = 0
    mode = 0
    model_snr = 12
   # test(path, mode, model_snr)

    train(path, train_SNR,use_pre_model)



