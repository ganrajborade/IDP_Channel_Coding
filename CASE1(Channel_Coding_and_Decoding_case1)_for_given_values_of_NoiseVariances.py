import numpy as np
from matplotlib import pyplot as plt
from pylab import *
img_array = np.load('mss.npy')

m = img_array.reshape(-1)  #Converting 400x300 array to 1D-array of size 120000.

T = 10**(-6)  #T = 1 microseconds
fc = 2*(10**(6)) # Carrier frequency = 2MHz

fs = 50*(10**(6)) # Sampling rate of fs =  4MHz.
nm = m.size
n=8
k=4
c_size = nm*n/k


#Channel Encoding :
C = np.zeros((30000,8))
G= [[1,1,1,1,0,0,0,0],[1,1,0,0,1,1,0,0],[1,0,1,0,1,0,1,0],[0,1,1,0,1,0,0,1]]
G_T = np.transpose(G)
for i in range(0,int(nm/(k))):
    A = np.zeros(4)
    for j in range(4):
        A[j] = m[4*i + j]
        S = np.matmul(G_T,A)
        C[i] = S
C = C.reshape(-1)
for i in range(240000):
    C[i] = C[i]%2  #making each element in C[i] as either 0 or 1.

b = C


def PixelsToCoordinates(pixels):
    x = np.zeros(len(pixels))
    for i in range(len(pixels)):
        if (pixels[i] == 0):
            x[i] = 1
        elif(pixels[i] == 1):
            x[i] = -1
    return x

def TransmittedSignal(x) :
    f_c = np.int(2*(1e6))
    f_s = np.int(50*(1e6))
    t_s = 1/f_s
    s = np.zeros(50*120000)


    for n in range(len(s)) :
        s[n] = (x[np.int(2*(math.floor(n/50)))] * np.cos(2*(np.pi)*f_c*(n)*t_s)) + (x[np.int((2*(math.floor(n/50)))+1)] * np.sin(2*(np.pi)*f_c*(n)*t_s))
    return s

# Adding noise..
def Noise(N_0,s) :
    w = np.random.normal(0, np.sqrt(((N_0) )), 50*120000)
    s = s + w
    return s


# Ideal Signals..
def IdealizedSignals(a1,a2,a3):
    f_c = np.int(2*(1e6))
    f_s = np.int(50*(1e6))
    t_s = 1/f_s

    s = np.zeros(50)
    for n in range(len(s)):
        s[n] = a1*(np.cos(2*(np.pi)*f_c*(50*(a3) + n)*t_s)) + a2*(np.sin(2*(np.pi)*f_c*(50*(a3) + n)*t_s))
    return s


#Mean_square_error
def MeanSquareError(a,b):
    sum = 0
    for e in range(50):
        sum = ((a[e] - b[e])**2) + sum
    return (sum**0.5)


#  Demodulation...
def MinDistance(r):
    sym = []
    a = np.zeros(50)
    b = np.zeros(50)
    f_c = np.int(2 * (1e6))
    f_s = np.int(50 * (1e6))
    t_s = 1 / f_s

    for k in range(120000):
        l = []
        s1 = IdealizedSignals(1, 1, k)
        s2 = IdealizedSignals(1, -1, k)
        s3 = IdealizedSignals(-1, 1, k)
        s4 = IdealizedSignals(-1, -1, k)
        a = r[(50*k):((50*k) + 50)]

        l.append(MeanSquareError(a,s1))
        l.append( MeanSquareError(a, s2))
        l.append(MeanSquareError(a, s3))
        l.append(MeanSquareError(a, s4))

        sym.append(l.index((min(l))))
    return sym


#Demodulation
def Getting1DMSS(sym):
    mss = np.zeros(240000)
    for i in range(120000):
        if(sym[i] == 0):
            mss[(2*i)] = 0
            mss[(2*i)+1] = 0

        elif (sym[i] == 1):
            mss[(2 * i)] = 0
            mss[(2 * i) + 1] = 1

        elif (sym[i] == 2):
            mss[(2 * i)] = 1
            mss[(2 * i) + 1] = 0

        elif (sym[i] == 3):
            mss[(2 * i)] = 1
            mss[(2 * i) + 1] = 1
    return mss

x = PixelsToCoordinates(b)

s = TransmittedSignal(x)

r1 = Noise(5,s)
r2= Noise(7,s)
r3 = Noise(12,s)
r4 = Noise(20,s)



sym1 = MinDistance(r1)
sym2 = MinDistance(r2)
sym3 = MinDistance(r3)
sym4 = MinDistance(r4)



mss1 = Getting1DMSS(sym1)
mss2 = Getting1DMSS(sym2)
mss3 = Getting1DMSS(sym3)
mss4 = Getting1DMSS(sym4)



y1 = mss1
y2 = mss2
y3 = mss3
y4 = mss4



n1 =0
m1 = np.zeros((16,4))
for i in range(2):
    for j in range(2):
        for k in range(2):
            for l in range(2):
                m1[n1] = np.array([i,j,k,l])
                n1 = n1+1
               
def Hamming_Distance(a,b):
    if(len(a) == len(b)):
        i = 0
        for j in range(len(a)):
            if(a[j] != b[j]):
                i = i+1
    return i

#Channel Decoding
def ChannelDecoder(y):
	m_reformed = np.zeros(120000)
	for t in range(30000):
	    w = np.zeros(16)
	    for i in range(16):
	        testing = np.matmul(G_T,m1[i])
	        testing = testing%2
	        w[i] = Hamming_Distance(testing,y[8*t:(8*t)+8])
	    
	    index = np.argmin(w)
	    m_reformed[4*t:(4*t)+4] = m1[index]
	return m_reformed



def bitErrorRate(m_reformed,m):
	if (m_reformed.size == m.size):


		total_error_bits = 0

		for i in range(m.size):
		        	if(m_reformed[i] != m[i]):

		           		total_error_bits +=1
		return total_error_bits/len(m)


g1 = ChannelDecoder(y1)
g2 = ChannelDecoder(y2)
g3 = ChannelDecoder(y3)
g4 = ChannelDecoder(y4)


ber1 = bitErrorRate(g1,m)
ber2 = bitErrorRate(g2,m)
ber3 = bitErrorRate(g3,m)
ber4 = bitErrorRate(g4,m)


print("Bit Error Rate at variance of noise = "+str(5)+" is ",ber1)
print("Bit Error Rate at variance of noise= "+str(7)+" is ",ber2)
print("Bit Error Rate at variance of noise = "+str(12)+" is ",ber3)
print("Bit Error Rate at variance of noise = "+str(20)+" is ",ber4)


#Final Image
image_array_received_1 = g1.reshape(400,300)
image_array_received_2 = g2.reshape(400,300)
image_array_received_3 = g3.reshape(400,300)
image_array_received_4 = g4.reshape(400,300)


subplot(2,3,1)
plt.imshow(image_array_received_1,'gray')
plt.title("Received Image at variance of noise= "+str(5))

subplot(2,3,2)
plt.imshow(image_array_received_2,'gray')
plt.title("Received Image at variance of noise = "+str(7))

subplot(2,3,3)
plt.imshow(image_array_received_3,'gray')
plt.title("Received Image at variance of noise= "+str(12))

subplot(2,3,4)
plt.imshow(image_array_received_4,'gray')
plt.title("Received Image at variance of noise = "+str(20))


subplot(2,3,5)
plt.semilogy([5,7,12,20],[ber1,ber2,ber3,ber4])
plt.xlabel('Noise Variance')
plt.ylabel('BER')
plt.grid()
plt.show()