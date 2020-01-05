def f(x):
    return x**3


class bin:
    def __init__(self, a, b, n = 0):
        self.a = a
        self.b = b
        self.sum = 0
        self.sum2 = 0
        self.n = n
    
    def smth(self, f, n = None):
        import random
        random.seed(1)

        if n == None:
            n = self.n
        else:
            self.n = n

        self.sum = 0
        self.sum2 = 0
        self.xval = []

        for _ in range(n): 
            x = random.uniform(self.a, self.b)
            self.xval.append(x)
            fv = f(x)
            self.sum += fv
            self.sum2 += fv ** 2
        
        self.avg = self.sum / self.n
        self.avg2 = self.sum2 / self.n
        self.var = self.avg2 - (self.avg ** 2)


#variables
a = 0
b = 10
N = 10 #number of bins
Nmax = 100 #number of bins max
Nint = 10 #number of points per bin

#A = bin(-1, -0.8)
#A.smth(f, n = 10)




import random
random.seed(1) #used for reproducibility
h = (b-a) / N
x = [a + i*h for i in range(N+1)]

Aavg = []
Avar = []

for i in range(len(x)-1):
    A = bin(x[i], x[i+1])
    A.smth(f, n = Nint)
    Aavg.append(A.avg)
    Avar.append(A.var)
    #print(x[i], x[i+1])
    

TotAvg = sum(Aavg) * (b-a)/N
print('TotAvg', TotAvg)


#can i make a bin class?
