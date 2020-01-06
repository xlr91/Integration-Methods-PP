#import heartrate; heartrate.trace(browser=True)
import time
start_time = time.time()

def f(x):
    return x**3


class Bin:
    def __init__(self, a, b, n = 0, f = None):
        self.a = a
        self.b = b
        self.sum = 0
        self.sum2 = 0
        self.n = n #no of points in this bin
        self.xval = []

        self.avg = 0 
        self.avg2 = 0
        self.var = 0
        self.f = f
    
    def __repr__(self):
        newtuple = tuple([self.a, self.b, self.n, self.f])
        classname = self.__class__.__name__
        return '{}{}'.format(classname, newtuple)
    
    def MC(self, f = None, n = None):
        import random
        random.seed(1)

        if n == None:
            n = self.n
        else:
            self.n = n

        if f == None:
            f = self.f
        else:
            self.f = f

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
        self.val = self.avg * (self.b-self.a)

    def Bisect(self):
        midpoint = (self.a + self.b) / 2
        A = Bin(self.a, midpoint, self.n, self.f)
        B = Bin(midpoint, self.b, self.n, self.f)
        A.MC()
        B.MC()
        return A,B




#variables
a = 0
b = 10
N = 10 #number of bins
Nmax = 10 #number of bins max
Nintcheck = 10 #used to estimate bin size
Nint = 100000 #number of points per bin
MaxVar = 10

#A = bin(-1, -0.8)
#A.smth(f, n = 10)




import random
random.seed(1) #used for reproducibility
h = (b-a) / N
x = [a + i*h for i in range(N+1)]

Aval = []
Avar = []
BinList = []

#generates first general bins
for i in range(len(x)-1):
    A = Bin(x[i], x[i+1], Nintcheck, f)
    A.MC()
    BinList.append(A)
    Aval.append(A.val)
    Avar.append(A.var)
    #print(x[i], x[i+1])

#does the stratified sampling bit
while max(Avar) > MaxVar:
    maxind = Avar.index(max(Avar))
    newbins = BinList[maxind].Bisect()
    BinList[maxind] = newbins[0]
    Avar[maxind] = newbins[0].var
    BinList.insert(maxind + 1, newbins[1])
    Avar.insert(maxind + 1, newbins[1].var)

j = 0
for i in BinList:
    i.MC(n = Nint)
    j += i.val


TotVal = sum(Aval)
print('TotVal', TotVal)
print('j', j)


#can i make a bin class?
time_taken = time.time()-start_time
print('timetaken', time_taken)