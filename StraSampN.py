import time
start_time = time.time()

def g(x, y):
    R = x**2 + y**2
    if R <= 1:
        return 1
    else:
        return 0

def f(x, y):
    value = (-(x**2 + y**2) + 1) * 100
    if value < 0:
        value = 0
    return value

class binN:
    def __init__(self, a, b, n = 0, f = None):
        self.a = a
        self.b = b
        self.sum = 0
        self.sum2 = 0
        self.n = n
        self.dim = len(a)
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

        for _ in range(n): #points
            coords = []
            for d in range(self.dim):
                x = random.uniform(self.a[d], self.b[d])
                coords.append(x)
                #print(x)
            
            newtuple = tuple(coords)
            fname = f.__name__
            fstring = '{}{}'.format(fname, newtuple)
            fval = eval(fstring)
            

            self.sum += fval
            self.sum2 += fval**2


        
        self.avg = self.sum / self.n
        self.avg2 = self.sum2 / self.n
        self.var = self.avg2 - (self.avg ** 2)
        diffs = [self.b[i] - self.a[i] for i in range(self.dim)]
        #print(diffs)
        self.val = self.avg
        for delta in diffs:
            self.val *= delta

    def Bisect(self):
        midpoint = [(self.a[i] + self.b[i]) / 2  for i in range(self.dim)]
    
        Vb = [0 for i in range(self.dim)]
        REC = rec(2, self.dim, Vb)
        newbins = []
        
        amb = [self.a, midpoint, self.b]

        '''
        a [-5, -5]
        m [0, 0]
        b [5, 5]

        make amb = [[-5, -5], [0, 0] ,[5, 5]]
        
        want [amb[rec[d]][d], amb[rec[d]][d]] and [amb[rec[d] + 1][d], amb[rec[d] + 1][d]] 
        [-5, -5], [0, 0] for [0, 0] => [amb[0][0], amb[0][1]] and [amb[1][0], amb[1][1]]
        [-5, 0], [0, -5] for [0, 1] => [amb[0][0], amb[1][1]] and [amb[1][0], amb[2][1]]
        [0, -5], [5, 0] for [1, 0] => [amb[1][0], amb[0][1]] and [amb[2][0], amb[1][1]]
        [0, 0], [5, 5] for [1, 1] => [amb[1][0], amb[1][1]] and [amb[2][0], amb[2][1]]
        '''

        try:
            while True:
                Vnow = next(REC)
                start = [amb[Vnow[d]][d] for d in range(self.dim)]
                fin = [amb[Vnow[d] + 1][d] for d in range(self.dim)]
                newbin = binN(start, fin, self.n, self.f)
                newbin.MC()
                newbins.append(newbin)
        except:
            pass

        return newbins




def rec(N, d, V, x = 0):
    '''
    N = number of points in each dimension
    d = number of dimensions (len(A))
    V = list of indexes at for this iteration
    '''
    if x != d-1:
        y = 1 * x
        x += 1    
        for i in range(N):
            V[y] = i
            yield from rec(N, d, V, x)
    else:
        for i in range(N):
            V[x] = i
            yield V
            #V[x] = i
            #can now stick a function in here related to V or make it a generator
            #print(V) #now i just need to do something with this
            #return V






#A = [-5, -5, -5]
#B = [5, 5, 5]


A = [-1, -1]
B = [1, 1]
V = [0, 0]


BIN = binN(A, B, 10000, f)
BIN.MC()

print(BIN.val)
#untested territory

#variables
N = 4 #number of bins per dim
Nmax = 10 #number of bins per dim max 
Nintcheck = 10 #used to estimate bin size
Nint = 1000 #number of points per bin
MaxVar = 10

#A = bin(-1, -0.8)
#A.smth(f, n = 10)


d = len(A)

import random
random.seed(1) #used for reproducibility

H = [(B[i]-A[i]) / N for i in range(d)]
X = [[A[i] + j*H[i] for j in range(N+1)] for i in range(d)]


Aval = []
Avar = []
BinList = []


V = [0 for i in range(d)]

REC = rec(N, d, V)

acoord = [X[0][0], X[1][0]] #[X[dimension][indexer]]
bcoord = [X[0][1], X[1][1]] 
totrecs = N**d #was d**N

for _ in range(totrecs):
    Vnow = next(REC)
    Acoord = []
    Bcoord = []
    for o in range(d):
        Acoord.append(X[o][Vnow[o]])
        Bcoord.append(X[o][Vnow[o] + 1])

    thebin = binN(Acoord, Bcoord, Nintcheck, f) 
    thebin.MC()
    BinList.append(thebin)
    Aval.append(thebin.val)
    Avar.append(thebin.var)




while max(Avar) > MaxVar:
    maxind = Avar.index(max(Avar))
    newbins = BinList[maxind].Bisect()
    newvars = [i.var for i in newbins]

    del BinList[maxind]
    del Avar[maxind]

    BinList += newbins
    Avar += newvars






finalval = 0
for i in BinList:
    i.MC(n = Nint)
    finalval += i.val


TotVal = sum(Aval)
print('TotVal', TotVal)
print('Final Value', finalval)


#can i make a bin class?
time_taken = time.time()-start_time
print('timetaken', time_taken)
