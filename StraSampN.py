def f(x, y):
    R = x**2 + y**2 
    if R <= 1:
        return 1
    else:
        return 0

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
        print(diffs)
        self.val = self.avg
        for delta in diffs:
            self.val *= delta

    def Bisect(self):
        midpoint = (self.a + self.b) / 2
        A = bin(self.a, midpoint, self.n, self.f)
        B = bin(midpoint, self.b, self.n, self.f)
        A.MC()
        B.MC()
        return A,B


A = [-1, -1]
B = [1, 1]

BIN = binN(A, B, 100000, f)
BIN.MC()
