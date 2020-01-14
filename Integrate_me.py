## check what the uncertainty is as a function of sampled points n
# demonstrate convergence, perform timing tests


class Bin:
    def __init__(self, a, b, n = 0, f = None):
        self.a = a
        self.b = b
        self.sum = 0
        self.sum2 = 0

        self.islist = False

        if type(self.a) == list:
            self.islist = True
            self.dim = len(a)
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


        if self.islist == False:
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

        else:
            for _ in range(n): #points
                coords = []
                for d in range(self.dim):
                    x = random.uniform(self.a[d], self.b[d])
                    coords.append(x)
                    #print(x)
                
                newtuple = tuple(coords)
                #fname = f.__name__
                fstring = 'self.f{}'.format(newtuple)
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
        if self.islist == False:
            midpoint = (self.a + self.b) / 2
            A = Bin(self.a, midpoint, self.n, self.f)
            B = Bin(midpoint, self.b, self.n, self.f)
            A.MC()
            B.MC()
            return A,B
        
        else:
            midpoint = [(self.a[i] + self.b[i]) / 2  for i in range(self.dim)]
        
            Vb = [0 for i in range(self.dim)]
            REC = Integrator().rec(2, self.dim, Vb)
            newbins = []
            
            amb = [self.a, midpoint, self.b]

            try:
                while True:
                    Vnow = next(REC)
                    start = [amb[Vnow[d]][d] for d in range(self.dim)]
                    fin = [amb[Vnow[d] + 1][d] for d in range(self.dim)]
                    newbin = Bin(start, fin, self.n, self.f)
                    newbin.MC()
                    newbins.append(newbin)
            except:
                pass

            return newbins


'''
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
        REC = Integrator().rec(2, self.dim, Vb)
        newbins = []
        
        amb = [self.a, midpoint, self.b]

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
'''

class Integrator:
    '''
    A class used to numerically integrate one or multi-dimensional functions using numerical methods
    ...

    Attributes
    ----------
    f : function 
        The function meant to be integrated
    i,j : int
        A diagnostic attribute. i and j represents the number of calls AdaptInt1 and AdaptInt2 does, respectively.
    '''
    def __init__(self, f = None):
        '''
        Initializer of the function. f is the function meant to be integrated. 
        '''
        self.f = f
        #self.i = 0
        self.j = 0
    
    def __repr__(self):
        '''
        returns an executable string that represents this class
        '''
        newtuple = tuple([self.f])
        classname = self.__class__.__name__
        return '{}{}'.format(classname, newtuple)

    def midpoint(self, a, b):
        '''
        Integrates a one dimensional function using the midpoint rule. 

        Attributes
        ----------
        a : float 
            lower limit of function
        b : float
            upper limit of function
        '''

        R = b-a
        return R*self.f((b+a)/2)

    def cmidpoint(self, a, b, N):
        '''
        Integrates a one dimensional function using the composite midpoint rule. 

        Attributes
        ----------
        a : float 
            lower limit of function
        b : float
            upper limit of function
        N : int
            Number of  intervals for the function to calculate
        '''
        f = self.f
        h = (b-a)/N
        x = [a + i*h for i in range(N+1)]
        value =  0

        #xin = [] #unused
        for i in range(N):
            xa = x[i]
            xb = x[i+1]
            value += (xb - xa) * f((xb+xa)/2)
            #xin.append((xb+xa)/2)
        
        return value

    def trapz(self, a, b):
        '''
        Integrates a one dimensional function using the trapezium rule. 

        Attributes
        ----------
        a : float 
            lower limit of function
        b : float
            upper limit of function
        '''

        f = self.f
        R = (b-a)/2
        return R*(f(a) + f(b))

    def ctrapz(self, a, b, N):
        '''
        Integrates a one dimensional function using the composite trapezium rule. 

        Attributes
        ----------
        a : float 
            lower limit of function
        b : float
            upper limit of function
        N : int
            Number of  intervals for the function to calculate
        '''
        f = self.f
        h = (b-a)/N
        x = [a + i*h for i in range(N+1)]
        value = 0
        value += h*f(x[0])/2

        for i in range(1, N):
            #print(i)
            value += f(x[i]) * h
        value += h*f(x[N])/2
        return value

    def simpsons(self, a, b):
        '''
        Integrates a one dimensional function using the simpsons rule. 

        Attributes
        ----------
        a : float 
            lower limit of function
        b : float
            upper limit of function
        '''
        f = self.f
        R = (b-a)/6
        return R * (f(a) + f(b) + 4*f((a+b)/2))
    
    def csimpsons(self, a, b, N):
        '''
        Integrates a one dimensional function using the composite simpsons rule. 

        Attributes
        ----------
        a : float 
            lower limit of function
        b : float
            upper limit of function
        N : int
            Number of  intervals for the function to calculate
        '''
        f = self.f
        h = (b-a)/N
        R = h/3
        x = [a + i*h for i in range(N+1)]
        weight = []
        for i in range(N):
            if i % 2 == 0:
                weight.append(2)
            else:
                weight.append(4)
        weight[0] = 1
        weight.append(1)
        
        value = 0 
        for i in range(N+1):
            value += R * weight[i] * f(x[i])
        return value

    def retanalysis(self, a, b, Nmax, Ndiffs, intmethod):
        '''
        Performs time, convergence, and accuracy test as a function of intervals used
        Uses the Newton-Cotes integration method

        Attributes
        ----------
        a : float 
            lower limit of function
        b : float
            upper limit of function
        Nmax : int
            Number of maximum intervals for the function to calculate
        Ndiffs: int
            Differences between successive intervals tested. eg. if Nmax = 5 and Ndiffs = 2, then those tested 
            are N = 1, 3, 5
        '''
        import time
        start_time = time.time()

        nPoints = [2]
        intvalue = [self.NCInt(a, b, nPoints[0], intmethod)]
        time_taken = [time.time()-start_time]

        for i in range(int(Nmax/Ndiffs)):
            npoint = (i+1)*Ndiffs
            
            start_time = time.time()
            tintvalue = self.NCInt(a, b, npoint, intmethod)
            timepast = time.time() - start_time

            nPoints.append(npoint)
            intvalue.append(tintvalue)
            time_taken.append(timepast)

        return [nPoints, intvalue, time_taken]

    def plotme(self, a, b, intmethod, Nmax = 500, Ndiffs = 10, realvalue = None):
        '''
        Performs time, convergence, and accuracy test as a function of intervals used
        Uses the Newton-Cotes integration method
        Plots results in two or three separate figures, depending on if realvalue was specified.

        Attributes
        ----------
        a : float 
            lower limit of function
        b : float
            upper limit of function
        intmethod : 0, 1, 2
            Method of choice for integration, per the NCInt method. 
        Ndiffs: int
            Differences between successive intervals tested. eg. if Nmax = 5 and Ndiffs = 2, then those tested 
            are N = 1, 3, 5
        realvalue: float
            The real value of the integral. Used to compare the numerical integration with the analytical one.
        '''
        import matplotlib.pyplot as plt 
        
        RetAn = self.retanalysis(a, b, Nmax, Ndiffs, intmethod)

        plt.figure(0)
        plt.title('Value of Integral as a function of sample points')
        plt.plot(RetAn[0], RetAn[1])
        plt.show()

        if realvalue != None:
            plt.figure(1)
            plt.title('Percentage Difference from Actual Value')
            RVlst = [(realvalue - i)*100 / realvalue for i in RetAn[1]]
            plt.plot(RetAn[0], RVlst)
            plt.show()
            
        plt.figure(2)
        plt.title('Timing Test')
        plt.plot(RetAn[0], RetAn[2])
        plt.show()
        #return nPoints

    def plotmeval(self, a, b, intmethod, Nmax = 500, Ndiffs = 10, realvalue = None):
        '''
        Performs time, convergence, and accuracy test as a function of intervals used
        Uses the Newton-Cotes integration method
        Plots results in one figure.

        Attributes
        ----------
        a : float 
            lower limit of function
        b : float
            upper limit of function
        intmethod : 0, 1, 2
            Method of choice for integration, per the NCInt method. 
        Ndiffs: int
            Differences between successive intervals tested. eg. if Nmax = 5 and Ndiffs = 2, then those tested 
            are N = 1, 3, 5
        realvalue: float
            The real value of the integral. Used to compare the numerical integration with the analytical one.
        '''
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec

        RetAn = self.retanalysis(a, b, Nmax, Ndiffs, intmethod)
        
        fig = plt.figure(tight_layout=True)
        gs = gridspec.GridSpec(2, 2)
        ax = fig.add_subplot(gs[0, :])
        ax.plot(RetAn[0], RetAn[2])
        ax.set_title('%s for intmethod %s' % (self.f.__name__, intmethod))
        ax.set_ylabel('Time taken (s)')
        ax.set_xlabel('Number of Sample Points')

        if realvalue != None:
            ax = fig.add_subplot(gs[1,1])
            RVlst = [(realvalue - i)*100 / realvalue for i in RetAn[1]]
            
            ax.plot(RetAn[0], RVlst)
            ax.set_title('Percentage Diffs')
            ax.set_ylabel('Percentage Difference (%)')
            ax.set_xlabel('Number of Sample Points')

            ax = fig.add_subplot(gs[1,0])
        else: 
            ax = fig.add_subplot(gs[1,:])
        
        ax.plot(RetAn[0], RetAn[1])
        ax.set_title('Integral Values')
        ax.set_ylabel('Value of Integral')
        ax.set_xlabel('Number of Sample Points')

        fig.align_labels()
        plt.show()

    def NCInt(self, a, b, N, kind):
        '''
        Integrates a one dimensional function using the Newton-Cotes Rule. 

        Attributes
        ----------
        a : float 
            lower limit of function
        b : float
            upper limit of function
        N : int
            Number of  intervals for the function to calculate
        kind : int
            Determines the weights used for the integration
            0 = Midpoint Rule
            1 = Trapezium Rule
            2 = Simpsons Rule
        '''
        value = 0
        h = (b-a) / N
        wlist = self.w(N, kind, h)
        x = [a + i*h for i in range(N+1)]

        #Dealing with midpoint rule, as x needs a change of coordinates
        if kind == 0:
            xnew = []
            for i in range(N):
                xa = x[i]
                xb = x[i+1]
                xnew.append((xb+xa)/2)
            x = xnew
            N -= 1

        for i in range(N+1):
            value += wlist[i] * self.f(x[i])
        
        return value

    def w(self, N, kind, h):
        '''
        Returns the weights used for the Newton-Cotes integration rule. Depends on the type of integration desired. 

        Attributes
        ----------
        N : int
            Number of  intervals for the function to calculate
        kind : int
            Determines the weights used for the integration
            0 = Midpoint Rule
            1 = Trapezium Rule
            2 = Simpsons Rule
        h : float
            The distance between two successive intervals
        '''

        #could also be done with a switch case
        if kind == 0: #Midpoint Rule
            wlist = [h for i in range(N)]

        elif kind == 1: #Trapezium Rule
            wlist = [h for i in range(N)]
            wlist[0] = h/2
            wlist.append(h/2)

        elif kind == 2: #Simpsons Rule
            wlist = []
            for i in range(N):
                if i % 2 == 0:
                    wlist.append(2 * h / 3)
                else:
                    wlist.append(4 * h / 3)
            wlist[0] = 1 * h/3
            wlist.append(1 * h/3)

        return wlist

    def AdaptInt(self, a, b, tau, intmeth):
        '''
        Integrate using the Adaptive Integration Method, using the Newton-Cotes Rule. 
        Uses the the difference between two divisions of the same level to calculate error

        Attributes
        ----------
        a : float 
            lower limit of function
        b : float
            upper limit of function
        tau : float
            Error Tolerance
        intmeth : int
            Integration method of choice, from NCInt
        Q0 : float
            Variable used for recursion, compares current value with value from previous iteration
        '''
        N = 1
        if intmeth == 2:
            N = 2
        self.j += 1

        val = self.NCInt(a, b, N, intmeth) #this needs to be fixed because it needs to take range
        m = (a+b)/2
        val1 = self.NCInt(a, m, N, intmeth)
        val2 = self.NCInt(m, b, N, intmeth)

        #err = val1 - val2
        err = val - (val1 + val2)
        if abs(err) > tau:
            if err > 0:
                val = self.AdaptInt(a, m, tau, intmeth)
                val += self.AdaptInt(m, b, tau, intmeth)
            else:
                val = self.AdaptInt(m, b, tau, intmeth)
                val += self.AdaptInt(a, m, tau, intmeth)
        return val
    
    def MonteCarlo(self, a, b, n = 1000):
        '''
        Method that performs the Monte Carlo Integration
    
        Attributes
        ----------
        a : float 
            lower limit of function
        b : float
            upper limit of function
        n : int
            Number of random points used
        '''

        import random
        random.seed(1) #used for reproducibility

        f = self.f
        fvalue = 0
        for _ in range(n):
            x = random.uniform(a, b)
            fvalue += f(x)

        favg = fvalue/n
        return favg * (b-a)

    def MonteCarloN(self, A, B, n = 1000):
        '''
        Method that performs the Monte Carlo Integration in N dimensions
    
        Attributes
        ----------
        a : float 
            lower limit of function
        b : float
            upper limit of function
        n : int
            Number of random points used
        '''
        MCN = Bin(A, B, n, f = self.f)
        MCN.MC()
        return MCN.val

    #Test me
    def rec(self, N, d, V, x = 0):
        '''
        Generator Function used for getting the multidimensional weight indexer.
        
        N = number of points in each dimension as a list
        d = number of dimensions (len(A))
        V = list of indexes at for this iteration
        '''

        if type(N) != list:
            Nnew = [N for i in range(len(A))]
            N = Nnew

        if x != d-1:
            y = 1 * x
            x += 1    
            for i in range(N[y]):
                V[y] = i
                yield from self.rec(N, d, V, x)
        else:
            for i in range(N[x]):
                V[x] = i
                yield V

    def NCIntN(self, A, B, N, kind):
        """
        Method that performs the Newton-Cotes Integration in multiple dimensions
        All lists must have the same number of elements, which correspond to the total dimensions 
        
        Attributes
        ----------
        A : list of floats 
            lower limits of function, each element corresponds to value at one dimension
        B : list of floats
            upper limits of function, each element corresponds to value at one dimension
        N : list of ints
            Number intervals in a specific dimension.
        kind : int
            Number corresponding to the type of Newton-Cotes integration needed.
        """

        def wfval(V, Wlist, f, X):
            '''
            Function that returns the value of a multidimensional function and the weights
            Given the indexer list V

            Attributes
            ----------
            V : List of ints
                Indexer list. Each integer shows which weight is used for that particular dim
                From the rec generator
            Wlist : List of Lists
                The Weight lists, from 
            '''
            thevalue = 1
            coords = []
            for o in range(len(V)):
                thevalue *= Wlist[o][V[o]]
                coords.append(X[o][V[o]])
            newtuple = tuple(coords)
            #fname = f.__name__
            fstring = 'self.f{}'.format(newtuple)
            #fval = eval(fstring)
            #return [thevalue, fval]
            return [thevalue, fstring]

        '''
        #Variables
        A = [0, 0, 0]
        B = [10, 10, 10]
        N = [100, 100, 100]
        kind = 2 #trapz 
        '''
        #Initialized number of dimensions and indexer V
        d = len(N)
        V = [0 for i in range(d)]
        f = self.f
        value = 0

        dmul = 1
        for i in range(d):
            dmul *= N[i]

        H = [(B[i]-A[i]) / N[i] for i in range(d)]
        Wlist = [self.w(N[i], kind, H[i]) for i in range(d)]
        X = [[A[i] + j*H[i] for j in range(N[i]+1)] for i in range(d)]

        #deals with the cmidpoint rule
        if kind == 0: 
            for i in range(d):
                x = X[i]
                xnew = []
                for j in range(N[i]):
                    xa = x[j]
                    xb = x[j+1]
                    xnew.append((xb+xa)/2)
                X[i] = xnew
                N[i] -= 1

        gen = self.rec(N, d, V)

        for _ in range(dmul):
            Vnow = next(gen)        
            valnow = wfval(Vnow, Wlist, f, X)
            value += valnow[0]*eval(valnow[1])
            
        return value

    def StratSamp(self, a, b, Nbin, Ninbin, Nintcheck = 10, MaxVar = 10):
        """
        Method that performs one-dimensional Stratified Sampling MC.
        
        Attributes
        ----------
        a : float
            lower limit
        b : float
            upper limit
        Nbin : int
            Number of initial bins used to divide the range
        Ninbin : int
            Number of points used to perform MC in each bin once bins have been divided
        Nintcheck : int
            Number of points used to perform MC used to estimate the variance in the bin
        MaxVar : float
            Maximum variance in each bin. If the bin's variance exceeds this, bin will be subdivided
        """ 
        N = Nbin
        Nint = Ninbin

        import random
        random.seed(1) #used for reproducibility
        '''
        #variables
        a = 0
        b = 10
        N = 10 #number of bins
        Nintcheck = 10 #used to estimate bin size
        Nint = 1000 #number of points per bin was 100000
        MaxVar = 10
        '''

        h = (b-a) / N
        x = [a + i*h for i in range(N+1)]

        Aval = []
        Avar = []
        BinList = []

        #generates initial bins
        for i in range(len(x)-1):
            A = Bin(x[i], x[i+1], Nintcheck, self.f)
            A.MC()
            BinList.append(A)
            Aval.append(A.val)
            Avar.append(A.var)
        

        #Stratified Sampling alg
        while max(Avar) > MaxVar:
            maxind = Avar.index(max(Avar))
            newbins = BinList[maxind].Bisect()
            BinList[maxind] = newbins[0]
            Avar[maxind] = newbins[0].var
            BinList.insert(maxind + 1, newbins[1])
            Avar.insert(maxind + 1, newbins[1].var)

        finalval = 0
        for i in BinList:
            i.MC(n = Nint)
            finalval += i.val

        return finalval

    def StratSampN(self, a, b, Nbin, Ninbin, Nintcheck = 10, MaxVar = 10):
        """
        Method that performs multi-dimensional Stratified Sampling MC.
        
        Attributes
        ----------
        A : List
            each element is the lower limits for a particular dimension
        B : List
            each element is the lower limits for a particular dimension
        Nbin : int
            Number of initial bins used to divide the range
        Ninbin : int
            Number of points used to perform MC in each bin once bins have been divided
        Nintcheck : int
            Number of points used to perform MC used to estimate the variance in the bin
        MaxVar : float
            Maximum variance in each bin. If the bin's variance exceeds this, bin will be subdivided
        """ 

        
        N = Nbin
        Nint = Ninbin

        '''
        #variables
        N = 4 #number of bins per dim
        Nintcheck = 10 #used to estimate bin size
        Nint = 1000 #number of points per bin
        MaxVar = 10
        '''


        d = len(A)
        V = [0 for i in range(d)]

        import random
        random.seed(1) #used for reproducibility

        H = [(B[i]-A[i]) / N for i in range(d)]
        X = [[A[i] + j*H[i] for j in range(N+1)] for i in range(d)]


        Aval = []
        Avar = []
        BinList = []

        REC = self.rec(N, d, V)
        
        #Initial division of bins
        totrecs = N**d
        for _ in range(totrecs):
            Vnow = next(REC)
            Acoord = []
            Bcoord = []
            for o in range(d):
                Acoord.append(X[o][Vnow[o]])
                Bcoord.append(X[o][Vnow[o] + 1])

            thebin = Bin(Acoord, Bcoord, Nintcheck, self.f) 
            thebin.MC()
            BinList.append(thebin)
            Aval.append(thebin.val)
            Avar.append(thebin.var)

        #Stratified Sampling Alg
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

        return finalval

    def AdaptIntN(self, A, B, tau, intmeth):
        dim= len(A)
        N = 1
        if intmeth == 2:
            N = 2
        self.j += 1

        Nproper = [N for i in range(dim)]
        val = self.NCIntN(A, B, Nproper, intmeth) #this needs to be fixed because it needs to take range
        #print('val', val)
        #Bisection
        M = [(A[i] + Nproper[i]) / 2  for i in range(dim)]
        
        Vb = [0 for i in range(dim)]
        REC = self.rec(2, dim, Vb)
        amb = [A, M, B]
        newvals = []
        startlist = []
        finlist = []
        #print(amb)
        try:
            while True:
                Vnow = next(REC)
                #print(Vnow)
                start = [amb[Vnow[d]][d] for d in range(dim)]
                #print('start', start)
                startlist.append(start)
                fin = [amb[Vnow[d] + 1][d] for d in range(dim)]
                finlist.append(start)
                
                newval = self.AdaptIntN(start, fin, tau, intmeth)
                newvals.append(newval)
        except:
            pass


        #err = max(newvals) - min(newvals)
        err = sum(newvals)
        #print('err', newvals)
        
        diff = val - (err)
        print('diff, tau', diff, tau)
        if abs(diff) > tau:
            val = 0
            for i in range(len(newvals)):
                val += self.AdaptIntN(startlist[i], finlist[i], tau, intmeth)

        return val




if __name__ == '__main__':
    def f1(x):
        return x**1

    def f2(x):
        return x**2

    def f3(x):
        return x**3 

    def f4(x):
        return x**4

    def f5(x):
        return x**5
    def g(x, y, z):
        return x + y + z


    def functiontester(xmin, xmax, VAL):
        flist = [f1, f2, f3, f4, f5]
        for func in flist:
            intme = Integrator(func)
            for i in range(3):
                ind = flist.index(func)
                Rval = VAL[ind]
                intme.plotmeval(xmin, xmax, i, realvalue = Rval )


    a = 0 
    b = 10
    VAL = [50, 1000/3, 2500, 20000, 500000/3]

    #functiontester(a, b, VAL)

    
    
    intme = Integrator(f4)
    ltest1 = intme.w(10, 1, 1)
    ltest2 = intme.w(10, 2, 1)
    import time

    print('small - small')
    start_time = time.time()
    print(intme.AdaptInt(0, 10, 1e-5, 2))
    print(time.time()-start_time)
    print('calls', intme.j)
    print('')
    print(intme.NCIntN([0], [10], [1000], 2))
    #print('large - small')
    #start_time = time.time()
    #print('value', intme.AdaptInt(0, 10, 1e-2, 2))
    #print('time taken', time.time()-start_time)
    #print('calls', intme.i)

    A = [0, 0, 0]
    B = [10, 10, 10]
    N = [100, 100, 100]
    kind = 2
    testme = Integrator(g)
    #NCINT = testme.NCIntN(A, B, N, 2)
    #print(NCINT)
    
    #ADAPTN = testme.AdaptIntN(A, B, 100, 2)
    #print(ADAPTN)