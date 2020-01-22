## check what the uncertainty is as a function of sampled points n
# demonstrate convergence, perform timing tests


class RSSBin:
    """
    A Bin class for Recursive Stratified Sampling Monte Carlo. 
    
    Attributes
    ----------
    a,b : float or list of floats
        The lower and upper limits of the bins
        If list, each element is the limit for that respective dimension
    n : int
        Number of sampling points used to Monte Carlo the bin
    f : function
        Function that is to be integrated in the bin
    sum : float
        The sum of the function calls after random sampling
    sum2 : float
        The sum squared of the function calls after random sampling, for variance calculation
    islist : Boolean
        Internal attribute to decide if the Bin is one dimensional or multi dimensional
    dim : int
        The number of dimensions the Bin has. Equal to length of a 
    xval : list
        List of x sampling points chosen, used for diagnostic purposes
    avg : float
        Average value of bin
    avg2 : float
        The squared average of bin <f(x)^2>
    var : float
        Variance of bin
    """
    def __init__(self, a, b, n = 0, f = None):
        """
        Initializes the class
        """
        self.a = a
        self.b = b
        self.n = n 
        self.f = f

        self.islist = False
        if type(self.a) == list:
            self.islist = True
            self.dim = len(a)
        
        self.xval = []
        self.sum = 0
        self.sum2 = 0
        self.avg = 0 
        self.avg2 = 0
        self.var = 0
        
    def __repr__(self):
        """
        Returns the representation of the class
        """
        newtuple = tuple([self.a, self.b, self.n, self.f])
        classname = self.__class__.__name__
        return '{}{}'.format(classname, newtuple)
    
    def MC(self, f = None, n = None):
        """
        Performs a Monte Carlo Integration for the specified function and specified Range
        Attributes
        ----------
        f : function
            The function that is meant to be integrated
        n : int
            Number of random sampling points used
        """
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

        #One dimensoinal Case
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

        #Multi Dimensional Case
        else:
            for _ in range(n):
                coords = []
                for d in range(self.dim):
                    x = random.uniform(self.a[d], self.b[d])
                    coords.append(x)
                newtuple = tuple(coords)
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
        """
        Returns a list of new Bins, resulted from dividing the current bin in two in all dimensions
        """
        #One dimensional case
        if self.islist == False:
            midpoint = (self.a + self.b) / 2
            A = RSSBin(self.a, midpoint, self.n, self.f)
            B = RSSBin(midpoint, self.b, self.n, self.f)
            A.MC()
            B.MC()
            return A,B
        
        #Multi Dimensional case
        else:
            midpoint = [(self.a[i] + self.b[i]) / 2  for i in range(self.dim)]
            Vb = [0 for i in range(self.dim)]
            REC = Integrator().rec(2, self.dim, Vb)
            newbins = []
            amb = [self.a, midpoint, self.b]

            #Loops over all possible bisection combinations
            try:
                while True:
                    Vnow = next(REC)
                    start = [amb[Vnow[d]][d] for d in range(self.dim)]
                    fin = [amb[Vnow[d] + 1][d] for d in range(self.dim)]
                    newbin = RSSBin(start, fin, self.n, self.f)
                    newbin.MC()
                    newbins.append(newbin)
            except:
                pass

            return newbins

class Integrator:
    """
    A class used to numerically integrate one or multi-dimensional functions using numerical methods
    ...

    Attributes
    ----------
    f : function 
        The function meant to be integrated
    j : int
        A diagnostic attribute. j is used to count how many recursions have been done in methods requiring recursion
    vaar : float
        A diagnostic attribute. The variance returned in the intergration after Monte Carlo was performed.
    stratcalls : int
        Number of stratified sampling calls used thus far 
    """
    def __init__(self, f = None):
        """
        Initializer of the function. f is the function meant to be integrated. 
        """
        self.f = f
        self.j = 0
        self.vaar = 0
    
    def __repr__(self):
        """
        Returns the representation of the class
        """
        newtuple = tuple([self.f])
        classname = self.__class__.__name__
        return '{}{}'.format(classname, newtuple)

    def midpoint(self, a, b):
        """
        Integrates a one dimensional function using the midpoint rule. 

        Attributes
        ----------
        a : float 
            lower limit of function
        b : float
            upper limit of function
        """

        R = b-a
        return R*self.f((b+a)/2)

    def cmidpoint(self, a, b, N):
        """
        Integrates a one dimensional function using the composite midpoint rule. 

        Attributes
        ----------
        a : float 
            lower limit of function
        b : float
            upper limit of function
        N : int
            Number of  intervals for the function to calculate
        """
        f = self.f
        h = (b-a)/N
        x = [a + i*h for i in range(N+1)]
        value =  0

        for i in range(N):
            xa = x[i]
            xb = x[i+1]
            value += (xb - xa) * f((xb+xa)/2)
        return value

    def trapz(self, a, b):
        """
        Integrates a one dimensional function using the trapezium rule. 

        Attributes
        ----------
        a : float 
            lower limit of function
        b : float
            upper limit of function
        """
        f = self.f
        R = (b-a)/2
        return R * (f(a) + f(b))

    def ctrapz(self, a, b, N):
        """
        Integrates a one dimensional function using the composite trapezium rule. 

        Attributes
        ----------
        a : float 
            lower limit of function
        b : float
            upper limit of function
        N : int
            Number of  intervals for the function to calculate
        """
        f = self.f
        h = (b-a)/N
        x = [a + i*h for i in range(N+1)]
        value = 0
        value += h*f(x[0])/2

        for i in range(1, N):
            value += f(x[i]) * h
        value += h*f(x[N])/2
        return value

    def simpsons(self, a, b):
        """
        Integrates a one dimensional function using the simpsons rule. 

        Attributes
        ----------
        a : float 
            lower limit of function
        b : float
            upper limit of function
        """
        f = self.f
        R = (b-a)/6
        return R * (f(a) + f(b) + 4*f((a+b)/2))
    
    def csimpsons(self, a, b, N):
        """
        Integrates a one dimensional function using the composite simpsons rule. 

        Attributes
        ----------
        a : float 
            lower limit of function
        b : float
            upper limit of function
        N : int
            Number of  intervals for the function to calculate
        """
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

    def NCInt(self, a, b, N, kind):
        """
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
        """
        value = 0
        h = (b-a) / N
        wlist = self.w(N, kind, h)
        x = [a + i*h for i in range(N+1)]

        #Midpoint rule requiring a change of coordinates
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
        """
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
        """

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
        """
        Integrate using the Adaptive Integration Method, using the Newton-Cotes Rule. 

        Attributes
        ----------
        a : float 
            lower limit of function
        b : float
            upper limit of function
        tau : float
            Error Tolerance of integration value
        intmeth : int
            Integration method of choice, from NCInt
        """

        N = 1
        if intmeth == 2:
            N = 2
        self.j += 1

        val = self.NCInt(a, b, N, intmeth)

        m = (a+b)/2
        val1 = self.NCInt(a, m, N, intmeth)
        val2 = self.NCInt(m, b, N, intmeth)
        err = val - (val1 + val2)

        if abs(err) > tau:
            if err > 0:
                val = self.AdaptInt(a, m, tau/2, intmeth)
                val += self.AdaptInt(m, b, tau/2, intmeth)
            else:
                val = self.AdaptInt(m, b, tau/2, intmeth)
                val += self.AdaptInt(a, m, tau/2, intmeth)

        return val
    
    def MonteCarlo(self, a, b, n = 1000):
        """
        Method that performs the Monte Carlo Integration
    
        Attributes
        ----------
        a : float 
            lower limit of function
        b : float
            upper limit of function
        n : int
            Number of random points used
        """

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
        """
        Method that performs the Monte Carlo Integration in N dimensions
    
        Attributes
        ----------
        a : float 
            lower limit of function
        b : float
            upper limit of function
        n : int
            Number of random points used
        """

        #utilises the RSSBin for convenience
        self.vaar = 0
        MCN = RSSBin(A, B, n, f = self.f)
        MCN.MC()
        self.vaar = MCN.var
        return MCN.val

    def rec(self, N, d, V, x = 0):
        """
        Generator function used to return the Cartesian products of all the indexes used in each dimension.
        
        Attributes
        ----------
        N : List of integers
            number of points in each dimension as a list (eg. x has 3 points, y has 4 points, and so on)
        d : int
            number of dimensions (len(N), eg. (x, y) has 2 dimensions))
        V : List 
            list of indexes at for this iteration. Starts at [0, ..., 0]
        x : int
            Internal variable for counting iterations
        """

        if type(N) != list:
            Nnew = [N for i in range(d)]
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
            """
            Function that returns the value of a multidimensional function and the weights
            Given the indexer list V

            Attributes
            ----------
            V : List of ints
                Indexer list. Each integer shows which weight is used for that particular dim
                From the rec generator
            Wlist : List of Lists
                The Weight lists, from 
            
            Returns : float, string
                Float is the value of the weight products, string is the string representation of the 
                function meant to be integrated
            """
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

        #Initialized number of dimensions and indexer V
        d = len(N)
        V = [0 for i in range(d)]
        f = self.f
        value = 0


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

        try:
            while True:
                Vnow = next(gen)        
                valnow = wfval(Vnow, Wlist, f, X)
                value += valnow[0]*eval(valnow[1])
        except:
            pass
            
        return value

    def StratSamp(self, a, b, Nbin = 4, Ninbin = 1000, Nintcheck = 10, MaxVar = 10, MaxIter = None):
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
        MaxIter : int
            Maximum divisions for the algorithm
        """ 
        N = Nbin
        Nint = Ninbin
        self.stratcalls = 0

        import random
        random.seed(1) 

        h = (b-a) / N
        x = [a + i*h for i in range(N+1)]

        Aval = [] #used for diagnostic purposes
        Avar = []
        BinList = []

        #generates initial bins
        for i in range(len(x)-1):
            A = RSSBin(x[i], x[i+1], Nintcheck, self.f)
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

            #MaxIter condition
            self.stratcalls += 1
            try:
                if self.stratcalls >= MaxIter:
                    break
            except:
                pass

        finalval = 0
        for i in BinList:
            i.MC(n = Nint)
            finalval += i.val

        return finalval

    def StratSampN(self, A, B, Nbin = 4, Ninbin = 1000, Nintcheck = 10, MaxVar = 1e-2, MaxIter = None):
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
        self.stratcalls = 0
        d = len(A)
        V = [0 for i in range(d)]

        import random
        random.seed(1) #used for reproducibility

        H = [(B[i]-A[i]) / N for i in range(d)]
        X = [[A[i] + j*H[i] for j in range(N+1)] for i in range(d)]

        Aval = [] #Diagnostic Purposes
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

            thebin = RSSBin(Acoord, Bcoord, Nintcheck, self.f) 
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
            self.stratcalls += 1
            
            #MaxIter Condition
            try:
                if self.stratcalls >= MaxIter:
                    break
            except:
                pass


        finalval = 0
        for i in BinList:
            i.MC(n = Nint)
            finalval += i.val

        return finalval

    def AdaptIntN(self, A, B, tau, intmeth, N = 2):
        """
        Integrate in multiple dimensions using the Adaptive Integration Method, using the Newton-Cotes Rule. 

        Attributes
        ----------
        A : list of floats
            lower limit of function
        B : list of floats
            upper limit of function
        tau : float
            Error Tolerance of integration value
        intmeth : int
            Integration method of choice, from NCInt
        N : int
            Number of points used for the Newton-Cotes integration
        """

        dim= len(A)
        self.j += 1

        Nproper = [N for i in range(dim)]
        val = self.NCIntN(A, B, Nproper, intmeth) 
        Vb = [0 for i in range(dim)]
        REC = self.rec(2, dim, Vb)

        M = [(A[i] + B[i]) / 2  for i in range(dim)]
        amb = [A, M, B]
        newvals = []
        startlist = []
        finlist = []

        #Division Algorithm
        try:
            while True:
                Vnow = next(REC)

                start = [amb[Vnow[d]][d] for d in range(dim)]
                startlist.append(start)

                fin = [amb[Vnow[d] + 1][d] for d in range(dim)]
                finlist.append(fin)

                newval = self.NCIntN(start, fin, Nproper, intmeth)
                newvals.append(newval)
        except:
            pass
        
        sums = sum(newvals)
        diff = val - (sums)
        if abs(diff) > tau:
            val = 0
            for i in range(len(newvals)):
                val += self.AdaptIntN(startlist[i], finlist[i], tau, intmeth, N)
        return val


class Integrator_Analysis:
    """
    A class used to analyse and generate plots of the Integrator Methods

    """
    def __init__(self):
        """
        Initialises the Algorithm
        """
        pass

    def plotme(self, RetAn, realvalue = None, ylim = None):
        """
        Plots the results of the supplied analysis 

        Attributes
        ----------
        RetAn : List
            Results of analysis from of of the multiple retanalysis functions
        realvalue : float
            The real (known) value of an integral. Used to compare the estimate with the theoretical value. 
        ylim : 2 element list
            The limits of the percentage difference graph, in the form [ymin, ymax]
        """
        import matplotlib.pyplot as plt 
        import numpy as np
        
        #Integral value vs sample points per dim
        plt.figure(0)
        plt.title('Value of Integral as a function of sample points')
        plt.plot(RetAn[0], RetAn[1])
        plt.grid()
        plt.xlabel('Sample Points n')
        plt.ylabel('Int. Value')


        #Percentage Diff vs sample points per dim
        if realvalue != None:
            plt.figure(1)
            plt.title('Percentage Difference from Actual Value')
            plt.grid()
            RVlst = [-(realvalue - i)*100 / realvalue for i in RetAn[1]]
            plt.plot(RetAn[0], RVlst, label='%diff')
            if RetAn[3] != None:
                plt.plot(RetAn[3][0], RetAn[3][1], label = 'upper bound', alpha=0.25, color = 'k')
                plt.plot(RetAn[3][0], -RetAn[3][1], label = 'lower bound', alpha=0.25, color = 'k')
            plt.ylim(-10, 10)
            
            if ylim != None:
                plt.ylim(ylim[0], ylim[1])
            
            plt.xlabel('Sample Points n')
            plt.ylabel('Error (%)')
            plt.legend()
                    
        #Timing Tests vs sample points per dim
        plt.figure(2)
        plt.grid()
        plt.title('Timing Test')
        plt.xlabel('Sample Points n per dimension')
        plt.ylabel('Time to Execute (s)')
        plt.plot(RetAn[0], RetAn[2])
        #plt.show()
        
        #Methods involving some form of recursion
        #Timing Tests vs Sample Points total and %diff vs calls
        if RetAn[4] != None:
            
            plt.figure(0)
            plt.clf()
            plt.title('Value of Integral as a function of Iterations')
            plt.grid()
            plt.plot(RetAn[4], RetAn[1])
            plt.xlabel('Iterations of Algorithm')
            plt.ylabel('Int. Value')  
        
            plt.figure(2)
            plt.clf()
            plt.grid()
            plt.plot(RetAn[4], RetAn[2])
            plt.title('Percentage Difference from Actual Value')
            plt.xlabel('Iterations of Algorithm')
            plt.ylabel('Time to Execute (s)')
            plt.title('Timing Test')
            
            if realvalue!= None:
                plt.figure(1)
                plt.clf()
                plt.grid()
                plt.ylim(-10, 10)
                plt.plot(RetAn[4], RVlst)
                plt.xlabel('Iterations of Algorithm')
                if ylim != None:
                    plt.ylim(ylim[0], ylim[1])
                plt.ylabel('Error (%)')
                plt.title('% Difference vs Number of Iterations of Alg')
        
        plt.show()

    def PLOTME(self, RETAN, LABEL, REALVALUE = None, ylim = None):
        """
        Plots the results of the supplied analysis. Used to compare the results of two or more methods

        Attributes
        ----------
        RETAN : List of lists
            Results of analysis from of of the multiple retanalysis functions
            Supplied in the form [RetAn1, RetAn2, ...] 

        LABEL : List of strings
            Labels of each of the analysis from the code
            Supplied in the form ['label for RetAn1', 'label for RetAn2', ...]

        REALVALUE : List of floats
            The real (known) value of an integral. Used to compare the estimate with the theoretical value. 
            Supplied in the form of [realvalue1, realvalue2, ...]
        ylim : 2 element list
            The limits of the percentage difference graph, in the form [ymin, ymax]
        """

        import matplotlib.pyplot as plt 
        import numpy as np
        
        
        #Integral value vs sample points per dim
        plt.figure(0)
        plt.title('Value of Integral as a function of sample points')
        for i in range(len(RETAN)):
            plt.plot(RETAN[i][0], RETAN[i][1], label= LABEL[i])
        plt.grid()
        plt.xlabel('Sample Points n')
        plt.ylabel('Int. Value')
        plt.legend()

        #Percentage Diff vs sample points per dim
        if REALVALUE != None:
            plt.figure(1)
            plt.title('Percentage Difference from Actual Value')
            plt.grid()
            RVlst = [[-(REALVALUE[j] - i)*100 / REALVALUE[j] for i in RETAN[j][1]] for j in range(len(REALVALUE))]
            for i in range(len(RETAN)):
                plt.plot(RETAN[i][0], RVlst[i], label= LABEL[i])
            plt.ylim(-10, 10)        
            if ylim != None:
                plt.ylim(ylim[0], ylim[1])
            plt.xlabel('Sample Points n per dimension')
            plt.ylabel('Error (%)')
            plt.legend()
        
        #Timing Tests vs sample points per dim
        plt.figure(2)
        plt.grid()
        plt.title('Timing Test')
        plt.xlabel('Sample Points n')
        plt.ylabel('Time to Execute (s)')
        for i in range(len(RETAN)):
            plt.plot(RETAN[i][0], RETAN[i][2], label= LABEL[i])
        plt.legend()
        
        #Timing Tests vs Sample Points total and %diff vs calls
        if RETAN[0][4] != None:
            plt.figure(0)
            plt.clf()
            
            plt.title('Value of Integral as a function of Iterations')
            plt.grid()
            
            for i in range(len(RETAN)):
                plt.plot(RETAN[i][4], RETAN[i][1], label= LABEL[i])
            plt.xlabel('Iterations of Algorithm')
            plt.legend()
            plt.ylabel('Int. Value')  
        
            plt.figure(2)
            plt.clf()
            plt.grid()
            for i in range(len(RETAN)):
                plt.plot(RETAN[i][4], RETAN[i][2], label= LABEL[i])
            plt.xlabel('Iterations of Algorithm')
            plt.ylabel('Time to Execute (s)')
            plt.legend()
            plt.title('Timing Test')
            
            if REALVALUE!= None:
                plt.figure(1)
                plt.clf()
                plt.grid()
                plt.ylim(-10, 10)
                for i in range(len(RETAN)):
                    plt.plot(RETAN[i][4], RVlst[i], label= LABEL[i])
                plt.legend()
                plt.xlabel('Iterations of Algorithm')
                if ylim != None:
                    plt.ylim(ylim[0], ylim[1])
                plt.ylabel('Error (%)')
                plt.title('% Difference vs Number of Iterations of Alg')
        
        plt.show()

    def retanalysis(self, C, a, b, intmethod, Nmax = 500, Ndiffs = 10):
        """
        Performs time, convergence, and accuracy test as a function of intervals used
        Uses the Newton-Cotes integration method

        Attributes
        ----------
        C: Integrator Object
            The object of integration containing the function.
        a : float 
            lower limit of function
        b : float
            upper limit of function
        intmethod : 0, 1, 2
            The method of NC integration as specified from NCInt
        Nmax : int
            Number of maximum intervals for the function to calculate
        Ndiffs: int
            Differences between successive intervals tested. eg. if Nmax = 5 and Ndiffs = 2, then those tested 
            are N = 1, 3, 5

        Returns : List
            [nPoints, intvalue, time_taken, error bounds [n, Y] ,j (adaptcalls)]
        """
        import time
        import numpy as np
        
        start_time = time.time()

        nPoints = [1]
        intvalue = [C.NCInt(a, b, nPoints[0], intmethod)]
        time_taken = [time.time()-start_time]

        for i in range(int(Nmax/Ndiffs)):
            npoint = (i+1)*Ndiffs

            start_time = time.time()
            tintvalue = C.NCInt(a, b, npoint, intmethod)
            timepast = time.time() - start_time

            nPoints.append(npoint)
            intvalue.append(tintvalue)
            time_taken.append(timepast)
        
        n = np.arange(4, nPoints[-1], 0.5)
        if intmethod == 2: 
            Y = ((b-a)**5) *100 / (180 * n**4)
        elif intmethod == 1:
            Y = ((b-a)**3) *100 / (12 * n**2)
        else:
            Y = ((b-a)**3) *100 / (24 * n**2)
        return [nPoints, intvalue, time_taken, [n,Y], None]

    def plotmecombi(self, C, a, b, Nmax = 500, Ndiffs = 10, realvalue = None):
        """
        Analyses and plots the results of all three NCInt method
        Uses the Simpsons Rule only

        Attributes
        ----------
        C: Integrator Object
            The object of integration containing the function.
        a : float 
            lower limit of function
        b : float
            upper limit of function
        Nmax : int
            Number of maximum intervals for the function to calculate
        Ndiffs: int
            Differences between successive intervals tested. eg. if Nmax = 5 and Ndiffs = 2, then those tested 
            are N = 1, 3, 5
        realvalue: float
            The real value of the integral. Used to compare the numerical integration with the analytical one.
        """
        import matplotlib.pyplot as plt 

        RetAn = [self.retanalysis(C, a, b, i, Nmax, Ndiffs) for i in range(3)]

        plt.figure(0)
        plt.grid()
        plt.title('Value of Integral as a function of sample points')
        plt.xlabel('Sample Points n')
        plt.ylabel('Int. Value')
        for i in range(3):
            plt.plot(RetAn[i][0], RetAn[i][1], label = 'Int. Approx. Order ' + str(i))
        plt.legend()
        plt.show()

        if realvalue != None:
            plt.figure(1)
            plt.grid()
            plt.title('Percentage Difference from Actual Value')
            for i in range(3):
                RVlst = [-(realvalue - i)*100 / realvalue for i in RetAn[i][1]]
                plt.plot(RetAn[i][0], RVlst, label = 'Int. Approx. Order ' + str(i))
            plt.legend()
            plt.xlabel('Sample Points n')
            plt.ylabel('Error (%)')
            plt.ylim(-10, 10)
            plt.show()

        plt.figure(2)
        plt.grid()
        plt.title('Timing Test')
        for i in range(3):
            plt.plot(RetAn[i][0], RetAn[i][2], label = 'Int. Approx. Order ' + str(i))
        plt.legend()
        plt.xlabel('Sample Points n')
        plt.ylabel('Time to Execute (s)')
        plt.show()
    
    def retanalysisadapt(self, C, a, b, Taumin10to, intmethod = 2):
        """
        Performs time, convergence, and accuracy test as a function of intervals used
        Uses the adaptive quadrature integration method

        Attributes
        ----------
        C: Integrator Object
            The object of integration containing the function.
        a : float 
            lower limit of function
        b : float
            upper limit of function
        Taumin10to : int
            The minimum tolerance for adaptive quadrature, 10^-(Taumin10to)
        intmethod : 0, 1 , 2
            Type of NC integration used
        
        Returns : List
            [nPoints, intvalue, time_taken, error bounds [n, Y] ,j (adaptcalls)]
        """
        
        import time
        start_time = time.time()
        
        C.j = 0
        nPoints = [1000]
        intvalue = [C.AdaptInt(a, b, nPoints[0], intmethod)]
        J = [C.j]
        time_taken = [time.time()-start_time]

        for i in range(int(3-Taumin10to)):
            npoint = 1000/(2**i)
            C.j = 0
            start_time = time.time()
            tintvalue = C.AdaptInt(a, b, npoint, intmethod)
            timepast = time.time() - start_time
            J.append(C.j)
            nPoints.append(npoint)
            intvalue.append(tintvalue)
            time_taken.append(timepast)

        return [nPoints, intvalue, time_taken, None, J]

    def retanalysismonteN(self, C, a, b, Nmax = 500, Ndiffs = 10):
        """
        Performs time, convergence, and accuracy test as a function of intervals used
        Uses the N dimensional monte carlo method

        Attributes
        ----------
        C: Integrator Object
            The object of integration containing the function.
        a : float 
            lower limit of function
        b : float
            upper limit of function
        Nmax : int
            Number of maximum intervals for the function to calculate
        Ndiffs: int
            Differences between successive intervals tested. eg. if Nmax = 5 and Ndiffs = 2, then those tested 
            are N = 1, 3, 5
        
        Returns : List
            [nPoints, intvalue, time_taken, error bounds [n, Y] ,j (adaptcalls)]
        """
        import time
        import numpy as np
        start_time = time.time()

        nPoints = [10]
        intvalue = [C.MonteCarloN(a, b, nPoints[0])]
        time_taken = [time.time()-start_time]
        std = [C.vaar ** (1/2)]

        for i in range(int(Nmax/Ndiffs)):
            npoint = (i+1)*Ndiffs

            start_time = time.time()
            tintvalue = C.MonteCarloN(a, b, npoint)
            timepast = time.time() - start_time

            nPoints.append(npoint)
            intvalue.append(tintvalue)
            time_taken.append(timepast)
            std.append(C.vaar**(1/2))

        return [nPoints, intvalue, time_taken, None, None]

    def retanalysisstrat(self, C, a, b, Ninbins = 1000, Nmax = 500, Ndiffs = 10):
        """
        Performs time, convergence, and accuracy test as a function of intervals used
        Uses the Recursive Stratified Sampling Monte Carlo method

        Attributes
        ----------
        C: Integrator Object
            The object of integration containing the function.
        a : float 
            lower limit of function
        b : float
            upper limit of function
        Ninbins : int
            Number of sampling points per bin once bins have finished subdividing
        Nmax : int
            Maximum number of iterations the algorithm can run
        Ndiffs : int
            Interval between successive number of algorithm the program will run    
        Returns : List
            [nPoints, intvalue, time_taken, error bounds [n, Y] ,j (adaptcalls)]
        """

        import time
        start_time = time.time()
        
        nPoints = [1]
        intvalue = [C.StratSamp(a, b, Ninbin = Ninbins, MaxVar = 1e-5, MaxIter = nPoints[0])]
        J = [C.stratcalls]
        time_taken = [time.time()-start_time]

        for i in range(int(Nmax/Ndiffs)):
            npoint = (i+1)*Ndiffs
            C.i = 0
            start_time = time.time()
            tintvalue = C.StratSamp(a, b, Ninbin = Ninbins, MaxVar = 1e-5, MaxIter = npoint)
            timepast = time.time() - start_time
            J.append(C.stratcalls)
            nPoints.append(npoint)
            intvalue.append(tintvalue)
            time_taken.append(timepast)

        return [nPoints, intvalue, time_taken, None, J]

    def retanalysisN(self, C, a, b, intmethod, Nmax = 30, Ndiffs = 10):
        """
        Performs time, convergence, and accuracy test as a function of intervals used
        Uses the Newton-Cotes integration method in N dimensions

        Attributes
        ----------
        C: Integrator Object
            The object of integration containing the function.
        a : list of floats 
            lower limit of functions
        b : list of floats
            upper limit of functions
        intmethod : 0, 1, 2
            The method of NC integration as specified from NCInt
        Nmax : int
            Number of maximum intervals for the function to calculate
        Ndiffs: int
            Differences between successive intervals tested. eg. if Nmax = 5 and Ndiffs = 2, then those tested 
            are N = 1, 3, 5

        Returns : List
            [nPoints, intvalue, time_taken, error bounds [n, Y] ,j (adaptcalls)]
        """
        import time
        import numpy as np
        start_time = time.time()
        
        nPoints = [2]
        Npoints = [nPoints[0] for i in range(len(a))]
        intvalue = [C.NCIntN(a, b, Npoints, intmethod)]
        time_taken = [time.time()-start_time]

        for i in range(int(Nmax/Ndiffs)):
            npoint = (i+1)*Ndiffs

            start_time = time.time()
            Npoint = [npoint for i in range(len(a))]
            tintvalue = C.NCIntN(a, b, Npoint, intmethod)
            timepast = time.time() - start_time

            nPoints.append(npoint)
            intvalue.append(tintvalue)
            time_taken.append(timepast)
            
        return [nPoints, intvalue, time_taken, None , None]

    def retanalysisadaptN(self, C, a, b, Taumin10to, intmethod = 2):
        """
        Performs time, convergence, and accuracy test as a function of intervals used
        Uses the adaptive quadrature integration method in N dimensions

        Attributes
        ----------
        C: Integrator Object
            The object of integration containing the function.
        a : list of floats 
            lower limit of functions
        b : list of floats
            upper limit of functions
        Taumin10to : int
            The minimum tolerance for adaptive quadrature, 10^-(Taumin10to)
        intmethod : 0, 1 , 2
            Type of NC integration used
        
        Returns : List
            [nPoints, intvalue, time_taken, error bounds [n, Y] ,j (adaptcalls)]
        """
        import time
        start_time = time.time()
        
        C.j = 0
        nPoints = [1000]
        intvalue = [C.AdaptIntN(a, b, nPoints[0], intmethod, N = 2)]
        J = [C.j]
        time_taken = [time.time()-start_time]

        for i in range(int(3-Taumin10to)):
            npoint = 1000/10**i
            C.j = 0
            start_time = time.time()
            tintvalue = C.AdaptIntN(a, b, npoint, intmethod, N = 2)
            timepast = time.time() - start_time
            J.append(C.j)
            nPoints.append(npoint)
            intvalue.append(tintvalue)
            time_taken.append(timepast)

        return [nPoints, intvalue, time_taken, None, J]

    def retanalysisstratN(self, C, A, B, Ninbins = 1000, Nmax = 100, Ndiffs = 10):
        """
        Performs time, convergence, and accuracy test as a function of intervals used
        Uses the Recursive Stratified Sampling Monte Carlo method in N dimensions

        Attributes
        ----------
        C: Integrator Object
            The object of integration containing the function.
        A : list of floats 
            lower limit of functions
        B : list of floats
            upper limit of functions
        Ninbins : int
            Number of sampling points per bin once bins have finished subdividing
        Nmax : int
            Maximum number of iterations the algorithm can run
        Ndiffs : int
            Interval between successive number of algorithm the program will run    
        Returns : List
            [nPoints, intvalue, time_taken, error bounds [n, Y] ,j (adaptcalls)]
        """
        import time
        start_time = time.time()
        
        nPoints = [1]
        intvalue = [C.StratSampN(A, B, Ninbin = Ninbins, MaxVar = 1e-5, MaxIter = nPoints[0])]
        J = [C.stratcalls]
        time_taken = [time.time()-start_time]

        for i in range(int(Nmax/Ndiffs)):
            npoint = (i+1)*Ndiffs
            C.i = 0
            start_time = time.time()
            tintvalue = C.StratSampN(A, B, Ninbin = Ninbins, MaxVar = 1e-5, MaxIter = npoint)
            timepast = time.time() - start_time
            J.append(C.stratcalls)
            nPoints.append(npoint)
            intvalue.append(tintvalue)
            time_taken.append(timepast)

        return [nPoints, intvalue, time_taken, None, J]


if __name__ == '__main__':
    def f1(x):
        return x 

    def f2(x):
        return x**2

    def f3(x):
        return x**3

    def f5(x):
        return x+x**2 - x**3 + x**5

    def f4(x, y):
        return x+y

    def Gauss(x):
        return 2.718281828459045**(-x**2) 

    def S(x, y, z, w):
        value = 0
        if x**2 + y**2 + z**2 + w**2 <=1:
            value = 1
        return value

    def circ(x, y, z, w ):
        return x**2 + y**2 + z**2 + w**2
    

    F1 = Integrator(f1) 
    F2 = Integrator(f2) 
    F3 = Integrator(f3) 
    F4 = Integrator(f4) 
    F5 = Integrator(f5)
    G = Integrator(Gauss)
    SPh = Integrator(S)
    Circ = Integrator(circ)

    IA = Integrator_Analysis()

    RetAnG = IA.retanalysis(G, -10, 10, 2, Ndiffs = 1, Nmax = 200)
    IA.plotme(RetAnG, realvalue = 1.7724538509)

    print('F3, NC Integration')    
    IA.plotmecombi(F3,0, 10, Nmax = 50, realvalue = 2500, Ndiffs = 2)

    print('F5, NC Integration')
    IA.plotmecombi(F5, -2, 6, Nmax = 50, realvalue = 7536, Ndiffs = 2)

    print('Gaussian, NC Integration')
    IA.plotmecombi(G,-5, 5, Nmax = 50, realvalue = 1.7724538509, Ndiffs = 1)

    RetAd1 = IA.retanalysisadapt(F5, -2, 6, -20, 2)
    RetAd2 = IA.retanalysisadapt(G, -10, 10, -24, 2)

    print('Combined Adaptive Integration')
    IA.PLOTME([RetAd1, RetAd2], LABEL = ['F5', 'Gaussian'], REALVALUE = [7536, 1.7724538509], ylim = [-6, 2.5])


    print('F5, Adaptive Quadrature')
    IA.plotme(RetAd1, realvalue = 7536, ylim = [-0.2, 2])

    print('Gauss Adaptive')
    IA.plotme(RetAd2, realvalue = 1.7724538509)

    RetMC1 = IA.retanalysismonteN(G, -5, 5, Nmax = 5000, Ndiffs = 10)
    RetMC2 = IA.retanalysismonteN(F5, -2, 6, Nmax = 20000)

    print('MC Gaussian')
    IA.plotme(RetMC1, realvalue = 1.7724538509)

    print('MC F5')
    IA.plotme(RetMC2, realvalue = 7536)

    RetSt1 = IA.retanalysisstrat(F5, -2, 6)
    RetSt2 = IA.retanalysisstrat(G, -10, 10)

    print('F5 stratified')
    IA.plotme(RetSt1, realvalue = 7536, ylim = [-0.2, 2])

    print('F5 gaussian')
    IA.plotme(RetSt2, realvalue = 1.7724538509, ylim = [-0.2, 2])

    IA.PLOTME([RetSt1, RetSt2], LABEL=['test', 'test'], REALVALUE= [7536, 1.7724538509])

    RetAnN1 = IA.retanalysisN(SPh,[-1, -1, -1, -1], [1, 1, 1, 1], 2, Ndiffs = 5, Nmax = 40)

    print('Multi Dimensional Integrals Now')

    print('SPh, N dimensional Integration')
    IA.plotme(RetAnN1,realvalue=4.9348022005)
    

    RetAnN2 = IA.retanalysisN(Circ, [-1, -1, -1, -1], [1, 1, 1, 1], 2, Nmax = 60, Ndiffs = 10)
     
    print('Circ, N dimensional Integration)')
    IA.plotme(RetAnN2, realvalue = 64/3, ylim = [-17.5, -1])

    print('Plotted together, N dimensional Integration')
    IA.PLOTME([RetAnN1, RetAnN2], LABEL=['SPh', 'Circ'], REALVALUE=[4.9348022005, 64/3], ylim = [-16, 2.5])

    RetAdtN1 = IA.retanalysisadaptN(SPh,[-1, -1, -1, -1], [1, 1, 1, 1], -4, intmethod = 2)

    print('SPh, Adaptive Quadrature N dimensions ')
    IA.plotme(RetAdtN1, realvalue = 4.9348022005, ylim = [-100, 100])

    RetAdtN2 = IA.retanalysisadaptN(Circ,[-1, -1, -1, -1], [1, 1, 1, 1],Taumin10to= -5, intmethod = 2)

    print('Circ, Adaptive Quadrature N dimensions N')
    IA.plotme(RetAdtN2, realvalue = 64/3, ylim = [-100, -29])

    print('Plotted Together, Adaptive Quadrature N dimensions')
    IA.PLOTME([RetAdtN1, RetAdtN2], LABEL = ['SPh', 'Circ'], REALVALUE = [4.9348022005, 64/3], ylim = [-75, -25] )

    RetMCN1 = IA.retanalysismonteN(SPh, [-1, -1, -1, -1], [1, 1, 1, 1], Nmax = 5000)

    print('SPh, Monte Carlo N dimensions')
    IA.plotme(RetMCN1, realvalue = 4.9348022005)

    RetMCN2 = IA.retanalysismonteN(Circ,[-1, -1, -1, -1], [1, 1, 1, 1], Nmax = 10000)

    print('Circ, Monte Carlo in N dimensions')
    IA.plotme(RetMCN2, realvalue = 64/3, ylim = [-6, 7.5])

    RetStratN1 = IA.retanalysisstratN(SPh, [-1, -1, -1, -1], [1, 1, 1, 1])

    print('SPh, Stratified Sampling N dimensions')
    IA.plotme(RetStratN1, realvalue=4.9348022005, ylim = [-0.3, 0])

    print('Circ, Stratified Sampling in N dimensions ')
    RetStratN2 = IA.retanalysisstratN(Circ, [-1, -1, -1, -1], [1, 1, 1, 1])
    IA.plotme(RetStratN2, realvalue = 64/3, ylim = [0, 0.25])
    IA.PLOTME([RetStratN1,RetStratN2], REALVALUE = [4.9348022005, 64/3], LABEL = ['SPh', 'Circ'], ylim = [-0.5, 0.5])

