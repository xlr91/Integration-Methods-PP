

## check what the uncertainty is as a function of sampled points n
# demonstrate convergence, perform timing tests

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
    def __init__(self, f):
        '''
        Initializer of the function. f is the function meant to be integrated. 
        '''
        self.f = f
        self.i = 0
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

        xin = [] #unused
        for i in range(N):
            xa = x[i]
            xb = x[i+1]
            value += (xb - xa) * f((xb+xa)/2)
            xin.append((xb+xa)/2)
        
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

    def AdaptInt1(self, a, b, tau, intmeth, Q0 = 10000):
        '''
        Integrate using the Adaptive Integration Method, using the Newton-Cotes Rule. 
        Uses the Big division - Small division to determine the error

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
        self.i += 1

        val = self.NCInt( a, b, N, intmeth) 
        Q1 = +val
        err = abs(val - Q0)

        if err > tau:
            m = (a+b)/2
            val = self.AdaptInt1(a, m, tau, intmeth, Q1) + self.AdaptInt1(m, b, tau, intmeth, Q1)
        
        return val
      
    def AdaptInt2(self, a, b, tau, intmeth):
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

        err = val1 - val2
        if abs(err) > tau:
            if err > 0:
                val = self.AdaptInt2(a, m, tau, intmeth)
                val += self.AdaptInt2(m, b, tau, intmeth)
            else:
                val = self.AdaptInt2(m, b, tau, intmeth)
                val += self.AdaptInt2(a, m, tau, intmeth)
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
        return favg * (a-b)

        
    def NCIntN(self, A, B, N, kind):
        #A = [a1, a2, a3, ...] the lower bounds in each dimension
        #B = [b1 ,b2, b3, ...] the upper bounds in each dimension
        #N = [N1, N2, N3, ...] the number of points in each dimension
        #all must be of the same length
        A = [0, 0]
        B = [10, 10]
        N = [10, 20]
        
        d = len(N)#number of dimensions
        #Add a check for length of A B and N maybe
        value = 0

        H = [(B[i]-A[i]) / N[i] for i in range(d)]
        #h = (b-a) / N
        Wlist = [self.w(N[i], kind, H[i]) for i in range(d)]
        #wlist = self.w(N, kind, h)
        X = [[A[i] + j*H[i] for j in range(N[i]+1)] for i in range(d)]
        

        if kind == 0: #deals with the cmidpoint rule
            for i in range(d):
                x = X[i]
                xnew = []
                for j in range(N[i]):
                    xa = x[j]
                    xb = x[j+1]
                    xnew.append((xb+xa)/2)
                X[i] = xnew
                N[i] -= 1
        #wlist = self.w(N, kind, h)
        
    
        #need to work on this part now
        for i in range(N+1):
            #print(i)
            #xi = self.xmin + i*h
            #print('coord', xi)
            value += Wlist[i] * self.f(x[i])
            #print(value)
        #return value
        return value




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



def functiontester(xmin, xmax, VAL):
    flist = [f1, f2, f3, f4, f5]
    for func in flist:
        intme = Integrator(func)
        for i in range(3):
            ind = flist.index(func)
            Rval = VAL[ind]
            intme.plotmeval(xmin, xmax, i, realvalue = Rval )

if __name__ == '__main__':

    a = 0 
    b = 10
    VAL = [50, 1000/3, 2500, 20000, 500000/3]

    functiontester(a, b, VAL)

    
    
    intme = Integrator(f1)
    ltest1 = intme.w(10, 1, 1)
    ltest2 = intme.w(10, 2, 1)
    import time

    print('small - small')
    start_time = time.time()
    print(intme.AdaptInt2(0, 10, 1e-5, 2))
    print(time.time()-start_time)
    print('calls', intme.j)
    print('')
    print('large - small')
    start_time = time.time()
    print('value', intme.AdaptInt1(0, 10, 1e-2, 2))
    print('time taken', time.time()-start_time)
    print('calls', intme.i)



    A = [0, 0]
    B = [10, 10]
    N = [10, 20]
    kind = 1
    
    d = len(N)#number of dimensions
    #Add a check for length of A B and N maybe
    value = 0

    H = [(B[i]-A[i]) / N[i] for i in range(d)]
    #h = (b-a) / N
    Wlist = [intme.w(N[i], kind, H[i]) for i in range(d)]
    #wlist = self.w(N, kind, h)
    X = [[A[i] + j*H[i] for j in range(N[i]+1)] for i in range(d)]
    

    if kind == 0: #deals with the cmidpoint rule
        for i in range(d):
            x = X[i]
            xnew = []
            for j in range(N[i]):
                xa = x[j]
                xb = x[j+1]
                xnew.append((xb+xa)/2)
            X[i] = xnew
            N[i] -= 1

    #need to be able to generalize the N for loops
    '''
    for i in range(N1):
        for j in range(N2):
            for k in range(N3)
                have Wlist[0][i]*Wlist[1][j]*Wlist[2][k]*f(xi, xj, xk)
    '''
    #for i in range(3):
        #intme.plotmeval(i, realvalue = 20000, Nmax = 25000, Ndiffs = 100)
    #    intme.plotmeval(i, realvalue = 20000)
    #intme.plotme(Nmax = 500, realvalue = 1000/3)
    
