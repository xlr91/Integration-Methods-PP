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
        
    def AdaptInt2N(self, A, B, tau, intmeth):

        N = 1
        if intmeth == 2:
            N = 2
        self.j += 1

        val = self.NCIntN(A, B, N, intmeth) #this needs to be fixed because it needs to take range
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