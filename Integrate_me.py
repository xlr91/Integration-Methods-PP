

## check what the uncertainty is as a function of sampled points n
# demonstrate convergence, perform timing tests

class Integrator:
    def __init__(self, f, xmin, xmax):
        self.f = f
        self.xmin = xmin
        self.xmax = xmax
    
    def __repr__(self):
        newtuple = tuple([self.f, self.xmax, self.xmin])
        classname = self.__class__.__name__
        return '{}{}'.format(classname, newtuple)

    def midpoint(self):
        a = self.xmin 
        b = self.xmax
        R = b-a
        return R*self.f((b+a)/2)

    def cmidpoint(self, N):
        a = self.xmin 
        b = self.xmax
        f = self.f
        h = (b-a)/N
        x = [a + i*h for i in range(N+1)]
        value =  0
        xin = []
        for i in range(N):
            xa = x[i]
            xb = x[i+1]
            value += (xb - xa) * f((xb+xa)/2)
            xin.append((xb+xa)/2)
        
        #return value
        return value

    def trapz(self):
        a = self.xmin 
        b = self.xmax
        f = self.f
        R = (b-a)/2
        return R*(f(a) + f(b))

    def ctrapz(self, N):
        a = self.xmin 
        b = self.xmax
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

    def simpsons(self):
        a = self.xmin 
        b = self.xmax
        f = self.f
        R = (b-a)/6
        return R * (f(a) + f(b) + 4*f((a+b)/2))
    
    def csimpsons(self, N):
        a = self.xmin 
        b = self.xmax
        f = self.f
        h = (b-a)/N
        R = h/3
        x = [a + i*h for i in range(N+1)]
        thelist = []
        for i in range(N):
            if i % 2 == 0:
                thelist.append(2)
            else:
                thelist.append(4)
        thelist[0] = 1
        thelist.append(1)
        #print(thelist)
        value = 0 
        for i in range(N+1):
            value += R * thelist[i] * f(x[i])
        return value


    def retanalysis(self, Nmax, Ndiffs, intmethod):
        import time
        
        nPoints = [2]
        start_time = time.time()
        intvalue = [self.NCInt(nPoints[0], intmethod)]
        time_taken = [time.time()-start_time]

        for i in range(int(Nmax/Ndiffs)):
            npoint = (i+1)*Ndiffs
            
            start_time = time.time()
            tintvalue = self.NCInt(npoint, intmethod)
            timepast = time.time() - start_time

            nPoints.append(npoint)
            intvalue.append(tintvalue)
            time_taken.append(timepast)

        return [nPoints, intvalue, time_taken]

    def plotme(self, intmethod, Nmax = 500, Ndiffs = 10, realvalue = None):
        
        import matplotlib.pyplot as plt 
        
        RetAn = self.retanalysis(Nmax, Ndiffs, intmethod)

        plt.figure(0)
        plt.title('Value of Integral as a function of sample points')
        #plt.xlim(0, xlim)
        #plt.ylim(0, 100)
        #plt.hlines(10, 0, self.z[len(self.z)-1])
        plt.plot(RetAn[0], RetAn[1])
        
        #plt.xlabel('Redshift (z)')
        #plt.ylabel('Percentage Difference ')
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

    def plotmeval(self, intmethod, Nmax = 500, Ndiffs = 10, realvalue = None):
        import matplotlib.pyplot as plt
        import numpy as np
        import matplotlib.gridspec as gridspec

        RetAn = self.retanalysis(Nmax, Ndiffs, intmethod)

        fig = plt.figure(tight_layout=True)
        #fig.suptitle('sup bruh ')
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

        fig.align_labels()  # same as fig.align_xlabels(); fig.align_ylabels()

        plt.show()

    def NCInt(self, N, kind):
        #does it well for simpsons and trapz
        #works for cmidpoint, more testing is required
        #THis is newton-cotes integration, basically its all of the 
        #three basic integrations all in one method (!!!)
        value = 0
        a = self.xmin
        b = self.xmax
        h = (b-a) / N
        wlist = self.w(N, kind, h)
        x = [a + i*h for i in range(N+1)]

        if kind == 0: #deals with the cmidpoint rule
            xnew = []
            for i in range(N):
                xa = x[i]
                xb = x[i+1]
                xnew.append((xb+xa)/2)
            x = xnew
            N -= 1
        #wlist = self.w(N, kind, h)
        for i in range(N+1):
            #print(i)
            #xi = self.xmin + i*h
            #print('coord', xi)
            value += wlist[i] * self.f(x[i])
            #print(value)
        #return value
        return value

    def w(self, N, kind, h=1):
        
        if kind == 0:
            wlist = [h for i in range(N)]

        elif kind == 1:
            wlist = [h for i in range(N)]
            wlist[0] = h/2
            wlist.append(h/2)

        elif kind == 2:
            wlist = []
            for i in range(N):
                if i % 2 == 0:
                    wlist.append(2 * h / 3)
                else:
                    wlist.append(4 * h / 3)
            wlist[0] = 1 * h/3
            wlist.append(1 * h/3)

        return wlist

    def AdaptInt(self, a, b, tau, intmeth, Q0 = 1000):
        
        N = 1
        if intmeth == 2:
            N = 2

        val = self.NCInt(N, intmeth) #this needs to be fixed because it needs to take range
        err = abs(val - Q0)
        if err >= tau:
            m = (a+b)/2
            val = self.AdaptInt(a, m, tau, intmeth, Q0)
            val += self.AdaptInt(m, b, tau, intmeth, Q0)
        return val
    



        #recursion bois
        
    




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
        intme = Integrator(func, xmin, xmax)
        for i in range(3):
            ind = flist.index(func)
            Rval = VAL[ind]
            intme.plotmeval(i, realvalue = Rval )

if __name__ == '__main__':

    a = 0 
    b = 10
    VAL = [50, 1000/3, 2500, 20000, 500000/3]

    #functiontester(a, b, VAL)

    
    
    intme = Integrator(f3, 0, 10)
    ltest1 = intme.w(10, 1)
    ltest2 = intme.w(10, 2)

    #for i in range(3):
        #intme.plotmeval(i, realvalue = 20000, Nmax = 25000, Ndiffs = 100)
    #    intme.plotmeval(i, realvalue = 20000)
    #intme.plotme(Nmax = 500, realvalue = 1000/3)
    
