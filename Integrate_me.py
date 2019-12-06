def f(x):
    
    return x**4



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
        for i in range(N):
            xa = x[i]
            xb = x[i+1]
            value += (xb - xa) * f((xb+xa)/2)
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

        for i in range(1, N-1):
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


    def Integrate(self, N, kind):
        value = 0
        a = self.xmin
        b = self.xmax
        h = (b-a) / N
        wlist = self.w(N, kind, h)
        x = [a + i*h for i in range(N+1)]
        #wlist = self.w(N, kind, h)
        for i in range(N+1):
            #print(i)
            #xi = self.xmin + i*h
            #print('coord', xi)
            value += wlist[i] * self.f(x[i])
            #print(value)
        return value


    def w(self, N, kind, h=1):
        'returns w list of values'
        if kind == 0:
            wlist = [2 for i in range(N)]

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
    




intme = Integrator(f, 0, 10)
ltest1 = intme.w(10, 1)
ltest2 = intme.w(10, 2)


