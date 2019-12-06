#Midpoint rule
def f(x):
    
    return -x**2 + 1



class Midpoint:
    def __init__(self, f, xmin, xmax):
        self.f = f
        self.xmin = xmin
        self.xmax = xmax
    def midpoint(self):
        a = self.xmin 
        b = self.xmax
        R = b-a
        return R*self.f((b+a)/2)
    



intme = Midpoint(f, -1, 1)