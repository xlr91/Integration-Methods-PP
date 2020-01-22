def circle(x,y):
    """
    Returns 1 if x,y are inside a unit circle, 0 otherwise.
    Parameters:
    x : scalar
        x-coordinate value
    y : scalar
        y-coordinate value
    """
    value = 0
    if x**2 + y**2 <=1:
        value = 1
    return value

def f(x, y):
    value = -(x**2 + y**2) + 1
    if value < 0:
        value = 0
    return value

class Integrator:
    """
    A class used for monte calro integration
    ...

    Attributes
    ----------
    f : Function
        The function that will be integrated using monte carlo
    xmin, xmax : floats
        the minimum and maximum values of x, respectively
    ymin, ymax : floats
        the minimum and maximum values of y, respectively
    """
    def __init__(self, f, xmin, xmax, ymin, ymax):
        """
        Initializes the class
        ...

        Attributes
        ----------
        f : Function
            The function that will be integrated using monte carlo
        xmin, xmax : floats
            the minimum and maximum values of x, respectively
        ymin, ymax : floats
            the minimum and maximum values of y, respectively
        """
        self.f = f
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
    
    def mc(self, n = 1000):
        """
        Method that performs the Monte Carlo Integration
        ...

        Attributes
        ----------
        n : int
            Number of random points used
        """
        import random
        random.seed(1) #used for reproducibility

        f = self.f
        xmin = self.xmin
        xmax = self.xmax
        ymin = self.ymin
        ymax = self.ymax

        fvalue = 0
        
        for _ in range(n):
            x = random.uniform(xmin, xmax)
            y = random.uniform(ymin, ymax)
            #print(x,y)
            fvalue += f(x,y)
            #print(f(x,y))

        
        favg = fvalue/n
        print('favg', favg)
        return favg * (xmax-xmin) * (ymax-ymin)

TEST = Integrator(f, -1, 1, -1, 1)
TEST.mc(10000)