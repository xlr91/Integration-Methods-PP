#import heartrate; heartrate.trace(browser=True)
import time
start_time = time.time()

def w(N, kind, h=1):
    
    if kind == 0:
        wlist = [h for i in range(N)]
        #check the https://www.value-at-risk.net/numerical-integration-multiple-dimensions/#exercise_2_25
        #cuz its [0, h]...

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

#GENERATOR WORKS
def rec(N, d, V, x = 0):
    if x != d-1:
        y = 1 * x
        x += 1    
        for i in range(N[y]):
            V[y] = i
            yield from rec(N, d, V, x)
    else:
        for i in range(N[x]):
            V[x] = i
            yield V
            #V[x] = i
            #can now stick a function in here related to V or make it a generator
            #print(V) #now i just need to do something with this
            #return V


def f(x, y, z):
    return x + y + z

def wlistval(V, Wlist):
    thevalue = 1
    for o in range(len(V)):
        thevalue *= Wlist[o][V[o]]
    return thevalue

def wfval(V, Wlist, f, X):
    thevalue = 1
    coords = []
    for o in range(len(V)):
        thevalue *= Wlist[o][V[o]]
        coords.append(X[o][V[o]])
    newtuple = tuple(coords)
    fname = f.__name__
    fstring = '{}{}'.format(fname, newtuple)
    fval = eval(fstring)
    return [thevalue, fval]

def fvalue(V, f, X):
    coords = []
    for o in range(len(V)):
        coords.append(X[o][V[o]])
    newtuple = tuple(coords)
    classname = f.__name__
    fstring = '{}{}'.format(classname, newtuple)
    return eval(fstring)




#Variables
A = [0, 0, 0]
B = [10, 10, 10]
N = [100, 100, 100]
kind = 2 #trapz 


d = len(N)#number of dimensions
V = [0 for i in range(d)]

dmul = 1
for i in range(d):
    dmul *= N[i]

#Add a check for length of A B and N maybe
value = 0


H = [(B[i]-A[i]) / N[i] for i in range(d)]

Wlist = [w(N[i], kind, H[i]) for i in range(d)]

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



gen = rec(N, d, V)


for _ in range(dmul):
    Vnow = next(gen)
    #value += wlistval(Vnow, Wlist)*fvalue(Vnow, f, X)
    
    valnow = wfval(V, Wlist, f, X)
    value += valnow[0]*valnow[1]
    

#V = [0, 0, 0]
#for all the Vs, its value += 
print(value)
time_taken = time.time()-start_time
print('timetaken', time_taken)
