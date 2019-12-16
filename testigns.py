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


A = [0, 0]
B = [10, 10]
N = [10, 20]
kind = 1

d = len(N)#number of dimensions
#Add a check for length of A B and N maybe
value = 0

H = [(B[i]-A[i]) / N[i] for i in range(d)]
#h = (b-a) / N
Wlist = [w(N[i], kind, H[i]) for i in range(d)]
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


def loop_rec(y, n):
    if n >= 1:
        for _ in range(y):
            loop_rec(y, n - 1)
    else:
       print(Wlist[y][n])


def loopme(I, n, wlistval):
    wlistval 
        



#maybe we can start with value
N = [2, 3, 4]
d = 3
V = [0, 0, 0]
for x in range(3):
    V[0] = x
    for y in range(3):
        V[1] = y
        for z in range(3):
            V[2] = z
            #print([x, y, z])
            print(V)
            #V = []

N = [2, 3, 4]
d = 3
V = [0, 0, 0]

def rec(d, V, n = 0, x = None):
    d -= 1
    if d >= 1:
        for x in range(N[d]):
            V[d] = x
            #print('x', x)
            rec(d, V, x)
    else:
        for x in range(N[d]):
            V[d] = x
            for i in range(len(N)):
                print(V)
                thevalue *= Wlist[i][V[i]]
    
'''
#thevalue = f(xi, xj, xk), say its 1 for ease
thevalue = 1 
so then what you wanna do is 
thevalue *= Wlist[0][i]*Wlist[1][j]*Wlist[2][k]*...

so its 
thevalue *= Wlist[0][i]
thevalue *= Wlist[1][j]
thevalue *= Wlist[2][k]

so just for that i can have a 


ijk = [0, 0, 0] #this can be easily made
#just need a way to change them all 
thevalue = 1
for o in range(d):
    thevalue *= Wlist[o][ijk[o]]
value += thevalue

ijk[0] +=1
if ijk[0] == something[0]:
    ijk[0] = 0
    ijk[0+1] += 1

now i need to find a way to loop over the correct ijks  
#need to be able to generalize the N for loops
for i in range(N1):
    for j in range(N2):
        for k in range(N3)
            have Wlist[0][i]*Wlist[1][j]*Wlist[2][k]*f(xi, xj, xk)
'''