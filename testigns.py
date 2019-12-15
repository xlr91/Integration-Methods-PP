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

'''
#need to be able to generalize the N for loops
for i in range(N1):
    for j in range(N2):
        for k in range(N3)
            have Wlist[0][i]*Wlist[1][j]*Wlist[2][k]*f(xi, xj, xk)
'''