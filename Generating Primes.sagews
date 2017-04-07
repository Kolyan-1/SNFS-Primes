import numpy
from time import time
####################
# Helper Functions
####################
#reference sagemath
def random_between(j,k):
    a=int(random()*(k-j+1))+j
    return a
#get prime efficiently (tested)
def get_prime(bits_q):
    q = random_prime(2^bits_q-1,False,2^(bits_q-1))
    while(not is_prime(q)):
        q = random_prime(2^bits_q-1,False,2^(bits_q-1))
    return q

#get f poly
def get_f(X,degree_f, bits_q):
    flag_irreducible = False
    while (not flag_irreducible):
        #f_vec = [ZZ.random_element(-int(2^(10)-1),int(2^(10)-1)) for _ in range(degree_f+1)]
        f_vec = [ZZ.random_element(-int(2^(bits_q/(2*(degree_f+1)))),int(2^(bits_q/(2*(degree_f+1))))) for _ in range(degree_f+1)]
        #f_vec = [ZZ.random_element(1,int(2^(bits_q/(2*(degree_f+1))))) for _ in range(degree_f+1)]

        norm_f = max(map(abs,f_vec))
        f_poly = X(list(f_vec))
        if f_poly.is_irreducible():
            flag_irreducible = True
    return f_poly,norm_f
#get g poly
def get_g(g1,x,bits_p,degree_f,norm_f):
    g0 = ZZ.random_element(-int(2^(bits_p/degree_f)/norm_f),int(2^(bits_p/degree_f)/norm_f))
    #g0 = ZZ.random_element(1,int(2^(bits_p/degree_f)/norm_f))
    g_poly = g1*x+g0
    return g_poly,g0
def get_G(f_poly,g_poly):
    G_temp = f_poly.sylvester_matrix(g_poly,variable=x)
    G_poly = G_temp.determinant()-1
    return G_poly

def get_their_poly(X):
    f_coefs =  list(reversed([1155,1090,440,531,-348,-223,-1385]))
    f_poly  = X(f_coefs)
    norm_f = 1385
    return f_poly,norm_f

    
######################################################################################
######################################################################################
######################################################################################
############           main script
######################################################################################
######################################################################################
######################################################################################

### Parameters ###
degree_f = 8
bits_q = 63 #
bits_p = bits_q+1 #2q+1

## Other Variables ##
count_primes = 0
NumIt = 500

### Main Execution
t0 = time()


# Generating required Rings
X.<x>     = ZZ['x']
G.<x,g1>  = ZZ['x,g1']
# Main Loop
for i in range(NumIt):

    #generating prime and Associated Field
    q = get_prime(bits_q)
    T.<g2>    = Integers(q)['g2']

    flag_roots = True
    #getting G with roots
    while flag_roots:
        f_poly,norm_f = get_f(X,degree_f,bits_q)
        #print(f_poly)
        #f_poly,norm_f = get_their_poly(X)
        g_poly,g0 = get_g(g1,x,bits_p,degree_f, norm_f)

        G_poly = get_G(f_poly,g_poly)
        temp = list(G_poly.coefficients())
        temp.reverse()
        T2 = T(temp)
        if len(T2.roots())>0:
            flag_roots=False
    # getting roots
    r    = T2.roots()
    root = r[0][0]
    rt   = int(root)
#    while(rt<int(2^(bits_p/degree_f)/norm_f)):
#        rt+=q
    p    = X([G_poly(1,rt)+1])


    if p<0:
        p = -p
        dummy123=1
        #print('swapping signs')
    else:
        dummy123=0#print('Not swapping signs')
    if is_prime(p):
        print('p is prime with bits:'+str(ceil(log(p,base=2))))
        if(mod(p-1,q)==0):
            count_primes+=1
            print(str(p)+','+str(degree_f)+','+str(bits_q)+','+str(ceil(log(p,base=2))))
        else:
            print(dummy123)
            print('But q does not divide p-1')

    else:
        pass
        #print('Iteration:%d'%(i))
print('Ratio of Primes Found: %f'%(count_primes/NumIt))
print(time()-t0)
