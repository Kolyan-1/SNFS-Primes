import numpy
from time import time
####################
# Helper Functions
####################
#reference sagemath
def print_times(arr):
    print('0! ' + str(arr[0]))
    print('1! ' + str(arr[1]))
    print('2! ' + str(arr[2]))
    print('3! ' + str(arr[3]))
    print('4A ' + str(arr[4]))
    print('4! ' + str(arr[5]))


def test_p(prime_p,G_poly,root):
    P.<y> = Integers(prime_p)['y']
    temp = list(G_poly.coefficients())
    temp.reverse()
    P2    = P(temp)
    rt1    = P(root)
    print(P2(rt1))
    if(P2(rt1)==prime_p-1):
        print('OK')
    else:
        print('NOT OK!')
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
        #f_vec = [ZZ.random_element(-int(2^(bits_q/(2*(degree_f+1)))),int(2^(bits_q/(2*(degree_f+1))))) for _ in range(degree_f+1)]
        f_vec = [ZZ.random_element(1,int(2^(bits_q/(2*(degree_f+1))))) for _ in range(degree_f+1)]

        norm_f = max(map(abs,f_vec))
        f_poly = X(list(f_vec))
        if f_poly.is_irreducible():
            flag_irreducible = True
    return f_poly,norm_f
#get g poly
def get_g(g1,x,bits_p,degree_f,norm_f):
    #g0 = ZZ.random_element(-int(2^(bits_p/degree_f)/norm_f),int(2^(bits_p/degree_f)/norm_f))
    g0 = ZZ.random_element(1,int(2^(bits_p/degree_f)/norm_f))
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
#################
#main script
#################

### Parameters ###
degree_f = 6
bits_q = 255 #
bits_p = bits_q+1 #2q+1

## Other Variables ##
NumIt = 1000


###### BENCHMARKING VECTORS
times        = numpy.zeros([NumIt,6])
count_primes = numpy.zeros(2)
bits_ps      = 0

### Main Execution
t0 = time()

t1 = time()
# Generating required Rings
X.<x>     = ZZ['x']
G.<x,g1>  = ZZ['x,g1']
#print('0! ' + str(time()-t1))
times[0,0] = time()-t1

# Main Loop
for i in range(NumIt):
    t1 = time()
    #generating prime and Associated Field
    q = get_prime(bits_q)
    T.<g2>    = Integers(q)['g2']
    #print('1! ' + str(time()-t1))
    times[i,1] += time()-t1


    flag_roots = True
    #getting G with roots
    t1 = time()
    while flag_roots:
        t2 = time()
        t3 = time()
        f_poly,norm_f = get_f(X,degree_f,bits_q)
        #print('2F ' + str(time()-t3))

        #print(f_poly)
        #f_poly,norm_f = get_their_poly(X)
        t3 = time()
        g_poly,g0 = get_g(g1,x,bits_p,degree_f, norm_f)
        #print('2g ' + str(time()-t3))

        t3 = time()
        G_poly = get_G(f_poly,g_poly)
        #print(G_poly)
        #print('2G ' + str(time()-t3))

        t3 = time()
        temp = list(G_poly.coefficients())
        temp.reverse()
        T2 = T(temp)
        #print('2Q ' + str(time()-t3))

        #print('2A ' + str(time()-t2))

        t2 = time()
        if len(T2.roots())>0:
            flag_roots=False
        #print('2B ' + str(time()-t2))
    #print('2! ' + str(time()-t1))
    times[i,2] += time()-t1





    t1 = time()
    # getting roots
    r    = T2.roots()
    root = r[0][0]
    rt   = int(root)
#    while(rt<int(2^(bits_p/degree_f)/norm_f)):
#        rt+=q
    p    = X([G_poly(1,rt)+1])
    #print('3! ' + str(time()-t1))
    times[i,3] += time()-t1


    if p<0:
        p = -p
        dummy123=1
        #print('swapping signs')
    else:
        dummy123=0#print('Not swapping signs')
    t1=time()
    if is_prime(p):
        #print('p is prime with bits:'+str(ceil(log(p,base=2))))
        #test_p(p,G_poly,rt)
        bits_ps += ceil(log(p,base=2))
        t2 = time()
        count_primes[0]+=1

        if(mod(p-1,q)==0):
            count_primes[1]+=1
            times[i,4] += time()-t2
            #print(str(p)+','+str(degree_f)+','+str(bits_q)+','+str(ceil(log(p,base=2))))
            #print('4A ' + str(time()-t2))

        else:
            #print(dummy123)
            #print('But q does not divide p-1')
            dummy123+=2
        #print('4! ' + str(time()-t1))
        times[i,5] += time()-t1
    else:
        pass
        #print('Iteration:%d'%(i))
print('Total Time: '+ str(time()-t0))
times = numpy.mean(times,axis=0) #correct axis already
bits_p_final = bits_ps/count_primes[0]
print_times(times)
print('Ratio of Primes Found: %f'%(count_primes[0]/NumIt))
print('Average Bit size p:' + str(bits_p_final))
