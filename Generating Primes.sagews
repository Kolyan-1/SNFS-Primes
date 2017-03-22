︠6d69cd82-4def-4784-bab8-85e896cf982as︠
import numpy
from time import time
####################
# Helper Functions
####################

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
        f_vec = vector(numpy.random.randint(1,high=2^(bits_q/(2*(degree_f+1))), size = (1,degree_f))[0])
        norm_f = max(f_vec)
        f_poly = X(list(f_vec))
        if f_poly.is_irreducible():
            flag_irreducible = True
    return f_poly,norm_f
#get g poly
def get_g(g1,x,bits_p,degree_f,norm_f):
    g0 = numpy.random.randint(0,high=(2^(bits_p/degree_f))/norm_f)
    g_poly = g1*x+g0
    return g_poly,g0
def get_G(f_poly,g_poly):
    G_temp = f_poly.sylvester_matrix(g_poly,variable=x)
    G_poly = G_temp.determinant()-1
    return G_poly

#################
#main script 
#################

### Parameters ###
degree_f = 6
bits_q   = 64
bits_p = bits_q #not testing for size p

## Other Variables ##
count_primes = 0
NumIt = 1000

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
        g_poly,g0 = get_g(g1,x,bits_p,degree_f, norm_f)

        G_poly = get_G(f_poly,g_poly)
        T2 = T(G_poly.coefficients())
        if len(T2.roots())>0:
            flag_roots=False
    # getting roots
    r    = T2.roots()
    root = r[0][0]
    rt         = X([root])
    new_g_poly = rt*x+g0
    p       = X([G_poly(1,rt)+1])

    if p<0:
        p = -p
    if is_prime(p):
        print(str(p)+','+str(degree_f)+','+str(bits_q)+','+str(ceil(log(p,base=2))))
        
        count_primes+=1
    else:
        pass
        #print('Iteration:%d'%(i))
print('Ratio of Primes Found: %f'%(count_primes/NumIt))
print(time()-t0)
︡4fbbeef1-a6f8-4a3f-9029-37ecd2d29788︡{"stdout":"92564209958231537167650984747444230876650374293570039178161348230228667530588085631159661183861,6,64,316"}︡{"stdout":"\n5897194432735738704530328586531397505223439258252493513647760457731692166569030649566678341837,6,64,312"}︡{"stdout":"\n3481429121996804280828326664255515121400745963831888213702327587228747506074067702106031531819573,6,64,321"}︡{"stdout":"\n6882178034041956672702582599873164541454547637560448740559566193958816703768705767430385077927,6,64,312"}︡{"stdout":"\n3890996701491561690919889212762632186725526951862171318312739291225982434491602464840384909339,6,64,311"}︡{"stdout":"\n"}︡{"stdout":"Ratio of Primes Found: 0.005000\n"}︡{"stdout":"33.9378631115\n"}︡{"done":true}︡









