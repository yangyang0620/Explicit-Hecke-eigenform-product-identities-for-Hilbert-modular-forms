︠61c62132-0244-480f-9108-94a01f13b346︠
##############################################################################
# Title: Numerical Verification of Eigenform Product Identities (Eisenstein-Eisenstein Case with Equal Weight)
# Description:
#   Based on Propositions 3.8, this script verifies the non-existence of
#   eigenform identities for the equal-weight case.
# The script consists of two parts:
#   1. Determining the upper bound for the weight k.
#   2. Iterating through discriminants D to find satisfying fields and
#      verifying the identities.

##############################################################################
# The results of this program reproduce Table 1 of the paper, helping to verify the first part of Theorem 1.
##############################################################################



# ==========================================
# Helper Functions
# ==========================================


# Computes the rational value of the Dedekind Zeta function zeta_K(s).
# Rational values of the Dedekind zeta function derived from generalized Bernoulli numbers.
def Zeta(K,s):
    n = 1-s
    deg = K.degree()
    cond = K.conductor()
    G = DirichletGroup(cond)
    chi = G.0^(ZZ(G.order()/deg))
    Lval = 1
    for i in (1..deg-1):
        Lval *= (-(chi^i).bernoulli(n)/n)
        return zeta(s) * Lval



# Verifies the Constant Term Identity (Eq 3.1):
#1 = | (zeta_F(1-k1) + zeta_F(1-k2)) * zeta_F(1-k1-k2) / (zeta_F(1-k1)*zeta_F(1-k2)) |, here k1=k2
def check_equation(a, b, c):
    return c * (a + b) == a * b



# Computes the boundary function C(D, k).
# Used to determine if eigenform identities can no longer exist for D > current.
def C(D, k):
    term1 = ((108/pi^6)^2 * D^(1/2) * k).n()
    return term1



# Filter potential quadratic fields with narrow class number 1 using the prime factorization of the discriminant
# Conditions:
#    1. Odd prime factors p must satisfy p = 1 mod 4 and have exponent 1.
#    2. The exponent of factor 2 must be 0 (n is odd) or 3 (n contains factor 8).

def check_field(n):
    if n <= 0:
        return False
    factors = factor(n)
    power_of_2 = 0
    for (p, e) in factors:
        if p == 2:
            power_of_2 = e
        else:
            if e != 1:
                return False
            if p % 4 != 1:
                return False
    if power_of_2 == 0 or power_of_2 == 3:
        return True
    return False



# ==========================================
# Main Program Logic
# ==========================================



# ------------------------------------------
# 1. Display Symbolic Expression
# ------------------------------------------

var('D k ')
print("="*100)
print("Symbolic Expression of C(D, k):")
expression_C = ((pi/2) * (36/pi^4)^3 * (D/(4 * pi^2))^(1/2) * k)
show(expression_C)
print("="*100)


# ------------------------------------------
# 2. Determine the Upper Bound for Weight Search
# We fix the minimum D=13 and increase k to find the critical point where CC > 1.
# ------------------------------------------

Maxnumber_k = 100   # Maximum search range for weights
Maxnumber_D = 2000  # Maximum search range for discriminant
k_start = 2 # Starting value for k (even integer)
D_start=5  # Starting value for D
k_max=0
found_k = False
found_D = False
print("\nDetermining the upper bound for the weight...")

for i in range(k_start, Maxnumber_k, 2):
    result = C(D_start,i).n()
    if result > 1:
        k_max = i - 2
        print(f"For the equal-weight case of two Eisenstein series, if the discriminant >={D_start}, the maximum possible valid weight is {k_max}.")
        found_k = True
        break
if not found_k:
    print(f"Warning: Cutoff for k not found within Maxnumber={Maxnumber_k}. Consider increasing range.")
print("="*100)


# ------------------------------------------
# 3. Determine the maximal discriminant for fixed k
# ------------------------------------------

if found_k:
    print("\nStarting traversal of weights and discriminants...")
    for i in range(2, k_max+1, 2):   # Iterate through every even weight k
        for D in range(D_start, Maxnumber_D, 1):  # Iterate through discriminants D
            if not check_field(D):
                continue     # Skip this iteration if D fails the necessary condition for a narrow class number of 1
            result = C(D,i).n()
            if result > 1:
                # === Boundary Triggered: Backtracking ===
                # The current D is too large. We need to backtrack to find the
                # largest D that satisfies specific conditions.
                for DD in range(D-1, D_start-1, -1):
                    if not check_field(DD):
                        continue
                    # by Prop.3.8, we only need to consider the case that (2) is inert, which means                                     # D\equiv 3 or -3 \mod 8.
                    if DD % 8 in [3, 5]:
                        K.<a> = QuadraticField(DD)
                        h = K.narrow_class_group().order()
                        if h == 1:
                            print(f"For weight {i}, the maximum valid D is {DD}.")
                            found_D = True
                            break
                break
            K.<a> = QuadraticField(D)
            aa = Zeta(K, 1-i).n()
            bb = Zeta(K, 1-i).n()
            cc = Zeta(K, 1-2*i).n()
            if check_equation(aa, bb, cc):
                print(f"For D={D}, k={i},the eigenform product identity hold.")
        if not found_D:
            print(f"Warning: Cutoff for D not found within Maxnumber={Maxnumber_D}. Consider increasing range.")
            break
        print(f"For k={i},the list of eigenform product identities is completed.")
        print("-"*60)
print(f"{'='*100}")
print("\nVerification process completed.")
︡40431c96-c497-48de-9930-6e288f4c8686︡{"stdout":"(D, k)\n"}︡{"stdout":"====================================================================================================\n"}︡{"stdout":"Symbolic Expression of C(D, k):\n"}︡{"html":"<div align='center'>$\\displaystyle \\frac{11664 \\, \\sqrt{D} k}{\\pi^{12}}$</div>"}︡{"stdout":"====================================================================================================\n"}︡{"stdout":"\nDetermining the upper bound for the weight...\n"}︡{"stdout":"For the equal-weight case of two Eisenstein series, if the discriminant >=5, the maximum possible valid weight is 34.\n"}︡{"stdout":"====================================================================================================\n"}︡{"stdout":"\nStarting traversal of weights and discriminants...\nFor D=5, k=2,the eigenform product identity hold.\nFor weight 2, the maximum valid D is 1549."}︡{"stdout":"\nFor k=2,the list of eigenform product identities is completed.\n------------------------------------------------------------\nFor weight 4, the maximum valid D is 389."}︡{"stdout":"\nFor k=4,the list of eigenform product identities is completed.\n------------------------------------------------------------\nFor weight 6, the maximum valid D is 173."}︡{"stdout":"\nFor k=6,the list of eigenform product identities is completed.\n------------------------------------------------------------\nFor weight 8, the maximum valid D is 61."}︡{"stdout":"\nFor k=8,the list of eigenform product identities is completed.\n------------------------------------------------------------\nFor weight 10, the maximum valid D is 61."}︡{"stdout":"\nFor k=10,the list of eigenform product identities is completed.\n------------------------------------------------------------\nFor weight 12, the maximum valid D is 37."}︡{"stdout":"\nFor k=12,the list of eigenform product identities is completed.\n------------------------------------------------------------\nFor weight 14, the maximum valid D is 29."}︡{"stdout":"\nFor k=14,the list of eigenform product identities is completed.\n------------------------------------------------------------\nFor weight 16, the maximum valid D is 13."}︡{"stdout":"\nFor k=16,the list of eigenform product identities is completed.\n------------------------------------------------------------\nFor weight 18, the maximum valid D is 13."}︡{"stdout":"\nFor k=18,the list of eigenform product identities is completed.\n------------------------------------------------------------\nFor weight 20, the maximum valid D is 13."}︡{"stdout":"\nFor k=20,the list of eigenform product identities is completed.\n------------------------------------------------------------\nFor weight 22, the maximum valid D is 5."}︡{"stdout":"\nFor k=22,the list of eigenform product identities is completed.\n------------------------------------------------------------\nFor weight 24, the maximum valid D is 5."}︡{"stdout":"\nFor k=24,the list of eigenform product identities is completed.\n------------------------------------------------------------\nFor weight 26, the maximum valid D is 5."}︡{"stdout":"\nFor k=26,the list of eigenform product identities is completed.\n------------------------------------------------------------\nFor weight 28, the maximum valid D is 5."}︡{"stdout":"\nFor k=28,the list of eigenform product identities is completed.\n------------------------------------------------------------\nFor weight 30, the maximum valid D is 5."}︡{"stdout":"\nFor k=30,the list of eigenform product identities is completed.\n------------------------------------------------------------\nFor weight 32, the maximum valid D is 5.\nFor k=32,the list of eigenform product identities is completed.\n------------------------------------------------------------\nFor weight 34, the maximum valid D is 5."}︡{"stdout":"\nFor k=34,the list of eigenform product identities is completed.\n------------------------------------------------------------\n"}︡{"stdout":"====================================================================================================\n"}︡{"stdout":"\nVerification process completed.\n"}︡{"done":true}









