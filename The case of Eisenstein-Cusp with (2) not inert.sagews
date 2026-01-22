︠c12b2183-c30a-4bb9-8aef-ff2e3746db54s︠
################################################################################
# Title: Numerical Verification of Eigenform Product Identities (Eisenstein-Cusp Case with (2) Not Inert)
# Description:
#   Based on Propositions 4.13, this script verifies the non-existence of
#   eigenform identities with (2) not inert.
# The script consists of two parts:
#   1. Determining the upper bound for the weight k1.
#   2. Iterating Iterate through each fixed k1 to find the upper bound for k2, and for each k2, find the upper bound for D
################################################################################


##############################################################################
# The results of this program reproduce Table 3 of the paper, helping to verify the second part of Theorem 1.
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


# These functions are derived from coefficient estimates.
# Formula (4.11)
def CC(k1):
    term1 = pi^5/18*(2*pi)^(2*(k1-1/2))*gamma(k1)^(-2)
    term2=2^(2*k1/2)-1-term1
    return term2.n()


# Formula (4.10)
def DD(k1,k2):
    term1 = pi^5/18*(2*pi)^(2*(k1-1/2))*gamma(k1)^(-2)
    term2=2^(2*(k2-1)/2)*(2^(2 * k1/2) - 1)-term1
    return term2.n()


# Formula (4.9)
def EE(k1,k2,D):
    term1 = pi^5/18*(2*pi)^(2*(k1-1/2))*gamma(k1)^(-2)/(D^(k1-1/2))
    term2=2^(2*(k2-1)/2)*(2^(2 * k1/2) - 1)-term1
    return term2.n()


# ==========================================
#  Display Symbolic Expression
# ==========================================


print(f"{'='*100}")
print("We begin dealing with the question of whether the eigenform product identity g=fh exists. In this case, we require one of f or h to be a cusp form and the other an Eisenstein series, and the ideal (2) to be not inert.")
print("="*100)

var('k1,k2,D')


term1=2^(2*k1/2)-1-pi^5/18*(2*pi)^(2*(k1-1/2))*gamma(k1)^(-2)
print("\nSymbolic Expression of CC(k1):\n")
show(term1)


term1 = pi^5/18*(2*pi)^(2*(k1-1/2))*gamma(k1)^(-2)
term2=2^(2*(k2-1)/2)*(2^(2 * k1/2) - 1)-term1
print("\nSymbolic Expression of DD(k1,k2):\n")
show(term2)

term1 = pi^5/18*(2*pi)^(2*(k1-1/2))*gamma(k1)^(-2)/(D^(k1-1/2))
term2=2^(2*(k2-1)/2)*(2^(2 * k1/2) - 1)-term1
print("\nSymbolic Expression of EE(k1,k2,D)\n")
show(term2)



# ==========================================
# Main Program Logic
# ==========================================


# Settings
Maxnumber_k1 = 100   # Maximum search range for weight k1
Maxnumber_k2 = 1000  # Maximum search range for weight k2
Maxnumber_D = 5000   # Maximum search range for Discriminant D
D_start = 8        #  Starting value for D
found_k1_cutoff = False
found_k2_cutoff = False
found_D_cutoff = False
found_k2_0 = False
count_k1 = 0
count_k2 = 0


# --- Step 1: Determine bound for k1 ---
print(f"{'='*100}")
print("\nDetermining the upper bound for the weight k1...\n")


for i in range(2, Maxnumber_k1+1, 2):
    result = CC(i)
    if result > 0:
        count_k1 = i - 2
        print(f"For the Eisenstein-Cusp case with (2) not inert, the maximum possible valid weight k1 of Eisenstein series is {count_k1}.")
        found_k1_cutoff = True
        break
if not found_k1_cutoff:
    print(f"Warning: Cutoff for k1 not found within Maxnumber={Maxnumber_k1}. Consider increasing range.")


# --- Step 2: Iterate over valid k1 ---
print(f"{'='*100}")
print("\nFor each valid k1, determine the max k2. For every pair (k1,k2), there is a max D, but we only report the global maximum D for each k1...\n")
if found_k1_cutoff:
    for i in range(2, count_k1 + 1, 2):
        found_k2_cutoff = False
        # ---  Determine bound for k2 given k1 ---
        for j in range(2, Maxnumber_k2 + 1, 2):
            result = DD(i,j)
            if result > 0:
                count_k2 = j - 2
                if count_k2 == 0:
                    found_k2_0 = True
                else:
                    print(f"\nWhen k1 = {i}，the maximal k2 is {count_k2}.")
                found_k2_cutoff = True
                break
        if found_k2_0:
            print(f"\nNo eigenform product identity exists for k1 = {i}，k2 = {count_k2}.")
            continue #  Skip this k1
        if not found_k2_cutoff:
            print(f"Warning: Cutoff for k2 not found within Maxnumber = {Maxnumber_k2}. Consider increasing range.")
            break  # No feasible k2 found. Increase search range. Skipping further execution.
        # --- Step 3: Determine bound for D ---
        # We use the smallest k2 (k=2) to find the global maximum for D.
        # Typically, larger k2 implies a tighter (smaller) bound on D.
        k=2
        # Iterate D
        found_D_cutoff = False
        for l in range(5, Maxnumber_D + 1,1):
            result = EE(i,k,l)
            if result > 0:
                # === Boundary Triggered: Backtracking ===
                # The current D is too large. We need to backtrack to find the
                # largest D that satisfies specific conditions.
                if l<=D_start:
                    print(f"For k1 = {i},k2 = {k}, the maximum valid D is {l}, we do not need to consider this case.")
                    found_D_cutoff = True
                    break

                for D_max in range(l-1, D_start-1, -1):
                    if not check_field(D_max):
                        continue
                    # by Prop.4.12, we first consider the case that (2) is not inert, which means                       # D\equiv 1 or --1\mod 8 or 2|D.
                    if D_max % 8 in [1, -1] or D_max % 2 ==0:
                        K.<a> = QuadraticField(D_max)
                        h = K.narrow_class_group().order()
                        if h == 1:
                            print(f"For k1= {i},k2= {k}, the maximum valid D is {D_max}.")
                            found_D_cutoff = True
                            break # End D search for this k1
                break
        if not found_D_cutoff:
            print(f"Warning: Cutoff for D not found within Maxnumber = {Maxnumber_D}. Consider increasing range.")
            break
print(f"{'='*100}")
print("\nVerification process completed.")
︡7416955d-4571-4f23-af0d-2c7ff1477cd9︡{"stdout":"====================================================================================================\n"}︡{"stdout":"We begin dealing with the question of whether the eigenform product identity g=fh exists. In this case, we require one of f or h to be a cusp form and the other an Eisenstein series, and the ideal (2) to be not inert.\n"}︡{"stdout":"====================================================================================================\n"}︡{"stdout":"(k1, k2, D)\n"}︡{"stdout":"\nSymbolic Expression of CC(k1):\n\n"}︡{"html":"<div align='center'>$\\displaystyle -\\frac{\\pi^{5} \\left(2 \\, \\pi\\right)^{2 \\, k_{1} - 1}}{18 \\, \\Gamma\\left(k_{1}\\right)^{2}} + 2^{k_{1}} - 1$</div>"}︡{"stdout":"\nSymbolic Expression of DD(k1,k2):\n\n"}︡{"html":"<div align='center'>$\\displaystyle -\\frac{\\pi^{5} \\left(2 \\, \\pi\\right)^{2 \\, k_{1} - 1}}{18 \\, \\Gamma\\left(k_{1}\\right)^{2}} + {\\left(2^{k_{1}} - 1\\right)} 2^{k_{2} - 1}$</div>"}︡{"stdout":"\nSymbolic Expression of EE(k1,k2,D)\n\n"}︡{"html":"<div align='center'>$\\displaystyle -\\frac{\\pi^{5} \\left(2 \\, \\pi\\right)^{2 \\, k_{1} - 1}}{18 \\, D^{k_{1} - \\frac{1}{2}} \\Gamma\\left(k_{1}\\right)^{2}} + {\\left(2^{k_{1}} - 1\\right)} 2^{k_{2} - 1}$</div>"}︡{"stdout":"====================================================================================================\n"}︡{"stdout":"\nDetermining the upper bound for the weight k1...\n\n"}︡{"stdout":"For the Eisenstein-Cusp case with (2) not inert, the maximum possible valid weight k1 of Eisenstein series is 12."}︡{"stdout":"\n"}︡{"stdout":"====================================================================================================\n"}︡{"stdout":"\nFor each valid k1, determine the max k2. For every pair (k1,k2), there is a max D, but we only report the global maximum D for each k1...\n\n"}︡{"stdout":"\nWhen k1 = 2，the maximal k2 is 10.\nFor k1= 2,k2= 2, the maximum valid D is 73."}︡{"stdout":"\n\nWhen k1 = 4，the maximal k2 is 14.\nFor k1= 4,k2= 2, the maximum valid D is 8.\n\nWhen k1 = 6，the maximal k2 is 14.\nFor k1 = 6,k2 = 2, the maximum valid D is 5, we do not need to consider this case.\n\nWhen k1 = 8，the maximal k2 is 12.\nFor k1 = 8,k2 = 2, the maximum valid D is 5, we do not need to consider this case.\n\nWhen k1 = 10，the maximal k2 is 8.\nFor k1 = 10,k2 = 2, the maximum valid D is 5, we do not need to consider this case.\n\nWhen k1 = 12，the maximal k2 is 2.\nFor k1 = 12,k2 = 2, the maximum valid D is 5, we do not need to consider this case.\n"}︡{"stdout":"====================================================================================================\n"}︡{"stdout":"\nVerification process completed.\n"}︡{"done":true}









