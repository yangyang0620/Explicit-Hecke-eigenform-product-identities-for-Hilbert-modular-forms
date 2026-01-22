︠50e2b4af-cbd2-4603-a49a-15aee447e4ab︠
##############################################################################
# Title: Numerical Verification of Eigenform Product Identities (Eisenstein-Eisenstein Case with Equal Weight)
# Description:
#   Based on Propositions 3.3-3.4, this script verifies the non-existence of
#   eigenform identities for the distinct-weight case.
#
#   Workflow:
#   1. Initialize Real Quadratic Field (e.g., D=8, 13, 17, 29, 37).
#   2. Define lower bound functions C1, C2, and  C_full_bound
#      derived from the g=fh coefficient formulas (Eq 3.2).
#   3. Phase 1: Determine a "Cutoff" for k2. For any k2 >= Cutoff, the
#      lower bound > 1, guaranteeing no identities exist.
#   4. Phase 2: Exhaustively search the range k2 < Cutoff.
#      - If lower bound > 1: Discard the pair.
#      - If lower bound <= 1: Verify using the constant term identity (Eq 3.1).
##############################################################################

##############################################################################
# The results of this program  help to verify the first part of Theorem 1.
##############################################################################

# ==========================================
# 1. Setup and Mathematical Definitions
# ==========================================

# Set the parameter for the Real Quadratic Field
d=8  # Modify this value (e.g.8, 13, 17, 29, 37)
K.<a> = QuadraticField(d)
Discriminant = K.discriminant()


print(f"{'='*100}")
print(f"Initializing Real Quadratic Field: Q(sqrt({Discriminant}))")
print("Discriminant:", Discriminant)
print("Is Totally Real:", K.is_totally_real())
print(f"{'='*100}\n")

# ------------------------------------------
# Helper Functions
# ------------------------------------------

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
#1 = | (zeta_F(1-k1) + zeta_F(1-k2)) * zeta_F(1-k1-k2) / (zeta_F(1-k1)*zeta_F(1-k2)) |
def check_equation(a, b, c):
    return c * (a + b) == a * b


# We determine the minimal k2 for which the term inside the absolute value is positive, using a lower bound function that is strictly increasing in k2, see line. -6 on page 7.
#
def C1(D,k2):
    term1 = (291600/pi^(12))
    term2 = (D * k2^2 / (4 * pi^2))^(2)
    return term1 * term2 - 1


#Simplified lower bound function dependent on D, k1, and k2.
def C2(D,k1,k2):
    term1 = (291600/pi^(12))
    term2 = (D /(4 * pi^2))^(k1 - k2)*(gamma(k1) / gamma(k2))^2
    return term1 * term2 - 1


#The precise lower bound function C(D, k1, k2) derived from Eq (3.2).
def C_full_bound(D,k1,k2):
    term = (zeta(4 * (k1 + k2))/zeta(k1 + k2)^2/zeta(k1)^2) * (D/(4 * pi^2))^(k2) * (gamma(k1 + k2)/gamma(k1))^2 * abs((zeta(4 *(k1))/zeta(k1)^2/zeta(k2)^2) * (D/(4 * pi^2))^(k1-k2) * (gamma(k1)/gamma(k2))^2-1)
    return term.n()



# Display symbolic expression to verify correctness
var('D k1 k2')
print("\n Preview of Symbolic Expression for C1:")
show(C1(D, k2))
print("\n Preview of Symbolic Expression for C2:")
show(C2(D, k1, k2))
print("\n Preview of Symbolic Expression for C(D,k1,k2):")
show(((zeta(4 * (k1 + k2))/zeta(k1 + k2)^2/zeta(k1)^2) * (D/(4 * pi^2))^(k2) * (gamma(k1 + k2)/gamma(k1))^2 * abs((zeta(4 *(k1))/zeta(k1)^2/zeta(k2)^2) * (D/(4 * pi^2))^(k1-k2) * (gamma(k1)/gamma(k2))^2-1)))





# ==========================================
# 2. Main Program Logic
# ==========================================


D_val = Discriminant
Maxnumber = 100 # Maximum search range for weights
k2_start = 2 # Starting value for k2 (even integer)

# ------------------------------------------
# Phase 1: Determine the Global Cutoff for k2
# ------------------------------------------
# Find the smallest k2 such that for all k1 > k2, the lower bound > 1.
# Since the function is increasing in k1, checking k1_min = k2 + 2 is sufficient.

print(f"{'='*100}")
print("Phase 1: Determining the cutoff threshold for k2...")
found_cutoff = False
for k2 in range(k2_start, Maxnumber + 2, 2):
    result_C1 = C1(D_val, k2)
    result_C = C_full_bound(D_val, k2 + 2, k2)
    print(f" k2 = {k2} yields C1 = {result_C1.n()}.")
    if result_C1.n() > 0 and result_C.n() > 0:
        print(f"The value of C_full_bound(D,k1,k2) at D = {D_val} k1 = {k2+2} k2 = {k2} is: {result_C.n()} > 0. Implies no eigenform identities exist for k1 > k2 = {k2}.")
        cutoff_k2 = k2
        found_cutoff=True
        print(f"{'='*100}")
        break
if not found_cutoff:
    print(f"Warning: Cutoff for k2 not found within Maxnumber={Maxnumber}. Consider increasing range.")




# ------------------------------------------
# Phase 2: Detailed Search Below the Cutoff
# ------------------------------------------
print(f"Phase 2: Exhaustive verification for k2 < {cutoff_k2}...")

# Iterate through k2 up to the cutoff
if found_cutoff:
    for k2 in range(k2_start, cutoff_k2, 2):
        # Inner loop: Iterate through k1
        # Since the bound increases with k1, if we hit a failure, we can break early.
        found_k1 = False
        for k1 in range(k2 + 2, Maxnumber + 2, 2):
            result_C2 = C2(D_val, k1, k2)
            result_CC = C_full_bound(D_val, k1, k2)
            print("\n Value of C2(D,k1,k2) at D=", Discriminant, "k1=", k1, "k2=", k2, ":", result_C2.n())
            if result_C2.n()>0 and result_CC.n()>0:
                print(f"\nValue of C(D,k1,k2) at D={D_val} k1={k1} k2={k2} is: {result_CC.n()} > 0. Implies no eigenform identities exist for k2={k2}, k1 >= {k1}.")
                found_k1 = True
                break
        if not found_k1:
            print(f"Warning: Cutoff of k1 not found within Maxnumber={Maxnumber}. Consider increasing range.")
        if found_k1:
            if k1 - k2 == 2:
                    break
            else:
                for count_k1 in range(k2 + 2, k1, 2):
                    result_CCC = C_full_bound(D_val, count_k1, k2)
                    if result_CCC.n() <= 0:
                         # Bound <= 1. This is a "Candidate Triple" that the bound cannot rule out.
                        # We must verify it using the constant term identity (Eq 3.1).
                        if check_equation(Zeta(1 - count_k1,1 - k2,1 - count_k1 - k2)):
                            print(f"For D={Discriminant}, k1={count_k1}, k2={k2}, the eigenform product identity holds.")
                    else:
                        print(f"For D={Discriminant}, k1={count_k1}, k2={k2}, C(D,k1,k2)={result_CCC.n()},the eigenform product identity does not hold.")
print(f"{'='*100}")
print("\nVerification process completed.")
︡11cb8c79-563b-4497-80ae-4de160062678︡{"stdout":"====================================================================================================\n"}︡{"stdout":"Initializing Real Quadratic Field: Q(sqrt(8))\n"}︡{"stdout":"Discriminant: 8\n"}︡{"stdout":"Is Totally Real: True\n"}︡{"stdout":"====================================================================================================\n\n"}︡{"stdout":"(D, k1, k2)\n"}︡{"stdout":"\n Preview of Symbolic Expression for C1:\n"}︡{"html":"<div align='center'>$\\displaystyle \\frac{18225 \\, D^{2} k_{2}^{4}}{\\pi^{16}} - 1$</div>"}︡{"stdout":"\n Preview of Symbolic Expression for C2:\n"}︡{"html":"<div align='center'>$\\displaystyle \\frac{291600 \\, \\left(\\frac{D}{4 \\, \\pi^{2}}\\right)^{k_{1} - k_{2}} \\Gamma\\left(k_{1}\\right)^{2}}{\\pi^{12} \\Gamma\\left(k_{2}\\right)^{2}} - 1$</div>"}︡{"stdout":"\n Preview of Symbolic Expression for C(D,k1,k2):\n"}︡{"html":"<div align='center'>$\\displaystyle \\frac{\\left(\\frac{D}{4 \\, \\pi^{2}}\\right)^{k_{2}} {\\left| \\frac{\\left(\\frac{D}{4 \\, \\pi^{2}}\\right)^{k_{1} - k_{2}} \\Gamma\\left(k_{1}\\right)^{2} \\zeta(4 \\, k_{1})}{\\Gamma\\left(k_{2}\\right)^{2} \\zeta(k_{1})^{2} \\zeta(k_{2})^{2}} - 1 \\right|} \\Gamma\\left(k_{1} + k_{2}\\right)^{2} \\zeta(4 \\, k_{1} + 4 \\, k_{2})}{\\Gamma\\left(k_{1}\\right)^{2} \\zeta(k_{1} + k_{2})^{2} \\zeta(k_{1})^{2}}$</div>"}︡{"stdout":"====================================================================================================\n"}︡{"stdout":"Phase 1: Determining the cutoff threshold for k2...\n"}︡{"stdout":" k2 = 2 yields C1 = -0.792714210254152."}︡{"stdout":"\n k2 = 4 yields C1 = 2.31657263593357.\nThe value of C_full_bound(D,k1,k2) at D = 8 k1 = 6 k2 = 4 is: 186576.531691153 > 0. Implies no eigenform identities exist for k1 > k2 = 4.\n====================================================================================================\n"}︡{"stdout":"Phase 2: Exhaustive verification for k2 < 4...\n"}︡{"stdout":"\n Value of C2(D,k1,k2) at D= 8 k1= 4 k2= 2 : -0.533606973071842\n\n Value of C2(D,k1,k2) at D= 8 k1= 6 k2= 2 : 6.66077206104477\n\nValue of C(D,k1,k2) at D=8 k1=6 k2=2 is: 532.503010432984 > 0. Implies no eigenform identities exist for k2=2, k1 >= 6.\nFor D=8, k1=4, k2=2, C(D,k1,k2)=7.22915264352208,the eigenform product identity does not hold.\n"}︡{"stdout":"====================================================================================================\n"}︡{"stdout":"\nVerification process completed.\n"}︡{"done":true}









