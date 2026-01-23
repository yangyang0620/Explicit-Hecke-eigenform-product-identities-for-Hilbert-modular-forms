# Explicit-Hecke-eigenform-product-identities-for-Hilbert-modular-forms

# Overview
This repository contains **SageMath** scripts and **Magma** code developed to verify the numerical results and validate Theorem 1 of the paper:
             **Explicit Hecke eigenform product identities for Hilbert modular forms**

# Methodology
The core strategy for verifying the product identities relies on the comparison of Fourier coefficients and their analytic properties. The verification process follows these logical steps:
1. Deriving Inequalities: By analyzing the equality relations of Fourier coefficients (for the identity $f⋅h=g$) and the analytic bounds of the associated functions, we construct necessary inequalities involving the weights 
$k_1,k_2$ and the discriminant $D$.
2. Establishing Bounds: These inequalities allow us to determine upper bounds for $k_1$, $k_2$, and $D$, effectively defining a finite search space.
3. Exhaustive Verification: The program traverses all possible tuples $(k_1,k_2,D)$ within these bounds.  
(1) **Eisenstein-Eisenstein Case**: Verification is performed by checking the exact rational values of the Fourier coefficients (specifically the constant term).  
(2) **Eisenstein-Cusp Case**: Verification is performed based on dimensional constraints of the space of cusp forms.  


# File Structure & Description
The codebase consists of four SageMath scripts corresponding to different cases of the modular forms, plus an independent Magma code for dimension computation.

1. The case of Eisenstein-Eisenstein of distinct weights  
Mathematical Basis: **Propositions 3.3 and 3.4**.  
Logic:  
Initializes Real Quadratic Fields (e.g., $D=8,13,…$).  
Defines lower bound functions derived from the coefficient formulas (Eq 3.2).  
  Phase 1: Determines a "Cutoff" for $k_2$. For any $k_1 > k_2 \geq$ Cutoff, the lower bound exceeds 1, implying no identities exist.  
  Phase 2: Exhaustively searches the range $k_2< $ Cutoff. If the bound allows, it verifies the existence using the constant term identity (Eq 3.1).  
Result: This script reproduces the proof of **Propositions 3.4**.  

2. The case of Eisenstein-Eisenstein of equal weights  
Mathematical Basis: **Proposition 3.8**.  
Logic:  
Determines the upper bound for the weight k.  
Iterates through discriminants $D$ to find satisfying fields and verifies the identities.  
Result: This script reproduces **Table 1** of the paper, verifying the first part of Theorem 1.  

3. The case of Eisenstein-Cusp with (2) inert  
Mathematical Basis: **Proposition 4.12**.  
Logic:  
Determines the upper bound for the weight $k_1$.
Iterates through each fixed $k_1$ to find the upper bound for $k_2$. For each pair $(k_1,k_2)$, finds the upper bound for $D$.  
Result: This script reproduces **Table 2** of the paper, verifying the second part of Theorem 1.  

4. The case of Eisenstein-Cusp with (2) not inert  
Mathematical Basis: **Proposition 4.13**.  
Logic:  
Similar to the inert case, determines the upper bound for $k_1$.
Iterates through $k_1$ to limit $k_2$, and subsequently limits $D$.  
Result: These results reproduce **Table 3** of the paper, verifying the second part of Theorem 1.  

5. Dimension of space of Hilbert cusp form (Magma)  
Platform: [Magma Online Calculator](http://magma.maths.usyd.edu.au/calc/)  
Description: It calculates the dimension of the space of Hilbert cusp forms $S_k(\Gamma)$ for real quadratic fields with full level.  
Usage: It is used to verify the **dimensional reasons** in the Eisenstein-Cusp cases (Section 4). Specifically, it confirms the bounds and non-existence results (e.g., checking where $\dim S_k(\Gamma) > 1$) referenced in the proofs for **Table 3**.




