# core printed outputs remain stable

    Call:
    geewa(formula = seizures ~ treatment + lnbaseline + lnage, data = epilepsy, 
        id = id, repeated = visit)
    
    Coefficients:
           (Intercept) treatmentprogabide         lnbaseline              lnage 
            -18.809872          -1.285207          10.205033           2.923632 
    
    Number of iterations : 1 
    Algorithm converged  : TRUE 

---

    
    call:
    geewa(formula = seizures ~ treatment + lnbaseline + lnage, data = epilepsy, 
        id = id, repeated = visit)
    
    Estimating Method   : gee 
    Number of iterations: 1 
    Algorithm converged : TRUE 
    
    Marginal Model
    Family       : gaussian 
    Link Function: identity 
    
    Coefficients:
                       Estimate Std. Error z value  Pr(>|z|)    
    (Intercept)        -18.8099    10.3183 -1.8230 0.0683084 .  
    treatmentprogabide  -1.2852     1.9281 -0.6666 0.5050445    
    lnbaseline          10.2050     2.8856  3.5366 0.0004054 ***
    lnage                2.9236     3.4958  0.8363 0.4029695    
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    Std.Errors are taken from the robust covariance matrix. 
    
    Dispersion Parameter: 102.5824 
    
    Association Structure: independence 
    Association Parameter: 0 

---

    Analysis of Wald Statistic Table
    
    Model: gaussian, link: identity
    
    Response: seizures
    
    Terms added sequentially (first to last)
    
    
               Df Resid. Df     Chi  Pr(>Chi)    
    NULL                235                      
    treatment   1       234  0.0490 0.8248146    
    lnbaseline  1       233 11.5155 0.0006902 ***
    lnage       1       232  0.6994 0.4029695    
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# stepwise printed outputs remain stable

    Single term additions using Wald test:
    
    Model:
    seizures ~ treatment
               Df     CIC    Chi  Pr(>Chi)    
    <none>         6.4700                     
    lnbaseline  1 17.0765 11.515 0.0006902 ***
    lnage       1  9.0164  0.245 0.6205917    
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

---

    Single term deletions using Wald test:
    
    Model:
    seizures ~ treatment + lnbaseline + lnage
               Df     CIC     Chi  Pr(>Chi)    
    <none>        17.9762                      
    treatment   1 16.3055  0.4443 0.5050445    
    lnbaseline  1  9.0164 12.5073 0.0004054 ***
    lnage       1 17.0765  0.6994 0.4029695    
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

