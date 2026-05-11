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

    
    Call:
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
    (Intercept)        -18.8099    12.0952 -1.5552 0.1199092    
    treatmentprogabide  -1.2852     2.0996 -0.6121 0.5404658    
    lnbaseline          10.2050     2.9740  3.4315 0.0006003 ***
    lnage                2.9236     3.9501  0.7401 0.4592099    
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    Std. Errors are taken from the bias-corrected covariance matrix. 
    
    Dispersion Parameter: 97.2764 
    
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
    <none>         6.5258                     
    lnbaseline  1 17.5279 11.515 0.0006902 ***
    lnage       1  9.2547  0.245 0.6205917    
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

---

    Single term deletions using Wald test:
    
    Model:
    seizures ~ treatment + lnbaseline + lnage
               Df     CIC     Chi  Pr(>Chi)    
    <none>        18.9567                      
    treatment   1 16.7365  0.4443 0.5050445    
    lnbaseline  1  9.2547 12.5073 0.0004054 ***
    lnage       1 17.5279  0.6994 0.4029695    
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

