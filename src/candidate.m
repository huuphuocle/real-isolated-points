(*  Rewrite a rational parametrization into a polynomial parametrization
    Input:  par     :   a rational parametrization of output-type of FGb
            vars    :   variables (x1,...,xn,u)
*)

rewriteParam:=proc(par, vars)
    local n, w, dw, V, i, r, s:
    n:=nops(vars):
    w:=par[-1]:
    dw:=diff(w,vars[-1])/(degree(w)*lcoeff(w)):
    gcdex(w,dw,vars[-1],'r','s');
    V:=[seq(0,i=1..n)]:
    V[-1]:=primpart(w):
    for i from 1 to n-1 do:
        V[i]:=-rem(coeff(par[i],vars[i],0)*s,w,vars[-1])/coeff(par[i],vars[i],1):
    end do:
    return V:
end proc:

(*  Computing the zero-dimensional parametrization of candidates as
    the union of limits of critical points w.r.t. a random projection
    Input:  f       :   polynomial
            vars    :   variables
            u       :   slack variable (for shape position)
            lf      :   linear form of type u + a1*x1 + ... + an*xn
*)
candidates:=proc(f,vars,u,lf,verb:=0)    
    local J,i,j,n,par_cand,sys,gbsolve,gbs1,par,iso,lfA,test_cand:
    description "Computes the candidates through limits of critical points.":
    n:=nops(vars):

    # Prepare a list of n random linear forms for projecting
    A:=LinearAlgebra[RandomMatrix](n,n,generator=rand(1..19));
    while LinearAlgebra[Determinant](A) = 0 do:
        A:=LinearAlgebra[RandomMatrix](n,n):
    od:
    lfA:=LinearAlgebra[Multiply](A,Vector(vars)):
    J:=[seq([],i=1..n)]:

    # Elimination step to compute limits
    for i from 1 to nops(J) do:
        J[i]:=FGb[fgb_gbasis_elim]([seq(diff(u*f-lfA[i],vars[j]),j=1..n)],0,[u],vars,{"verb"=verb}):
    od:

    # Compute the parametrization of candidates
    gbsolve:=FGb[fgb_matrixn_radical]([f,seq(op(J[i]),i=1..nops(J)),lf],0,[op(vars),u],0,{"verb"=verb}):
    if gbsolve=[] then:         # No candidate
        return [[],[]]:
    else:                       # Rewrite the candidates
        iso := RootFinding[Isolate](gbsolve[-1],u,output='interval'):   # Isolate u-coordinates of candidates
        if nops(iso) = 0 then:  # No candidate
            return [[],[]]:
        else:                   
            par_cand:=rewriteParam(gbsolve,[op(vars),u]):   # The parametrization of candidates
            # Compute the test set for candidates
            gbs1:=FGb[fgb_matrixn_radical]([op(J[1]),f,lf],0,[op(vars),u],0,{"verb"=verb}):
            test_cand:=rewriteParam(gbs1,[op(vars),u]):
            return [par_cand,iso,test_cand]:
        end if:
    end if:
end proc: