# rewrite the parametrization 
rewriteParam:=proc(par, vars)
    local n, w, dw, V, i, r, s:
    n:=nops(vars):
    w:=par[-1]:
    dw:=diff(w,vars[-1])/(degree(w)*lcoeff(w)):
    gcdex(w,dw,vars[-1],'r','s');
    V:=[seq(0,i=1..n)]:
    V[-1]:=w:
    for i from 1 to n-1 do:
        V[i]:=-rem(coeff(par[i],vars[i],0)*s,w,vars[-1])/coeff(par[i],vars[i],1):
    end do:
    return V:
end proc:

# Compute the candidates using critical point method
candidates:=proc(f,vars,lf,verb:=0)
    local J, grad, i, n, cand, u, sys, gbsolve, par, iso:
    description "This function computes the candidates.":
    printf("Start computing candidates:\n");
    n:=nops(vars):
    J:=[seq([],i=1..n)]:
    grad:=[seq(diff(f,vars[i]),i=1..n)]:

    for i from 1 to n do:
        sys:=grad:
        sys[i]:=u*grad[i]-1:
        # printf("Eliminate for %d\n",i);
        J[i]:=FGb[fgb_gbasis_elim](sys,0,[u],vars,{"verb"=verb}):
    end do:
    # printf("Solve the system:\n");
    # gbsolve:=FGb[fgb_gbasis]([f,seq(op(J[i]),i=1..n),u-lf],0,[],[op(vars),u],{"verb"=verb}):
    gbsolve:=[]:
    gbsolve:=FGb[fgb_matrixn_radical](gbsolve,0,[op(vars),u],0,{"verb"=verb}):
        
    if gbsolve=[] then:
        # printf("There is no candidate.\n");
        return []:
    else:
        # printf("Rewrite the candidates.\n");
        iso := RootFinding[Isolate](gbsolve[-1],u,output='interval'):
        if nops(iso) = 0 then:
            # printf("There is no candidate.\n");
            return []:
        else:
            cand:=rewriteParam(gbsolve,[op(vars),u]):
            # printf("Done.\n");
            return [cand,iso]:
        end if:
    end if:
end proc:

# Computing the zero-dimensional parametrization of the union 
# of limits of critical points w.r.t. all projections on x[i]
candidates2:=proc(f, vars, lf, verb:=0)
    local J, grad, i, n, cand, u, sys, iso, gbsolve:
    
    n:=nops(vars):
    J:=[seq([],i=1..n)]:
    grad:=[seq(diff(f,vars[i]),i=1..n)]:

    for i from 1 to n do:
        sys:=grad:
        sys[i]:=u*grad[i]-1:
        J[i]:=FGb[fgb_gbasis_elim](sys,0,[u],vars,{"verb"=verb}):
        gbsolve:=FGb[fgb_matrixn_radical]([op(J[i]),f,u-lf],0,[op(vars),u],0,{"verb"=verb}):
        iso:=RootFinding[Isolate](gbsolve[-1],output='interval'):
    end do:
    return iso:
end proc: