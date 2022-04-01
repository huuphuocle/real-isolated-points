(*  Elimination variables in computing limits of
    critical points of delta restricted to V(f-eps)

    And a trick to simplify the elimination ideal
*)
elim_delta:=proc(f,vars,cand,u,a,verb:=0)
    local n, gbelim, i, v, syselim:
    n:=nops(vars):
    v:=[seq(cat('v',i),i=1..n)]:
    gbelim:=FGb[fgb_gbasis_elim]([seq(u*diff(f,vars[i])-2*a[i]*(vars[i]-v[i]),i=1..n)],0,[u],[op(vars[1..n]),op(v)],{"verb"=verb}):
    syselim:=[]:
    for i from 1 to nops(gbelim) do:
        if(degree(gbelim[-i],v[1]) = 1) then:
            syselim:=[op(syselim),gbelim[-i]]:
        end if:
    end do:
    syselim:=map(expand,subs([seq(v[i]=cand[i],i=1..n)],syselim)):
    return syselim:
end proc:

(*  Compute the value of e0 by computing critical values of
    the distance function delta restricted to V(f-eps)

    Input:  f       :   a polynomial
            vars    :   variables x1,...,xn
            cand    :   the parametrization of candidates
            u       :   slack variable in cand
            a       :   distance function coefficients
    Output: e0
*)
computeE0:=proc(f,vars,cand,u,a,verb:=0)
    local e, syselim, i, iso, delta, gbe, e0:
    # Distance function delta
    delta:=add(a[i]*(vars[i]-cand[i])^2,i=1..nops(vars)):

    # Elimination step
    syselim:=elim_delta(f,vars,cand,u,a,verb):

    # Solve the final system to get critical values of delta
    # gbe:=FGb[fgb_matrixn_radical]([op(syselim),f,e-delta,cand[-1]],0,[op(vars),:
    gbe:=FGb[fgb_gbasis_elim]([op(syselim),f,e-delta,cand[-1]],0,[op(vars),u],[e],{"verb"=verb}):
    
    # Isolate real roots of the univariate polynomial encoding e
    iso:=map(rhs,RootFinding[Isolate](gbe[1],e,output='interval')):
    if member([0,0],iso,'i') then:
        if nops(iso) > i then:
            e0:=iso[i+1][1]:
            if e0 >= 1 then:
                return 1:
            else:
                return e0: # return something with smaller bit-size
            end if:
        else:
            return 1:
        end if:
    else:
        return []:
    end if:
end proc: