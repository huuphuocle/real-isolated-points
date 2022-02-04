IsolatedPoints:=module()
description "Computing real isolated points of a given hypersurface";
option package;
export isolatedPoints:
local isolated_, # Isolated.mpl
rewriteParam, candidates, candidates2, # candidate.m
computeE0, # e0.m
boxsize, boxIntersect, approximations, # approx.m
elim_delta, 
Test1, curveLimit: # opti.m

$include "candidate.m":
$include "e0.m"
$include "approx.m":

# the exported function
# output: [cand (parametrization), boxes]
isolatedPoints:=proc(f,verb:=0)
    local vars:
    vars:=convert(indets(f),list):
    return isolated_(f,vars,verb):
end proc:

isolated_:=proc(f,vars,verb:=0)
    local u, n, cand, e0, a, appr, i, I_iso:
    n:=nops(vars):
    # define a linear form to use for generic projection
    lf:=u+add(rand()*vars[i],i=1..(n+1)) mod 11:
    # start by computing candidates
    cand:=candidates(f,vars,lf,verb)[1]:
    if cand = [] then:
        # if there is no candidate, we return
        return [[],[]]:
    else:
        # choose coefficients a for the distance function
        a:=[seq(1,i=1..n)]: # a:=[seq(rand(),i=1..n)]:
        # optimization should go here
        appr:=approximations(cand,u):
        heuristic(f,appr[1],vars):
        # if optimizations failed, we compute e0 (super slow)
        printf("Start computing e0:");
        e0:=computeE0(f,[op(vars),u],cand,a,verb):
        printf("Finish computing e0. Move to identification.");
        # compute the approximation points
        appr:=approximations(cand,u,a,e0):
        boxes:=[]:
        for box in appr[1] do:
            if boxIntersect(f,box,vars) then:
                boxes:=[op(boxes),box]:
            end if:
        end do:
        return [cand, box]:
    end if:
end proc:

end module:

