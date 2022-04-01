IsolatedPoints:=module()
description "Computing real isolated points of a given hypersurface";
option package;
export isolatedPoints:
local isolated_, # Isolated.mpl
rewriteParam, candidates, # candidate.m
elim_delta, computeE0, # e0.m
boxsizecheck, boxIntersect, approximations, verifyCandidates: # approx.m


$include "src/candidate.m":
$include "src/approx.m":
$include "src/e0.m":

(*  Exported function: Wrapper of isolated_
    Input:  f   :   a polynomial in x1,...,xn  
    Output: [parametrization_cand, boxes]
*)
isolatedPoints:=proc(f,verb:=0)
    local vars:
    vars:=convert(indets(f),list):
    return isolated_(f,vars,verb):
end proc:

isolated_:=proc(f,vars,verb:=0)
    local u,lf,n,cand,a,appr,l,verified_candidates,e0,boxes,i,st,roll:
    n:=nops(vars):
    roll:=rand(2..19):
    # define a linear form to use for generic projection
    lf:=u+add(roll()*vars[i],i=1..n):
    # start by computing candidates

    # cand = [par_cand,iso,test_cand]
    cand:=candidates(f,vars,u,lf,verb):

    if cand[1] = [] then: # No candidate
        return [[],[]]:
    else:
        # Coefficients a for the distance function
        a:=[seq(1,i=1..n)]: # a:=[seq(rand(),i=1..n)]:
        
        # Heuristic tests
        appr:=approximations(cand[1],u,a):  # Boxes for heuristic tests
        l:=nops(appr[2]):

        # If a box does not intersect f, 
        # that box corresponds to an isolated point
        verified_candidates:=verifyCandidates(f,appr[1],vars):

        if nops(verified_candidates) = l then:
            return [cand[1],appr[1]]:   # Expected to finish here
        end if:

        # If heuristic failed, we compute e0 (super slow)
        e0:=computeE0(f,vars,cand[1],u,a,verb):

        # Deterministic identification
        appr:=approximations(cand[1],u,a,e0): # Boxes with known e0
        verified_candidates:=verifyCandidates(f,appr[1],vars):
        boxes:=[seq(appr[1][i],i in verified_candidates)]:
        return ['param' = cand[1], 'boxes' = boxes]:
    end if:
end proc:

end module:

