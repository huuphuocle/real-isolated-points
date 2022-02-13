IsolatedPoints:=module()
description "Computing real isolated points of a given hypersurface";
option package;
export isolatedPoints:
local isolated_, # Isolated.mpl
rewriteParam, candidates, candidates2, # candidate.m
elim_delta, computeE0, # e0.m
boxsize, boxIntersect, approximations, verifyCandidates: # approx.m


$include "src/candidate.m":
$include "src/approx.m":
$include "src/e0.m":

# the exported function
# output: [cand (parametrization), boxes]
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
    #st:=time():
    cand:=candidates(f,vars,u,lf,verb)[1]:
    #print(time()-st);
    if cand = [] then:
        # if there is no candidate, we return
        return [[],[]]:
    else:
        # choose coefficients a for the distance function
        a:=[seq(1,i=1..n)]: # a:=[seq(rand(),i=1..n)]:
        # optimization should go here
        # approximations needs to compute boxes suitable for heuristic tests
        appr:=approximations(cand,u):
        l:=nops(appr[2]):
        # if a box does not intersect f, 
        # that box corresponds to an isolated point
        verified_candidates:=verifyCandidates(f,appr[1],vars):
        if nops(verified_candidates) = l then:
            return [cand,appr[1]]:
        end if:
        # if optimizations failed, we compute e0 (super slow)
        #printf("Start computing e0:");
        e0:=computeE0(f,vars,cand,u,a,verb):
        #printf("Finish computing e0. Move to identification.");
        # compute the approximation points with known e0
        appr:=approximations(cand,u,a,e0):
        verified_candidates:=verifyCandidates(f,appr[1],vars):
        boxes:=[seq(appr[1][i],i in verified_candidates)]:
        return [cand, boxes]:
    end if:
end proc:

end module:

