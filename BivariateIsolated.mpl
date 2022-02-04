# limit of a curve
IsolateBivariate:=proc(w)
    local gb, isol:
    gb:=fgb_matrixn_radical([w,diff(w,x),diff(w,y)],0,[x,y],0):
    RootFinding[Isolate](gb[-1]):
    return isol:
end proc:

CurveLimit:=proc(f, vars)
    local n, gbe, gb, l, u, v:
    n:=nops(vars):
    gbe:=FGb[fgb_gbasis_elim]([seq(diff(f,vars[i]),i=3..n),l*(u*diff(f,vars[1])-1),(1-l)*(v*diff(f,vars[2])-1)],0,[l,u,v],vars):
    gb:=FGb[fgb_gbasis]([op(gbe),f],0,vars[1..-3],vars[-2..-1]):
    return gb:
    # gb1:=FGb[fgb_gbasis]([op(gbe1),],0,vars[1..-3],vars[-2..-1]):
    # gbe1:=FGb[fgb_gbasis_elim]([seq(diff(f,vars[i]),i=3..n],0,[u],vars):
    # gb1:=FGb[fgb_gbasis]([op(gbe1),],0,vars[1..-3],vars[-2..-1]):
    # gbe2:=FGb[fgb_gbasis_elim]([seq(diff(f,vars[i]),i=3..n],0,[u],vars):
    # gb2:=FGb[fgb_gbasis]([op(gbe2),op()],0,vars[1..-3],vars[-2..-1]):
end proc:

CurveLimit:=proc(f, vars)
    local n, gbe, gbe0, gb, u:
    n:=nops(vars):
    gbe:=FGb[fgb_gbasis_elim]([seq(diff(f,vars[i]),i=3..n),u*diff(f,vars[2])-1],0,[u],vars,{"verb"=3}):
    gbe:=[seq(u*gbe[i],i=1..nops(gbe))]:
    gbe0:=FGb[fgb_gbasis_elim]([seq(diff(f,vars[i]),i=2..n),u*diff(f,vars[1])-1],0,[u],vars,{"verb"=3}):
    gbe0:=[seq(u*gbe0[i],i=1..nops(gbe0))]:
    gb:=FGb[fgb_gbasis_elim]([op(gbe),op(gbe0),f],0,[u],vars,{"verb"=3}):
    gb:=FGb[fgb_gbasis]([op(gbe),op(gbe0),f],0,vars[1..-3],vars[-2..-1]):
    return gb:
    # gb1:=FGb[fgb_gbasis]([op(gbe1),],0,vars[1..-3],vars[-2..-1]):
    # gbe1:=FGb[fgb_gbasis_elim]([seq(diff(f,vars[i]),i=3..n],0,[u],vars):
    # gb1:=FGb[fgb_gbasis]([op(gbe1),],0,vars[1..-3],vars[-2..-1]):
    # gbe2:=FGb[fgb_gbasis_elim]([seq(diff(f,vars[i]),i=3..n],0,[u],vars):
    # gb2:=FGb[fgb_gbasis]([op(gbe2),op()],0,vars[1..-3],vars[-2..-1]):
end proc:

RewriteParam:=proc(par, vars)
    local n, w, dw, V, i, r, s:
    n:=nops(vars):
    w:=par[-1]:
    V:=[seq(0,i=1..n)]:
    V[-1]:=w:
    for i from 1 to n-1 do:
        V[i]:=-coeff(par[i],vars[i],0)/coeff(par[i],vars[i],1):
    end do:
    return V:
end proc:

BiCandidates:=proc(f, vars, verb:=0)
    local J1, J2, i, cand, u, gbsolve, par:

    J1:=FGb[fgb_gbasis_elim]([diff(f,vars[1]),u*diff(f,vars[2])-1],0,[u],vars,{"verb"=verb}):
    J2:=FGb[fgb_gbasis_elim]([diff(f,vars[2]),u*diff(f,vars[1])-1],0,[u],vars,{"verb"=verb}):

    gbsolve:=FGb[fgb_matrixn_radical]([f,op(J1),op(J2)],0,vars,0,{"verb"=verb}):

    cand:=RewriteParam(gbsolve,vars):
    if nops(RootFinding[Isolate](cand[-1],vars[-1])) = 0 then:
        return []:
    else:
        return cand:
    end if:
end proc:


IsolatedBivariate:=proc(f,vars,verb:=0)
    local u, t, cand, dist, gbe:
    cand:=subs(vars[-1]=t,BiCandidates(f,vars,verb)):
    gbelim:=FGb[fgb_gbasis_elim]([u*diff(f,vars[1])-(vars[1]-vt),u*diff(f,vars[2])-(vars[2]-t)],0,[u],[op(vars),vt,t],{"verb"=verb}):
    gbelim:=expand(subs(vt = cand[1],gbelim)):
    gbe:=FGb[fgb_matrixn_radical]([e-(vars[1]-cand[1])^2-(vars[2]-t)^2,cand[-1],op(gbelim),f],0,[op(vars),t,e],0):
    e:=RootFinding[Isolate](gbe[-1],output='interval'):
    return e:
end proc:

vars:=[x1,x2]:
f:=randpoly(vars,degree=3):
gbe:=fgb_gbasis_elim([t*f,(1-t)*vars[1],(1-t)*vars[2]],0,[t],vars):
f:=expand(add(gbe[i]^2,i=1..nops(gbe))):