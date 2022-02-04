elim_delta:=proc(f,vars,cand,a,verb:=0)
    local n, u, gbelim, i, v, syselim:
    n:=nops(vars)-1:
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

computeE0:=proc(f,vars,cand,a,verb:=0)
    local n, u, e, syselim, i, iso, delta, gbe, e0:
    n:=nops(vars)-1:

    syselim:=elim_delta(f,vars,cand,a,verb):
    # gbe:=FGb[fgb_matrixn_radical]([op(syselim),f,e-delta,cand[-1]],0,[op(vars),:
    delta:=add(a[i]*(vars[i]-cand[i])^2,i=1..n):
    gbe:=FGb[fgb_gbasis_elim]([op(syselim),f,e-delta,cand[-1]],0,vars,[e],{"verb"=verb}):
    
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