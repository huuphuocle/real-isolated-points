(*  Check if the square of a box's diameter is smaller than 
    a given parameter e0^2 
    Input:  boxes   :   a list of hypercubes in (x1,...,xn)
            e0      :   the constant e0
*)
boxsizecheck:=proc(boxes,e0,a)
    local box,val,i,e02:
    e02:=e0^2:
    for box in boxes do:
        val:=add((box[i][1]-box[i][2])^2,i=1..nops(box)):
        if(val >= e02) then:
            # error "Need more precision!":
            return false:
        end if:        
    end do:
    return true:    
end proc:

(*  Decide the emptiness of f = 0 with a box B
    Input:  f       :   a polynomial
            box     :   a box
            vars    :   variables
*) 
boxIntersect:=proc(f,box,vars)
    local n,new_f,new_vars,res,b,ineqs,i:
    
    # Remove degenerate sides from the box
    new_vars:=[]:
    new_f:=f:
    for i from 1 to nops(vars) do:
        if box[i][1] = box[i][2] then:
            new_f:=subs(vars[i] = box[i][1],new_f):
        else:
            new_vars:=[op(new_vars),i]:
        end if:
    end do:
    n:=nops(new_vars):
    new_f:=primpart(new_f):

    # b is the set of inequalities defining the simplified box
    b:=[seq(vars[i] > box[i][1],i in new_vars),seq(vars[i] < box[i][2],i in new_vars)]:

    # We check the intersection of f with each side of the box
    for i from 1 to n do:
        ineqs:=[op(b[1..(i-1)]),op(b[(i+1)..(n+i-1)]),op(b[(n+i+1)..-1])]:
        res:=RAG[HasRealSolutions]([new_f=0,vars[new_vars[i]] = box[new_vars[i]][1],op(ineqs)]):
        if res <> [] then:
            return true:
        end if:
        res:=RAG[HasRealSolutions]([new_f=0,vars[new_vars[i]] = box[new_vars[i]][2],op(ineqs)]):
        if res <> [] then:
            return true:
        end if:
    end do:
    return false:
end proc:

(*  Compute boxes that isolate the candidates
    - Without e0    (for the heursitic)
    - With e0       (for the deterministic)
    Input:  cand    :   parametrization
            u       :   slack variable
            a       :   coefficients of distance function
            e0      :   yes here it is
*) 
#   TODO:   Optimize the else part
approximations:=proc(cand,u,a,e0:=0)
    local w,iso,boxes,i,l,e02,prec:
    prec:=10:
    w:=cand[-1]:
    iso:=RootFinding[Isolate](w, digits = prec, constraints = cand[1..-2], output = 'interval'):
    l:=nops(iso[1]):
    boxes:=[seq(map(rhs,iso[2][i]),i=1..l)]:
    if e0 = 0 then:
        return [boxes,iso[1]]:
    fi:
    while not boxsizecheck(box,e0) do:
        prec:=prec*2:
        iso:=RootFinding[Isolate](w, digits = prec, constraints = cand[1..-2], output = 'interval'):
        l:=nops(iso[1]):
        boxes:=[seq(map(rhs,iso[2][i]),i=1..l)]:
    od:
    return [boxes,iso[1]]:
end proc:

(*  Decide for each box from a list boxes if it intersects f = 0
    Input:  f       :   a polynomial
            boxes   :   a list of boxes
            vars    :   variables
*)
verifyCandidates:=proc(f,boxes,vars)
    local i,lindex:
    lindex:=[]:
    for i from 1 to nops(boxes) do:
        if not boxIntersect(f,boxes[i],vars) then:
            lindex:=[op(lindex),i]:
        end if:
    end do:
    return lindex:
end proc:

# HasRealSolutions:=proc(f,vars,ineqs)
#     local n,lf,gbe,gbsolve:
#     n:=nops(vars):
#     if n = 1 then:
#         # need a function to compute for one variable (maybe use subresultant sequence)
#         return true:
#     end if:
#     lf:=add(rand()*vars[i],i=1..n) mod 11:
#     gbe:=FGb[fgb_gbasis_elim]([seq(diff(f-u*lf,vars[i]),i=1..n)],0,[u],vars,{"verb"=3}):
#     gbsolve:=FGb[fgb_matrixn_radical]([op(gbe),f],0,vars,0,{"verb"=3}):
#     # need a function havepoints to check if a zero-dimensional system has solutions in a hypercube
#     if not havepoints(gbsolve,ineqs) then:
#         for i from 1 to n do:
#             HasRealSolutions(subs(vars[i]=,f),vars,ineqs):
#             HasRealSolutions(subs(vars[i]=,f),vars,ineqs):
#         end do:
#     end if:
#     return true:
# end proc: