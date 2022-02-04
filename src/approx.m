# check if the diameter of a box is smaller than e02 
boxsize:=proc(boxes,e0,a)
    local box,val,i,e02:
    e02:=e0^2:
    for box in boxes do:
        val:=add((box[i][1]-box[i][2])^2,i=1..nops(box)):
        if(val >= e02) then:
            printf("Need more precision!"):
            return false:
        end if:        
    end do:
    return true:    
end proc:

# compute boxes that isolate the solutions of cand
approximations:=proc(cand,u,a,e0)
    local w,iso,boxes,i,l:
    w:=primpart(cand[-1]):
    iso:=RootFinding[Isolate](w, constraints = cand[1..-2], output = 'interval'):
    l:=nops(iso[1]):
    boxes:=[seq(map(rhs,iso[2][i]),i=1..l)]:
    return [boxes,iso[1]]:
end proc:

# a function to decide emptiness of f = 0 with a box B
boxIntersect:=proc(f,box,vars)
    local n,new_f,new_vars,res,b,ineqs,i:
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
    new_f:=primpart(expand(new_f)):
    # b contains the boundaries of the box
    b:=[seq(vars[i] > box[i][1],i in new_vars),seq(vars[i] < box[i][2],i in new_vars)]:
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

# a function that takes f and a list of boxes and 
# check if those boxes intersect f = 0
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