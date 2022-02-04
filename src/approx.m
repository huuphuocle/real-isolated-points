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
# and check their diagonals
approximations:=proc(cand,u)
    local w,iso,box,boxes,i:
    w:=primpart(cand[-1]):
    iso:=RootFinding[Isolate](w, constraints = cand[1..-2], output = 'interval'):
    boxes:=[seq(map(rhs,iso[2][i]),i=1..nops(iso[1]))]:
    return [boxes,iso[1]]:
end proc:

# a function to decide emptiness of f = 0 with a box B
boxIntersect:=proc(f,box,vars)
    local n, res, eqs, ineqs, i:
    n:=nops(vars):
    # b contains the boundaries of the box
    b := [seq(vars[i] > box[i][1],i=1..n),seq(vars[i] < box[i][2],i=1..n)]:
    for i from 1 to n do:
        eqs:= [vars[i] = box[i][1]]:
        ineqs:=[vars[j] > B[j][1] , vars[j] < B[j][2]]:
        res:=RAG[HasRealSolutions]([f=0,op(eqs),op(ineqs)]):
        if res <> [] then:
            return true:
        end if:
    end do:
    return false:
end proc: