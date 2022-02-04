# This file contains the heuristic tests for 
# avoiding the heavy computation (e0).

# if a box does not intersect f, 
# that box corresponds to an isolated point
heuristic:=proc(f,boxes,vars)
    local iso,box,i,j,b,c,points:
    for box in boxes do:
        points:=[]:
        # check if box intersects f
        b := [seq(vars[i] > box[i][1],i=1..n),seq(vars[i] < box[i][2],i=1..n)]:
        for i from 1 to n do:
            c := [op(b[1..(i-1)]),op(b[(i+1)..(n+i-1)]),op(b[(n+i+1)..-1])]:
            points:=[op(points),op(RAG[HasRealSolutions]([f = 0, vars[i] = box[i][1], op(c)]))]:
        end do:
    end do:
    if points = [] then:
        return true:
    return points:
end proc:

# compute a parametrization for the limit of crit(pi_{1,2})
curveLimit:=proc(f, vars)
    local n, gbe1, gbe2, gbe, u, v, gb:
    n:=nops(vars):
    gbe1:=FGb[fgb_gbasis_elim]([seq(diff(f,vars[i]),i=3..n),u*diff(f,vars[1])-1],0,[u],vars):
    gbe2:=FGb[fgb_gbasis_elim]([seq(diff(f,vars[i]),i=3..n),u*diff(f,vars[2])-1],0,[u],vars):
    gbe:=FGb[fgb_gbasis_elim]([seq(u*v,v in gbe1),seq((1-u)*v,v in gbe2)],0,[u],vars):
    gb:=FGb[fgb_gbasis_elim]([op(gbe),f],0,vars[1..-3],vars[-2..-1]):
    return gb:
end proc: