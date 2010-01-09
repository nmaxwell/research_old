
def order( tuples, bins ):
    
#   symmetries:
#   for a,b in tuples,
#       if a[0] == b[1] and a[1] == b[0] then a == b
#       if a[2] == b[3] and a[3] == b[2] then a == b
#       if a[0] == b[2] and a[1] == b[3] and a[2] == b[0] and a[3] == b[1] then a == b

    popped = False
    
    for n,a in enumerate(C):
        
        if popped:
            break
        
        for m,b in enumerate(C):
            
            if popped:
                break
            
            if a[0] == b[1] and a[1] == b[0] and not a[0] == b[0] and not n == m:
                C_.append( C.pop(m) )
                popped = True
            
            if a[2] == b[3] and a[3] == b[2] and not a[2] == b[2] and not n == m:
                C_.append( C.pop(m) )
                popped = True
            
            if a[0] == b[2] and a[1] == b[3] and not a[2] == b[2] and not n == m:
                C_.append( C.pop(m) )
                popped = True
            
    return popped




N = [1,2,3]

C = []
D = []

for i in N:
    for j in N:
        for k in N:
            for l in N:
                C.append( (i,j,k,l) )


while rm( C, D ):
    continue

print len(C),C
print "\n",len(D),D





