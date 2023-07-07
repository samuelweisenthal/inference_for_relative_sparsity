# norms
# note that it is possible to 
# use an approximation, hardcoded here
# it gave almost identical results 
# to the real norm in practice

l1 = function(x){
 sum(abs(x))
}

l1.app = function(x){sum(sqrt(x**2+1e-50))}
#l1 = l1.app

l2 = function(x){
  sum(x**2)
}



