import numpy

test1=numpy.array([[1,1,1,0,0,0,], [0,0,0,1,1,1], [1,1,1,0,0,0]])
test2=numpy.ascontiguousarray(test1).view(numpy.dtype((numpy.void, test1.dtype.itemsize*test1.shape[1])))
_,idx=numpy.unique(test2, return_index=True)
test3=test1[idx]
print test1
print test3

