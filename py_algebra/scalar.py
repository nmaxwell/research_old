
"""
test:



from scalar import *
a = scalar()
b = scalar()



"""








class scalar:
    
    def __init__(self ):
        self.scalar_sum = None
        self.scalar_product = None
        
        self.scalar_value = None
        self.scalar_name = None
    
    def copy():
        foo = scalar()
        if rhs.scalar_sum:
            foo.scalar_sum = 
    
    def __add__(self, rhs ):
        
        if self.scalar_sum:
            if rhs.scalar_sum:
                return self.scalar_sum + rhs.scalar_sum
            else:
                return self.scalar_sum + [rhs]
        else:
            if rhs.scalar_sum:
                return [ self ]  + rhs.scalar_sum
            else return [ self, rhs]
        
    
    












