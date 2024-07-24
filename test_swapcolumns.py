
import numpy as np

def swap_columns(array, column1, column2):
        if array.ndim == 1:
            temp = array[column1].copy()
            array[column1]  = array[column2]
            array[column2] =  temp
        elif array.ndim == 3:
            array = np.swapaxes(array,1,2)
        return array

test = np.array(['a', 'b', 'c'])


swap_columns(test,1,2)

print(test)