# importing library to have the number pi #####################
import math

# define get_sphere_volume() function #########################
def get_sphere_volume(radius):
    '''Returns the volume of the sphere with the specified radius'''
    volume = 4/3*math.pi*radius**3
    return(volume)


# define factorial with recursivity ######################

def recursive_factorial(n):
    '''Returns the factorial of the input number'''

    if n == 1:
        return 1
    
    else:
        value = n*recursive_factorial(n-1)
        return(value)
    

# define factorial without recursivity ####################

def factorial(n):
    '''Returns the factorial of the input number'''
    i = 1
    value = n
    while i != n-1:
        value = value*(n-i)
        i = i + 1
    return(value)

# define count_up with recursivity #######################

def recursive_count_up(n , odd = False , counter = 0):
    '''Counts from 0 to n. If odd = True, show only odd numbers'''
    count = counter

    # base case: if counter gets to n
    if count == n:
        if odd == True and n%2 == 0:
            
            return()

        else:
            print(count)
            return(count)

    # if not in base case:
    else:
        if odd == False:
            print(count)          # print every number if odd = False
            odd_value = False
        
        else:
            if (count%2 == 1):    # print odd numbers if odd = True
                print(count)
            odd_value = True
        
        count = counter + 1       # update counter
        recursive_count_up(n, odd = odd_value, counter = count) 


# define count_up without recursivity ###################

def count_up(n,odd = False):
    '''Counts from 0 to n. If odd = True, show only odd numbers'''
    if odd == True:
        for i in range(1 , n+1 , 2):
            print(i)
    
    else:
        for i in range(0,n+1):
            print(i)


# debug code and create function get_final_price 

def get_final_price(price , discount_percentage = 10):
    '''Return the final price after applying the discount percentage'''
    return(price - (price*discount_percentage/100))