from phys_const import Gyr_per_Gpc as ggpc

def ttf(x):
    '''Converts number x to label
    if time > 1000 Gyr divide by 1000 and change units to Tyr '''
    unit="Gyr"
    if x > 1000:
        x=x/1000
        unit="Tyr"
    xstr = r'${0:.2f} {1}\,'.format(x,unit)
    return xstr+r'h^{-1}$'
