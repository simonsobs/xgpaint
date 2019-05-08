import numpy              as     np
import scipy              as     sp

#--------------------------------------------------
#               Interpolation
#--------------------------------------------------
def tableinterp2D_serial(func,x1,x2,y1,y2,nx,ny):
    '''
    Table for functions of two variables.  
        Since in general functions are not 
    vectorizable, a for loop is required
    '''
    x = np.linspace(x1,x2,nx)
    y = np.linspace(y1,y2,ny)
    f = np.zeros([nx,ny])
    for i in range(nx):
        xval = xt[i]
        for j in range(ny):
            yval = yt[j]
            f[i][j] = func(xval,yval)

    return sp.interpolate.RectBivariateSpline(x,y,f)

def tableinterp1D_serial(func,x1,x2,n):
    '''
    Table for functions of one variable.  
        Since in general functions are not 
    vectorizable, a for loop is required
    '''
    x = np.linspace(x1,x2,n)
    f = np.zeros(nx)
    for i in range(nx):
        xval = x[i]
        f[i] = func(xval)

    return sp.interpolate.interp1d(x,f)

def function2table_2D(func,x1,x2,y1,y2,nx,ny):

    '''
    Parameters
    ------------------------------------------------------
    func: 
        The function to have in the table
    x1:
        Minimum value of table in x-axis
    x2:
        Maximum value of table in x-axis
    y1:
        Minimum value of table in y-axis
    y2:
        Maximum value of table in y-axis
    nx:
        Dimension of table in x-axis
    ny:
        Dimension of table in y-axis

    Returns
    ------------------------------------------------------
    t:
       Table with values of the function in it 	    
    xt/yt:
        Either z or M array depending on return_option
        specified in params file. 

    '''

    t   = np.zeros([nx,ny])
    xt  = np.linspace(x1,x2,nx)
    yt  = np.linspace(y1,y2,ny)
    yt  = yt[::-1]
    x,y = np.meshgrid(xt,yt)

    for i in range(len(x)):
        for j in range(len(y)):
            t[i][j] = func(x[i][j],y[i][j])
    if   params.return_option == "z":
        return t,yt
    elif params.return_option == "M":
        return t,xt
