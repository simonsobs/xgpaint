To install:

    rm *.c */*.c *pyc */*.pyc # if you modify cython code
    python setup.py install [--user]
    
To run example:

    cd example
    python makecatalog.py
    python makemaps.py
    python showmap.py
