from base import *
suite = loader().discover('.',pattern='*Test.py')
tests(verbosity=2).run(suite)