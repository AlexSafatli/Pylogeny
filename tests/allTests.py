from base import *
suite = loader().discover('tests',pattern='*Test.py')
tests(verbosity=2).run(suite)
