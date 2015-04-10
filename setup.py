from setuptools import setup
files = ["template/*",
         "halley/*"]

setup(name='cpolymer',
      version='0.3',
      description='Creating initial configuration for polymer (to work with lammps)',
      url='https://github.com/jeammimi/cpolymer',
      author='Jean-michel Arbona',
      author_email='jeanmichel.arbona@gmail.com',
      license='MIT',
      packages=['cpolymer'],
      zip_safe=False)
