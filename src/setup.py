from distutils.core import setup
setup(name='pymoc',
      version='1.0',
      package_dir={'pymoc': 'src/modules'},
      py_modules=[
          'pymoc.column', 'pymoc.equi_column', 'pymoc.psi_SO',
          'pymoc.psi_thermwind', 'pymoc.SO_ML'
      ])
