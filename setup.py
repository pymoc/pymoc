from distutils.core import setup
setup(
    name='pymoc',
    version='1.0',
    package_dir={'pymoc': 'src'},
    py_modules=[
        'pymoc.modules.column',
        'pymoc.modules.equi_column',
        'pymoc.modules.psi_SO',
        'pymoc.modules.psi_thermwind',
        'pymoc.modules.SO_ML',
        'pymoc.utils.check_numpy_version',
        'pymoc.utils.gridit',
        'pymoc.utils.make_array',
        'pymoc.utils.make_func',
    ]
)
