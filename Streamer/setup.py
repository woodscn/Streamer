def configuration(parent_package='',top_path=None,compile_type='debug'):
    from numpy.distutils.misc_util import Configuration
    opts={'debug':'--debug --noopt --noarch','run':'--opt=-O3'}
    config = Configuration('Streamer',parent_package,top_path)
    config.add_extension('Godunov','Godunov.f90')
    config.add_extension('BoundaryConditionsStuff','BoundaryConditionsStuff.f90')
    config.add_extension('STLA_IO','stla_io.f90')
    config.add_extension('TimeAdvancementStuff','TimeAdvancementStuff.f90')
    config.add_extension('CGNS_Interface','cgns_interface.f90',include_dirs=
                         '/usr/local/include',libraries=['cgns'],library_dirs=
        '/usr/local/lib')
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
