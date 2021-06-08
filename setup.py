import os
try:
    from setuptools import setup
    from setuptools.command.install import install
except ImportError:
    from distutils.core import setup
    from distutils.command.install import install



class CheckPresenceOfExecutables(install):
    """
    Customized setuptools install command - checks for presence of
    necessary executables in path.
    """

    needed_executables = {
        'mafft': 'https://mafft.cbrc.jp/alignment/software/'
    }

    class ExecutableException(Exception):
        pass

    def which(self, program):
        """
        Check if an executable exists.

        Args:
            program: An executable name or path.

        Returns:
            The full path to the program on success, else None.
        """

        def is_exe(fpath):
            return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

        fpath, fname = os.path.split(program)
        if fpath:
            if is_exe(program):
                return program
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                path = path.strip('"')
                exe_file = os.path.join(path, program)
                if is_exe(exe_file):
                    return exe_file

        return None


    def run(self):
        for executable in self.needed_executables.keys():
            if self.which(executable) is None:
                raise self.ExecutableException(
                    "Executable " + executable + " is not present. " \
                    "Please download and install it from "
                    + self.needed_executables[executable]
                )
        install.run(self)

setup(
    name='intactness-pipeline',
    description='Test for proviral intactness.',
    author='Imogen Wright',
    author_email='imogen@hyraxbio.co.za',
    version='0.6',
    packages=['intact', 'util'],
    package_data={
                  'intact': ['data/*'],
                  'util': ['subtype_alignments/*']
                  },
    scripts=[
        'bin/proviral'
    ],
    cmdclass={
        'install': CheckPresenceOfExecutables
    },
    install_requires=[
        'appdirs>=1.4.3',
        'biopython>=1.71',
        'click>=6.7',
        'scipy>=1.6.0',
        'numpy>1.19.5'
    ],
    setup_requires=[
        'appdirs>=1.4.3',
        'biopython>=1.71',
        'click>=6.7',
        'scipy>=1.6.0',
        'numpy>=1.19.5'
    ]
)

