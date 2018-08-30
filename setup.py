from distutils.core import setup

setup(
      name='p-SCNAClonal',
      version='0.0.1',
      description='p-scnaclonal: a improved tool for inferring tumor subclonal populations',
      author='Chu Yanshuo, Nie Chenxi',
      author_email='yanshuochu@hit.edu.cn, 1160300507@stu.hit.edu.cn',
      url='https://github.com/Billy-Nie/p-scnaclonal',
      package_dir={'': 'lib'},
      packages=['pSCNAClonal',
                'pSCNAClonal.preprocess',
                'pSCNAClonal.model',
                'pSCNAClonal.postprocess'],
      scripts=['run.py']
     )
