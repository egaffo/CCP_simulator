import os, sys

def SymLink(target, source, env):
    os.symlink(os.path.abspath(str(source[0])), os.path.abspath(str(target[0])))

env = Environment(ENV=os.environ, SHELL = '/bin/bash')
ccp_bin_dir = os.path.join(env['ENV']['CCPSIM_HOME'], 'bin')

## ciri simulator
#wget https://sourceforge.net/projects/ciri/files/CIRI_simulator.pl -O bin/CIRI_simulator.pl
# CIRI SIMULATOR
CIRISIM_script = 'CIRI_simulator.pl'
CIRISIM_link = 'https://sourceforge.net/projects/ciri/files/' + CIRISIM_script
CIRISIM_target = [os.path.join(ccp_bin_dir, CIRISIM_script)]
CIRISIM = env.Command(CIRISIM_target, [], ['wget -O ${TARGETS[0]} ' + CIRISIM_link])

# ## gffread from Cufflinks
# wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz -O bin/cufflinks-2.2.1.Linux_x86_64.tar.gz
# tar xf bin/cufflinks-2.2.1.Linux_x86_64.tar.gz
# ln -s cufflinks-2.2.1.Linux_x86_64/gffread

## Polyester

env['ENV']['R_LIBS'] = os.path.join(ccp_bin_dir, "R_libs")

polyester_cmd = '''R -e 'if (!requireNamespace("BiocManager", quietly=TRUE)){install.packages("BiocManager")}; BiocManager::install("polyester")' '''
polyester_targets = [os.path.join(ccp_bin_dir, 'R_libs', 'polyester', 'R', 'polyester')]
polyester = env.Command(polyester_targets, 
                        [], 
                        polyester_cmd)

optparse_cmd = '''R -e 'install.packages("optparse")' '''
optparse_targets = [os.path.join(ccp_bin_dir, 'R_libs', 'optparse', 'R', 'optparse')]
optparse = env.Command(optparse_targets, 
                       [], 
                       optparse_cmd)

