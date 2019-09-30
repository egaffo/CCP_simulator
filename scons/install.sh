## ciri simulator
wget https://sourceforge.net/projects/ciri/files/CIRI_simulator.pl -O bin/CIRI_simulator.pl

## gffread from Cufflinks
wget http://cole-trapnell-lab.github.io/cufflinks/assets/downloads/cufflinks-2.2.1.Linux_x86_64.tar.gz -O bin/cufflinks-2.2.1.Linux_x86_64.tar.gz
tar xf bin/cufflinks-2.2.1.Linux_x86_64.tar.gz
ln -s cufflinks-2.2.1.Linux_x86_64/gffread
