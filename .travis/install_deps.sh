#see  https://github.com/tseemann/shovill/blob/master/.travis/install_deps.sh

HERE="$PWD"
WGET="wget --quiet"
MAKE="make --silent -j"
UNTAR="tar xf"

HTSLIBVER=1.9
HTSLIB=htslib-$HTSLIBVER
echo "* $HTSLIB"
$WGET https://github.com/samtools/htslib/releases/download/$HTSLIBVER/$HTSLIB.tar.bz2
$UNTAR $HTSLIB.tar.bz2
(cd $HTSLIB && ./configure --prefix=$HERE/$HTSLIB && $MAKE install)
PATH=$HERE/$HTSLIB/bin:$PATH


SAMTOOLSVER=1.9
SAMTOOLS=samtools-$SAMTOOLSVER
echo "* $SAMTOOLS"
$WGET https://github.com/samtools/samtools/releases/download/$SAMTOOLSVER/$SAMTOOLS.tar.bz2
$UNTAR $SAMTOOLS.tar.bz2
(cd $SAMTOOLS && ./configure --prefix=$HERE/$SAMTOOLS && $MAKE install)
PATH=$HERE/$SAMTOOLS/bin:$PATH

BCFTOOLSVER=1.9
BCFTOOLS=bcftools-$BCFTOOLSVER
echo "* $BCFTOOLS"
$WGET https://github.com/samtools/bcftools/releases/download/$BCFTOOLSVER/$BCFTOOLS.tar.bz2
$UNTAR $BCFTOOLS.tar.bz2
(cd $BCFTOOLS && ./configure --prefix=$HERE/$BCFTOOLS && $MAKE install)
PATH=$HERE/$BCFTOOLS/bin:$PATH

PATH=$PATH:$HERE

echo $PATH
export PATH
