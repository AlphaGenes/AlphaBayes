mkdir AlphaBayes

# Assumes that the program and manual have both been built.

# To build the program run:
# NOTE: Binaries should be moved to the "binaries" folder and uploaded to bitbucket after builds.
#cmake . ; make

# to build the manual using Sphinx:
#( cd doc ; make latexpdf )

cp -r example AlphaBayes

# Copy in the documentation.
cp doc/build/latex/AlphaBayes.pdf AlphaBayes/AlphaBayes.pdf

if [ $? != 0 ]; then                   # last command: echo
    echo "The manual needs to be built." # last command: [
    exit 1
fi


# Copy in the binaries
cp binaries/* AlphaBayes

# Create a version file

version=`git describe --tags --abbrev=0`
commit=`git rev-parse --short HEAD`

echo Version: $version > AlphaBayes/version.txt
echo Commit: $commit >> AlphaBayes/version.txt

zip -r AlphaBayes.zip AlphaBayes
