# Usage: makedoc.sh <path_to_gh-pages>
version=`python -c "from __version__ import VERSION; print VERSION"`
echo "Pylogeny ver. $version documentation being created..."
epydoc --pdf --name Pylogeny --url "http://AlexSafatli.github.io/Pylogeny" --graph all *.py -o ../docs
find ../docs/ -type f -not -name api.pdf -not -name makedoc.sh -delete
cd .. ; sphinx-apidoc -F -e -H Pylogeny -A "Alex Safatli" -V $version -R $version -o docs pylogeny ; cd - ; rm *.pyc
cd ../docs/ ; make html ; cp -r _build/html/* $1 ; find ../docs/ -type f -not -name api.pdf -not -name makedoc.sh -delete ; cd -
cd $1 ; git add * ; git commit -m "Version $1 documentation uploaded."; git push origin gh-pages ; cd -
