# Usage: makedoc.sh <version> <path_to_gh-pages>
epydoc --pdf --name Pylogeny --url "http://AlexSafatli.github.io/Pylogeny" --graph all *.py -o ../docs
find ../docs/ -type f -not -name api.pdf -not -name makedoc.sh -delete
cd .. ; sphinx-apidoc -F -e -H Pylogeny -A "Alex Safatli" -V $1 -R $1 -o docs pylogeny ; cd - ; rm *.pyc
cd ../docs/ ; make html ; cp -r _build/html/* $2 ; find ../docs/ -type f -not -name api.pdf -not -name makedoc.sh -delete ; cd -
cd $2 ; git add * ; git commit -m "Version $1 documentation uploaded."; git push origin gh-pages ; cd -
