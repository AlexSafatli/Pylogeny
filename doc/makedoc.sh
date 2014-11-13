epydoc --pdf --name Pylogeny --url "http://AlexSafatli.github.io/Pylogeny" --graph all *.py -o ../doc
find ../doc/ -type f -not -name api.pdf -not -name makedoc.sh -delete
epydoc --html --name Pylogeny --url "http://AlexSafatli.github.io/Pylogeny" --graph all *.py -o /githubpages/Pylogeny/
rm *.pyc