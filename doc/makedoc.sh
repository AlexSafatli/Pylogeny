epydoc --pdf --name Pylogeny --url "http://www.alexsafatli.me" --graph all *.py -o ../doc
find ../doc/ -type f -not -name api.pdf -not -name makedoc.sh -delete
