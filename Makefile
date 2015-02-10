NAME=Pylogeny
PKG=pylogeny
DOCS=docs
FILES=$(shell find ${PKG} -type f -iname "*.py")
VER=$(shell python -c "from ${PKG}.__version__ import VERSION; print VERSION")
DIST=build sdist upload
EDOC=epydoc
SDOC=sphinx-apidoc
PY=python
.PHONY: check_environment

all: ${FILES} setup.py
	${PY} setup.py ${DIST}

docs: ${FILES} check_environment
	${EDOC} --pdf --name ${NAME} --url "http://AlexSafatli.github.io/Pylogeny" --graph all ${FILES} -o ${DOCS}
	${SDOC} -F -e -H ${NAME} -A "Alex Safatli" -V ${VER} -R ${VER} -o ${DOCS} ${PKG}
	cd ${DOCS} && make html
	cp -r ${DOCS}/_build/html/* $(GH_DOC_REPO)
	cd $(GH_DOC_REPO) && git add * && git commit -m "Version ${VER} documentation." && git push origin gh-pages

clean:
	-find ${DOCS} -type f -not -name api.pdf -delete

check_environment:
ifndef GH_DOC_REPO
	$(error GH_DOC_REPO is undefined)
endif
